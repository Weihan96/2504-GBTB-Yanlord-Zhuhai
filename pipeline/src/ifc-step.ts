import { createHash } from "node:crypto";

export interface IfcElementRecord {
  stepId: number;
  globalId: string | null;
  name: string | null;
  objectType: string | null;
  tag?: string | null;
  longName?: string | null;
  predefinedType?: string | null;
  overallHeight?: number | null;
  overallWidth?: number | null;
  fillsOpening?: boolean;
  reference?: string | null;
}

export interface IfcSnapshot {
  source: {
    path: string;
    byteSize: number;
    sha256: string;
    headerTimestamp: string | null;
  };
  schema: string | null;
  projectName: string | null;
  entityCounts: Record<string, number>;
  integrity: {
    stepEntityCount: number;
    duplicateStepIds: number[];
    guidLikeCount: number;
    duplicateGlobalIds: string[];
    hasValidFooter: boolean;
  };
  annotations: {
    total: number;
    byObjectType: Record<string, number>;
    drawings: Array<{ globalId: string | null; name: string | null }>;
    elevations: Array<{ globalId: string | null; name: string | null }>;
  };
  documents: Array<{ location: string | null; identification: string | null }>;
  grids: { axes: string[] };
  spaces: {
    total: number;
    missingName: number;
    missingLongName: number;
    provisional: number;
    pendingDelete: number;
    containedInStorey: number;
    aggregatedUnderStorey: number;
    records: IfcElementRecord[];
  };
  doors: { total: number; missingTag: number; missingOverallSize: number; missingOpeningFill: number; records: IfcElementRecord[] };
  windows: { total: number; missingTag: number; missingOverallSize: number; missingOpeningFill: number; records: IfcElementRecord[] };
}

interface ParsedEntity {
  id: number;
  type: string;
  args: string[];
}

export function splitStepArgs(source: string): string[] {
  const result: string[] = [];
  let start = 0;
  let depth = 0;
  let quoted = false;

  for (let index = 0; index < source.length; index += 1) {
    const character = source[index];
    if (character === "'") {
      if (quoted && source[index + 1] === "'") {
        index += 1;
        continue;
      }
      quoted = !quoted;
      continue;
    }
    if (quoted) continue;
    if (character === "(") depth += 1;
    if (character === ")") depth -= 1;
    if (character === "," && depth === 0) {
      result.push(source.slice(start, index).trim());
      start = index + 1;
    }
  }
  result.push(source.slice(start).trim());
  return result;
}

export function decodeStepString(value: string): string | null {
  const trimmed = value.trim();
  if (trimmed === "$" || trimmed === "*") return null;
  if (!trimmed.startsWith("'") || !trimmed.endsWith("'")) return trimmed;
  let content = trimmed.slice(1, -1).replaceAll("''", "'");
  content = content.replace(/\\X2\\([0-9A-F]+)\\X0\\/gi, (_, hex: string) => {
    let decoded = "";
    for (let index = 0; index < hex.length; index += 4) {
      decoded += String.fromCharCode(Number.parseInt(hex.slice(index, index + 4), 16));
    }
    return decoded;
  });
  content = content.replace(/\\X4\\([0-9A-F]+)\\X0\\/gi, (_, hex: string) => {
    let decoded = "";
    for (let index = 0; index < hex.length; index += 8) {
      decoded += String.fromCodePoint(Number.parseInt(hex.slice(index, index + 8), 16));
    }
    return decoded;
  });
  return content;
}

function stepEnum(value: string | undefined): string | null {
  if (!value || value === "$" || value === "*") return null;
  const match = /^\.([^.]*)\.$/.exec(value.trim());
  return match ? match[1] : decodeStepString(value);
}

function stepNumber(value: string | undefined): number | null {
  if (!value || value === "$" || value === "*") return null;
  const parsed = Number(value);
  return Number.isFinite(parsed) ? parsed : null;
}

function stepTypedValue(value: string | undefined): string | null {
  if (!value || value === "$" || value === "*") return null;
  const match = /^[A-Z0-9_]+\((.*)\)$/.exec(value.trim());
  return decodeStepString(match?.[1] ?? value);
}

function parseEntityLine(line: string): ParsedEntity | null {
  const match = /^#(\d+)=([A-Z0-9_]+)\((.*)\);$/.exec(line);
  if (!match) return null;
  return {
    id: Number.parseInt(match[1], 10),
    type: match[2],
    args: splitStepArgs(match[3]),
  };
}

function extractReferences(value: string | undefined): number[] {
  if (!value) return [];
  return [...value.matchAll(/#(\d+)/g)].map((match) => Number.parseInt(match[1], 10));
}

function isCompressedGuid(value: string | null): value is string {
  return Boolean(value && /^[0-3][0-9A-Za-z_$]{21}$/.test(value));
}

function increment(record: Record<string, number>, key: string): void {
  record[key] = (record[key] ?? 0) + 1;
}

function createRecord(entity: ParsedEntity): IfcElementRecord {
  return {
    stepId: entity.id,
    globalId: decodeStepString(entity.args[0]),
    name: decodeStepString(entity.args[2]),
    objectType: decodeStepString(entity.args[4]),
  };
}

export function analyzeIfcText(text: string, path = "model.ifc", byteSize = Buffer.byteLength(text), sha256?: string): IfcSnapshot {
  const counts: Record<string, number> = {};
  const stepIds = new Set<number>();
  const duplicateStepIds = new Set<number>();
  const globalIds = new Set<string>();
  const duplicateGlobalIds = new Set<string>();
  const annotationTypes: Record<string, number> = {};
  const drawings: Array<{ globalId: string | null; name: string | null }> = [];
  const elevations: Array<{ globalId: string | null; name: string | null }> = [];
  const documents: Array<{ location: string | null; identification: string | null }> = [];
  const axes: string[] = [];
  const spaces: IfcElementRecord[] = [];
  const doors: IfcElementRecord[] = [];
  const windows: IfcElementRecord[] = [];
  const fillElementIds = new Set<number>();
  const containedRelations: number[][] = [];
  const aggregateRelations: number[][] = [];
  const propertySingleValues = new Map<number, { name: string | null; value: string | null }>();
  const propertySets = new Map<number, { name: string | null; propertyIds: number[] }>();
  const propertyRelations: Array<{ relatedIds: number[]; propertySetId: number | null }> = [];
  let projectName: string | null = null;

  let start = 0;
  for (let index = 0; index <= text.length; index += 1) {
    if (index !== text.length && text.charCodeAt(index) !== 10) continue;
    const end = index > start && text.charCodeAt(index - 1) === 13 ? index - 1 : index;
    const line = text.slice(start, end);
    start = index + 1;
    const entity = parseEntityLine(line);
    if (!entity) continue;

    increment(counts, entity.type);
    if (stepIds.has(entity.id)) duplicateStepIds.add(entity.id);
    stepIds.add(entity.id);

    const possibleGuid = decodeStepString(entity.args[0]);
    if (isCompressedGuid(possibleGuid)) {
      if (globalIds.has(possibleGuid)) duplicateGlobalIds.add(possibleGuid);
      globalIds.add(possibleGuid);
    }

    if (entity.type === "IFCPROJECT") projectName = decodeStepString(entity.args[2]);
    if (entity.type === "IFCANNOTATION") {
      const objectType = decodeStepString(entity.args[4]) ?? "<missing>";
      increment(annotationTypes, objectType);
      if (objectType === "DRAWING") drawings.push({ globalId: possibleGuid, name: decodeStepString(entity.args[2]) });
      if (objectType === "ELEVATION") elevations.push({ globalId: possibleGuid, name: decodeStepString(entity.args[2]) });
    }
    if (entity.type === "IFCPROPERTYSINGLEVALUE") {
      propertySingleValues.set(entity.id, {
        name: decodeStepString(entity.args[0]),
        value: stepTypedValue(entity.args[2]),
      });
    }
    if (entity.type === "IFCPROPERTYSET") {
      propertySets.set(entity.id, {
        name: decodeStepString(entity.args[2]),
        propertyIds: extractReferences(entity.args[4]),
      });
    }
    if (entity.type === "IFCRELDEFINESBYPROPERTIES") {
      propertyRelations.push({
        relatedIds: extractReferences(entity.args[4]),
        propertySetId: extractReferences(entity.args[5])[0] ?? null,
      });
    }
    if (entity.type === "IFCDOCUMENTREFERENCE") {
      documents.push({ identification: decodeStepString(entity.args[0]), location: decodeStepString(entity.args[1]) });
    }
    if (entity.type === "IFCGRIDAXIS") {
      const axis = decodeStepString(entity.args[0]);
      if (axis) axes.push(axis);
    }
    if (entity.type === "IFCSPACE") {
      spaces.push({
        ...createRecord(entity),
        longName: decodeStepString(entity.args[7]),
        predefinedType: stepEnum(entity.args[9]),
      });
    }
    if (entity.type === "IFCDOOR" || entity.type === "IFCWINDOW") {
      const record: IfcElementRecord = {
        ...createRecord(entity),
        tag: decodeStepString(entity.args[7]),
        overallHeight: stepNumber(entity.args[8]),
        overallWidth: stepNumber(entity.args[9]),
        predefinedType: stepEnum(entity.args[10]),
      };
      if (entity.type === "IFCDOOR") doors.push(record);
      else windows.push(record);
    }
    if (entity.type === "IFCRELFILLSELEMENT") {
      const relatedElement = extractReferences(entity.args[5])[0];
      if (relatedElement) fillElementIds.add(relatedElement);
    }
    if (entity.type === "IFCRELCONTAINEDINSPATIALSTRUCTURE") containedRelations.push(extractReferences(entity.args[4]));
    if (entity.type === "IFCRELAGGREGATES") aggregateRelations.push(extractReferences(entity.args[5]));
  }

  for (const record of [...doors, ...windows]) record.fillsOpening = fillElementIds.has(record.stepId);
  const spaceByStepId = new Map(spaces.map((space) => [space.stepId, space]));
  for (const relation of propertyRelations) {
    const propertySet = relation.propertySetId ? propertySets.get(relation.propertySetId) : undefined;
    if (propertySet?.name !== "Pset_SpaceCommon") continue;
    const reference = propertySet.propertyIds
      .map((propertyId) => propertySingleValues.get(propertyId))
      .find((property) => property?.name === "Reference")?.value;
    if (!reference) continue;
    for (const relatedId of relation.relatedIds) {
      const space = spaceByStepId.get(relatedId);
      if (space) space.reference = reference;
    }
  }
  const spaceIds = new Set(spaces.map((space) => space.stepId));
  const countRelatedSpaces = (relations: number[][]) => {
    const related = new Set<number>();
    for (const relation of relations) {
      for (const reference of relation) if (spaceIds.has(reference)) related.add(reference);
    }
    return related.size;
  };

  const schemaMatch = /FILE_SCHEMA\(\('([^']+)'\)\);/.exec(text);
  const timestampMatch = /FILE_NAME\('[^']*','([^']+)'/.exec(text);
  const byName = (left: IfcElementRecord, right: IfcElementRecord) => (left.name ?? "").localeCompare(right.name ?? "");
  spaces.sort(byName);
  doors.sort(byName);
  windows.sort(byName);

  return {
    source: {
      path,
      byteSize,
      sha256: sha256 ?? createHash("sha256").update(text).digest("hex"),
      headerTimestamp: timestampMatch?.[1] ?? null,
    },
    schema: schemaMatch?.[1] ?? null,
    projectName,
    entityCounts: Object.fromEntries(Object.entries(counts).sort(([left], [right]) => left.localeCompare(right))),
    integrity: {
      stepEntityCount: stepIds.size,
      duplicateStepIds: [...duplicateStepIds].sort((left, right) => left - right),
      guidLikeCount: globalIds.size,
      duplicateGlobalIds: [...duplicateGlobalIds].sort(),
      hasValidFooter: /ENDSEC;\s*END-ISO-10303-21;\s*$/.test(text),
    },
    annotations: { total: counts.IFCANNOTATION ?? 0, byObjectType: annotationTypes, drawings, elevations },
    documents,
    grids: { axes: [...new Set(axes)].sort() },
    spaces: {
      total: spaces.length,
      missingName: spaces.filter((space) => !space.name).length,
      missingLongName: spaces.filter((space) => !space.longName).length,
      provisional: spaces.filter((space) => space.objectType === "PROVISIONAL_GRID_CELL").length,
      pendingDelete: spaces.filter((space) => space.longName?.includes("待删除")).length,
      containedInStorey: countRelatedSpaces(containedRelations),
      aggregatedUnderStorey: countRelatedSpaces(aggregateRelations),
      records: spaces,
    },
    doors: {
      total: doors.length,
      missingTag: doors.filter((door) => !door.tag).length,
      missingOverallSize: doors.filter((door) => !door.overallHeight || !door.overallWidth).length,
      missingOpeningFill: doors.filter((door) => !door.fillsOpening).length,
      records: doors,
    },
    windows: {
      total: windows.length,
      missingTag: windows.filter((window) => !window.tag).length,
      missingOverallSize: windows.filter((window) => !window.overallHeight || !window.overallWidth).length,
      missingOpeningFill: windows.filter((window) => !window.fillsOpening).length,
      records: windows,
    },
  };
}

export async function snapshotIfc(path: string): Promise<IfcSnapshot> {
  const buffer = Buffer.from(await Bun.file(path).arrayBuffer());
  const sha256 = createHash("sha256").update(buffer).digest("hex");
  return analyzeIfcText(buffer.toString("utf8"), path, buffer.byteLength, sha256);
}
