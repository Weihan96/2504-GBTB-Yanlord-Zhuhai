"""Cache a full baseline schema scan keyed to unchanged formal IFC bytes."""
import json,re,ifcopenshell,ifcopenshell.validate
from collections import Counter
import migrate_tab02 as task

def errors(model):
    logger=ifcopenshell.validate.json_logger();ifcopenshell.validate.validate(model,logger)
    values=Counter((x.get('attribute',''),re.sub(r'#\d+','#STEP',re.sub(r'\s+',' ',x['message'])),task.pkg.fingerprint(x['instance'])) for x in logger.statements)
    return [{'key':list(k),'count':v} for k,v in sorted(values.items())]

if __name__=='__main__':
    data={'formal_sha256':task.pkg.sha256(task.FORMAL),'errors':errors(ifcopenshell.open(str(task.FORMAL)))}
    task.write(task.OUT/'formal-schema-baseline.json',data)
    print({'baseline_errors':sum(x['count'] for x in data['errors'])})
