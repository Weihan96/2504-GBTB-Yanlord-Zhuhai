#!/usr/bin/env bash
set -e

echo "=== Blender Shapely Installer (macOS) ==="

# ----------------------------
# 1. Find Blender.app
# ----------------------------
BLENDER_APP=$(find /Applications -maxdepth 1 -name "Blender*.app" | head -n 1)

if [ -z "$BLENDER_APP" ]; then
    echo "❌ Blender.app not found in /Applications"
    echo "Please move Blender into /Applications"
    exit 1
fi

echo "Blender found:"
echo "  $BLENDER_APP"

# ----------------------------
# 2. Locate blender python
# ----------------------------
PYTHON_BIN=$(find "$BLENDER_APP/Contents/Resources" \
    -path "*/python/bin/python3*" \
    -type f | head -n 1)

if [ ! -f "$PYTHON_BIN" ]; then
    echo "❌ Blender python not found"
    exit 1
fi

echo "Blender Python:"
echo "  $PYTHON_BIN"

# ----------------------------
# 3. Enable pip
# ----------------------------
echo "Setting up pip..."

"$PYTHON_BIN" -m ensurepip --upgrade || true
"$PYTHON_BIN" -m pip install --upgrade pip

# ----------------------------
# 4. Install shapely
# ----------------------------
echo "Installing shapely..."

"$PYTHON_BIN" -m pip install --upgrade shapely

# ----------------------------
# 5. Test
# ----------------------------
echo "Testing shapely..."

"$PYTHON_BIN" - <<EOF
import shapely
print("✅ Shapely version:", shapely.__version__)
EOF

echo ""
echo "🎉 Done! Blender can now import shapely."
