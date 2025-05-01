#!/bin/bash

# Usage check
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <source_mesh_output_path> <target_mesh_output_path>"
    echo "Example: $0 ./source_mesh ./target_mesh"
    exit 1
fi

# Get output paths from arguments
SOURCE_PATH="$1"
TARGET_PATH="$2"

# Create the parent directories if they don't exist
mkdir -p "$(dirname "$SOURCE_PATH")"
mkdir -p "$(dirname "$TARGET_PATH")"

# Set up environment variables for Gmsh and Omega_h
export LD_LIBRARY_PATH=/lore/paudea/Gmsh/gmsh_install/lib64:/lore/paudea/build/ADA89-meshField/omega_h-meshField/install/lib64:$LD_LIBRARY_PATH
/lore/paudea/omega_h/install/lib:$LD_LIBRARY_PATH
# Path to executables
GMSH_PATH="/lore/paudea/Gmsh/gmsh_install/bin/gmsh"
OMEGA_H_CONVERT="/lore/paudea/omega_h/install/bin/osh_convert"

# Check if Gmsh executable exists
if [ ! -f "$GMSH_PATH" ]; then
    echo "Error: Gmsh executable not found at $GMSH_PATH"
    echo "Please update the GMSH_PATH variable in the script."
    exit 1
fi

# Create a temporary directory for the .geo files
TEMP_DIR=$(mktemp -d)
trap 'rm -rf "$TEMP_DIR"' EXIT

# Create source mesh .geo file - fine mesh
cat > "$TEMP_DIR/source_mesh.geo" << EOF
// Source mesh - finer mesh with smaller elements
L = 1.0;
cl = 0.03;  // Smaller characteristic length for source mesh

Point(1) = {0, 0, 0, cl};
Point(2) = {L, 0, 0, cl};
Point(3) = {L, L, 0, cl};
Point(4) = {0, L, 0, cl};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Mesh.Algorithm = 8;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.MeshSizeMin = cl;
Mesh.MeshSizeMax = cl;
Mesh.ElementOrder = 1;
EOF

# Create target mesh .geo file - coarser mesh
cat > "$TEMP_DIR/target_mesh.geo" << EOF
// Target mesh - coarser mesh with larger elements
L = 1.0;
cl = 0.06;  // Larger characteristic length for target mesh

Point(1) = {0, 0, 0, cl};
Point(2) = {L, 0, 0, cl};
Point(3) = {L, L, 0, cl};
Point(4) = {0, L, 0, cl};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Mesh.Algorithm = 8;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.MeshSizeMin = cl;
Mesh.MeshSizeMax = cl;
Mesh.ElementOrder = 1;
EOF

echo "============================================"
echo "Generating source mesh..."
echo "============================================"

# Generate source mesh
$GMSH_PATH "$TEMP_DIR/source_mesh.geo" -2 -o "$TEMP_DIR/source_mesh.msh" -format msh2
if [ $? -ne 0 ]; then
    echo "Error: Failed to generate source mesh"
    exit 1
fi

# Convert to Omega_h format if converter exists
if [ -f "$OMEGA_H_CONVERT" ]; then
    $OMEGA_H_CONVERT "$TEMP_DIR/source_mesh.msh" "$SOURCE_PATH.osh"
    if [ $? -eq 0 ]; then
        echo "Successfully created source mesh: $SOURCE_PATH.osh"
    else
        echo "Error: Failed to convert source mesh to Omega_h format"
        exit 1
    fi
else
    # Just copy the .msh file if converter not available
    cp "$TEMP_DIR/source_mesh.msh" "$SOURCE_PATH.msh"
    echo "Source mesh created: $SOURCE_PATH.msh"
    echo "Note: Omega_h converter not found, using Gmsh format directly"
fi

echo "============================================"
echo "Generating target mesh..."
echo "============================================"

# Generate target mesh
$GMSH_PATH "$TEMP_DIR/target_mesh.geo" -2 -o "$TEMP_DIR/target_mesh.msh" -format msh2
if [ $? -ne 0 ]; then
    echo "Error: Failed to generate target mesh"
    exit 1
fi

# Convert to Omega_h format if converter exists
if [ -f "$OMEGA_H_CONVERT" ]; then
    $OMEGA_H_CONVERT "$TEMP_DIR/target_mesh.msh" "$TARGET_PATH.osh"
    if [ $? -eq 0 ]; then
        echo "Successfully created target mesh: $TARGET_PATH.osh"
    else
        echo "Error: Failed to convert target mesh to Omega_h format"
        exit 1
    fi
else
    # Just copy the .msh file if converter not available
    cp "$TEMP_DIR/target_mesh.msh" "$TARGET_PATH.msh"
    echo "Target mesh created: $TARGET_PATH.msh"
    echo "Note: Omega_h converter not found, using Gmsh format directly"
fi

echo "============================================"
echo "Mesh generation complete!"
echo "============================================"
echo "Source mesh: $SOURCE_PATH.osh (or .msh)"
echo "Target mesh: $TARGET_PATH.osh (or .msh)"
echo ""
echo "To use these meshes in your particle2mesh_map code:"
echo "  Omega_h::Mesh source_mesh(&comm);"
echo "  Omega_h::Mesh target_mesh(&comm);"
echo "  Omega_h::meshfile::read(\"$SOURCE_PATH.osh\", source_mesh);"
echo "  Omega_h::meshfile::read(\"$TARGET_PATH.osh\", target_mesh);"
