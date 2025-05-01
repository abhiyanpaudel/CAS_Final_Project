#!/bin/bash

# Directory paths
SOURCE_MESH_DIR="/lore/paudea/projects/particle2mesh_map/create_mesh/source_mesh"
TARGET_MESH_DIR="/lore/paudea/projects/particle2mesh_map/create_mesh/target_mesh"

# Set up environment variables for Gmsh and Omega_h
export LD_LIBRARY_PATH=/lore/paudea/Gmsh/gmsh_install/lib64:/lore/paudea/build/ADA89-meshField/omega_h-meshField/install/lib64:$LD_LIBRARY_PATH

# Path to executables
GMSH_PATH="/lore/paudea/Gmsh/gmsh_install/bin/gmsh"

OMEGA_H_CONVERT="/lore/paudea/build/ADA89-meshField/omega_h-meshField/install/bin/msh2osh"

# Check if Gmsh executable exists
if [ ! -f "$GMSH_PATH" ]; then
    echo "Error: Gmsh executable not found at $GMSH_PATH"
    echo "Please update the GMSH_PATH variable in the script."
    exit 1
fi

# Check if Omega_h converter exists
if [ ! -f "$OMEGA_H_CONVERT" ]; then
    echo "Error: osh_convert tool not found at $OMEGA_H_CONVERT"
    echo "Please update the OMEGA_H_CONVERT variable in the script."
    exit 1
fi

# Function to process .geo files in a directory
process_geo_files() {
    local directory=$1
    local dir_name=$(basename "$directory")
    
    echo "============================================"
    echo "Processing .geo files in $dir_name..."
    echo "============================================"
    
    # Find all .geo files in the directory
    for geo_file in "$directory"/*.geo; do
        if [ ! -f "$geo_file" ]; then
            echo "No .geo files found in $directory"
            return
        fi
        
        base_name=$(basename "$geo_file" .geo)
        msh_file="${directory}/${base_name}.msh"
        osh_file="${directory}/${base_name}.osh"
        
        echo "Processing $geo_file..."
        
        # Generate mesh with Gmsh
        $GMSH_PATH "$geo_file" -2 -o "$msh_file" -format msh2
        if [ $? -ne 0 ]; then
            echo "Error: Failed to generate mesh for $geo_file"
            continue
        fi
        
        # Convert to Omega_h format
        $OMEGA_H_CONVERT "$msh_file" "$osh_file"
        if [ $? -eq 0 ]; then
            echo "Successfully created Omega_h mesh: $osh_file"
        else
            echo "Error: Failed to convert mesh to Omega_h format"
            continue
        fi
    done
}

# Process source mesh directory
if [ -d "$SOURCE_MESH_DIR" ]; then
    process_geo_files "$SOURCE_MESH_DIR"
else
    echo "Source mesh directory not found: $SOURCE_MESH_DIR"
fi

# Process target mesh directory
if [ -d "$TARGET_MESH_DIR" ]; then
    process_geo_files "$TARGET_MESH_DIR"
else
    echo "Target mesh directory not found: $TARGET_MESH_DIR"
fi

echo "============================================"
echo "Mesh processing complete!"
echo "============================================"
echo "To use these meshes in your code:"
echo "  Omega_h::Mesh mesh(&comm);"
echo "  Omega_h::binary::read(\"path/to/mesh.osh\", world, &mesh);"
