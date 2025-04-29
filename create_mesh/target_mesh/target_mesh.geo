// square.geo - A 2D square mesh with uniform triangular elements

// Define the size of the square (side length)
L = 1.0;

// Define the characteristic length (controls element size)
cl = 0.15;  // Adjust this value to change element size

// Define the points
Point(1) = {0, 0, 0, cl};
Point(2) = {L, 0, 0, cl};
Point(3) = {L, L, 0, cl};
Point(4) = {0, L, 0, cl};

// Define the lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Define the surface loop
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Set mesh algorithm to ensure more uniform triangles
Mesh.Algorithm = 8;  // Delaunay triangulation

// Enable mesh optimization to improve element quality
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

// For even more uniform elements, you can set:
Mesh.MeshSizeMin = cl;
Mesh.MeshSizeMax = cl;

// Set element order (1 = linear elements)
Mesh.ElementOrder = 1;

// Set 2D meshing
Mesh 2;

// Save the mesh
Save "source.msh";
