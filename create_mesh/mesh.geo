// Simple 2D mesh with 6 triangles

// Define points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};
Point(5) = {0.5, 0, 0, 1.0};
Point(6) = {1, 0.5, 0, 1.0};
Point(7) = {0.5, 1, 0, 1.0};
Point(8) = {0, 0.5, 0, 1.0};
Point(9) = {0.5, 0.5, 0, 1.0};

// Define lines
Line(1) = {1, 5};
Line(2) = {5, 2};
Line(3) = {2, 6};
Line(4) = {6, 3};
Line(5) = {3, 7};
Line(6) = {7, 4};
Line(7) = {4, 8};
Line(8) = {8, 1};
Line(9) = {5, 9};
Line(10) = {6, 9};
Line(11) = {7, 9};
Line(12) = {8, 9};

// Define line loop and surface
Line Loop(1) = {1, 9, -12, 8};
Line Loop(2) = {2, 3, 10, -9};
Line Loop(3) = {4, 5, 11, -10};
Line Loop(4) = {6, 7, 12, -11};

// Define surfaces
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

// Now generate 2 more triangles manually (optional)
Line Loop(5) = {9, 10, -11, -12};
Plane Surface(5) = {5};

// Mesh settings (optional)
Mesh 2;
