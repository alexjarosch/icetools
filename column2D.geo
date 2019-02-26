// This .geo file defines a simple column geometry for the icetools_2d_demo.py 
// file. You can mesh the geometry by simply running the command below in your 
// shell:
// ~$ gmsh -2 column2D.geo
// which leaves you with a file called column2D.msh. This can be translated into
// the dolfin mesh format with:
// ~$ dolfin-convert column2D.msh column2D.xml
// The file column2D.xml.bak is not needed, so you can remove it. It is custom 
// within icetools and dolfin to gzip the mesh to save disk space, so you can
// just run
// ~$ gzip column2D.xml
// to produce the final mesh

// corresponds to cs_length in icetools_2d_demo.py
cs_length = 40.0;
// corresponds to cs_height in icetools_2d_demo.py
cs_height = 120.0;
// cl is the mesh charateristic length at the mesh point
cl = 10.0;

Point(1) = {0, 0, 0, cl};
Point(2) = {cs_length, 0, 0, cl};
Point(3) = {cs_length, cs_height, 0, cl};
Point(4) = {0, cs_height, 0, cl};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {2, 3, 4, 1};
Plane Surface(6) = {5};
