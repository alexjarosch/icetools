// This .geo file defines a simple column geometry for the icetools_3d_demo.py 
// file. You can mesh the geometry by simply running the command below in your 
// shell:
// ~$ gmsh -3 column3D.geo
// which leaves you with a file called column3D.msh. This can be translated into
// the dolfin mesh format with:
// ~$ dolfin-convert column3D.msh column3D.xml
// The file column3D.xml.bak is not needed, so you can remove it. It is custom 
// within icetools and dolfin to gzip the mesh to save disk space, so you can
// just run
// ~$ gzip column3D.xml
// to produce the final mesh

// corresponds to cs_length in icetools_3d_demo.py
cs_length = 40.0;
// corresponds to cs_height in icetools_3d_demo.py
cs_height = 120.0;
// cl is the mesh charateristic length at the mesh point
cl = 10.0;

General.Terminal = 1;

P1 = newp; Point(P1) = {0, 0, 0, cl};
P2 = newp; Point(P2) = {cs_length, 0, 0, cl};
P3 = newp; Point(P3) = {cs_length, 0, cs_height, cl};
P4 = newp; Point(P4) = {0, 0, cs_height, cl};
P5 = newp; Point(P5) = {0, cs_length, 0, cl};
P6 = newp; Point(P6) = {cs_length, cs_length, 0, cl};
P7 = newp; Point(P7) = {cs_length, cs_length, cs_height, cl};
P8 = newp; Point(P8) = {0, cs_length, cs_height, cl};

L01 = newl; Line(L01) = {P1, P2};
L02 = newl; Line(L02) = {P2, P3};
L03 = newl; Line(L03) = {P3, P4};
L04 = newl; Line(L04) = {P4, P1};
L05 = newl; Line(L05) = {P2, P6};
L06 = newl; Line(L06) = {P6, P7};
L07 = newl; Line(L07) = {P7, P3};
L08 = newl; Line(L08) = {P6, P5};
L09 = newl; Line(L09) = {P5, P8};
L10 = newl; Line(L10) = {P8, P7};
L11 = newl; Line(L11) = {P5, P1};
L12 = newl; Line(L12) = {P4, P8};

LL1 = newll; Line Loop(LL1) = { L01,  L02,  L03,  L04};
LL2 = newll; Line Loop(LL2) = { L05,  L06,  L07, -L02};
LL3 = newll; Line Loop(LL3) = { L08,  L09,  L10, -L06};
LL4 = newll; Line Loop(LL4) = { L11, -L04,  L12, -L09};
LL5 = newll; Line Loop(LL5) = { L01,  L05,  L08,  L11};
LL6 = newll; Line Loop(LL6) = {-L03, -L07, -L10, -L12};

// The more complex code below is to generate a mesh usable for periodic boundary
// conditions in Fenics
// the code is adapted from http://www.mail-archive.com/gmsh@geuz.org/msg01649.html
// thanks!

Geometry.AutoCoherence = 0;

S1 = news; Plane Surface(S1) = {LL1};
out[] = Extrude {0., cs_length, 0.} { Surface{S1}; Layers{{1}, {1}}; };
S3 = out[0];
Delete {Volume {out[1]};}
Delete {Surface {out[2]};}
Delete {Surface {out[3]};}
Delete {Surface {out[4]};}
Delete {Surface {out[5]};}

S4 = news; Plane Surface(S4) = {LL4};
out[] = Extrude {cs_length, 0., 0.} { Surface{S4}; Layers{{1}, {1}}; };
S2 = out[0];
Delete {Volume {out[1]};}
Delete {Surface {out[2]};}
Delete {Surface {out[3]};}
Delete {Surface {out[4]};}
Delete {Surface {out[5]};}

S5 = news; Plane Surface(S5) = {LL5};
out[] = Extrude {0., 0., cs_height} { Surface{S5}; Layers{{1}, {1}}; };
S6 = out[0];
Delete {Volume {out[1]};}
Delete {Surface {out[2]};}
Delete {Surface {out[3]};}
Delete {Surface {out[4]};}
Delete {Surface {out[5]};}

Geometry.AutoCoherence = 1;
Coherence;

Surface Loop(88) = {87, 19, 65, 64, 41, 42};
Volume(89) = {88};
