//+
/*
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Transfinite Curve {4, 1, 2, 3} = 101 Using Progression 1;
//+
Transfinite Surface {1};
Recombine Surface {1};
//+
Physical Curve("left") = {4};
//+
Physical Curve("right") = {2};
//+
Physical Curve("top") = {3};
//+
Physical Curve("bottom") = {1};
//+
Physical Surface("fluid") = {1};
*/

Mesh.Format = 1; 
Mesh.MshFileVersion = 2.2;
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.5, 0, 0, 1.0};
//+
Point(3) = {0, 0.5, 0, 1.0};
//+
Point(4) = {0, -0.5, 0, 1.0};
//+
Point(5) = {-0.5, 0, 0, 1.0};
//+
Circle(1) = {2, 1, 3};
//+
Circle(2) = {3, 1, 5};
//+
Circle(3) = {5, 1, 4};
//+
Circle(4) = {4, 1, 2};
//+
Point(6) = {13, 0, 0, 1.0};
//+
Point(7) = {-13, 0, 0, 1.0};
//+
Point(8) = {0, 13, 0, 1.0};
//+
Point(9) = {0, -13, 0, 1.0};
//+
Circle(5) = {6, 1, 8};
//+
Circle(6) = {8, 1, 7};
//+
Circle(7) = {7, 1, 9};
//+
Circle(8) = {9, 1, 6};
//+
Line(9) = {3, 8};
//+
Line(10) = {2, 6};
//+
Line(11) = {4, 9};
//+
Line(12) = {5, 7};
//+
Curve Loop(1) = {10, 5, -9, -1};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {-10, -4, 11, 8};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {-11, -3, 12, 7};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {9, 6, -12, -2};
//+
Plane Surface(4) = {4};
//+
Transfinite Curve {10, 9, 12, 11} = 250 Using Progression 1;
//+
Transfinite Curve {4, 8, 1, 5, 2, 6, 3, 7} = 10 Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};

Recombine Surface {1,2,3,4};
//+
Physical Curve("cyl") = {1,2,3,4};
//+
Physical Curve("outercyl") = {6, 5, 8, 7};
//+
Physical Curve("radial") = {10};
//+
Physical Surface("fluid") = {4, 1, 2, 3};
