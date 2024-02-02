//+
Point(1) = {-2, 2, 0, 1.0};
//+
Point(2) = {2, 2, 0, 1.0};
//+
Point(3) = {2, -2, 0, 1.0};
//+
Point(4) = {-2, -2, 0, 1.0};
//+
Point(5) = {-300, 0, 0, 1.0};
//+
Point(6) = {300, 0, 0, 1.0};
//+
Point(7) = {0, 300, 0, 1.0};
//+
Point(8) = {0, -300, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Point(9) = {0, 0, 0, 1.0};
//+
Point(10) = {180, 240, 0, 1.0};
//+
Point(11) = {180, -240, 0, 1.0};
//+
Point(12) = {-180, -240, 0, 1.0};
//+
Point(13) = {-180, 240, 0, 1.0};
//+
Line(7) = {4, 12};
//+
Line(8) = {3, 11};
//+
Line(9) = {2, 10};
//+
Line(10) = {1, 13};
//+
Point(14) = {300, 0, 0, 1.0};
//+
Line(11) = {9, 14};
//+
Circle(12) = {12, 9, 13};
//+
Circle(13) = {13, 9, 10};
//+
Circle(14) = {10, 9, 11};
//+
Circle(15) = {11, 9, 12};
//+
Curve Loop(1) = {4, 10, -12, -7};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {-10, 1, 9,  -13};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {2, 8, -14, -9};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {3, 7, -15, -8};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {-4, -3, -2, -1};
//+
Plane Surface(5) = {5};
//+
Transfinite Curve {4, 12, 1, 13, 2, 14, 3, 15} = 11 Using Progression 1;
//+
Transfinite Curve {7, 10, 9, 8} = 10 Using Progression 1;
//+
Transfinite Surface {1} = {13, 1, 4, 12};
//+
Transfinite Surface {2} = {13, 10, 2, 1};
//+
Transfinite Surface {3} = {2, 10, 11, 3};
//+
Transfinite Surface {4} = {3, 11, 12, 4};
//+
Transfinite Surface {5} = {4, 1, 2, 3};
//+
Recombine Surface {1, 2, 3, 4, 5};
//+
Physical Curve("NRBC", 16) = {12, 13, 14, 15};
//+
Physical Curve("Outer", 17) = {4, 1, 2, 3};
//+
Physical Curve("Test", 18) = {11};
//+
Physical Surface("fluid", 19) = {2, 3, 4, 1, 5};

Mesh.Format = 1;
Mesh.MshFileVersion = 2.2 ;
