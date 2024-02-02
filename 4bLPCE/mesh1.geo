//+
Point(1) = {100,100, 0, 1.0};
//+
Point(2) = {-100, -100, 0, 1.0};
//+
Point(3) = {-100, 100, 0, 1.0};
//+
Point(4) = {100, -100, 0, 1.0};
//+
Point(5) = {100, 0, 0, 1.0};
//+
Point(6) = {-100, 0, 0, 1.0};
//+
Point(7) = {0, 100, 0, 1.0};
//+
Point(8) = {0, -100, 0, 1.0};
//+
Point(9) = {0, 0, 0, 1.0};
//+
Point(10) = {150, 0, 0, 1.0};
//+
Point(11) = {-150, 0, 0, 1.0};
//+
Point(12) = {0, 150, 0, 1.0};
//+
Point(13) = {0, -150, 0, 1.0};
//+
Line(1) = {3, 6};
//+
Line(2) = {6, 9};
//+
Line(3) = {9, 7};
//+
Line(4) = {7, 3};
//+
Line(5) = {6, 2};
//+
Line(6) = {2, 8};
//+
Line(7) = {8, 9};
//+
Line(8) = {8, 4};
//+
Line(9) = {4, 5};
//+
Line(10) = {5, 9};
//+
Line(11) = {5, 1};
//+
Line(12) = {1, 7};
//+
Line(17) = {5, 10};
//+
Line(18) = {8, 13};
//+
Line(19) = {6, 11};
//+
x = 150/(2^0.5);
//+
Point(14) = {x, x, 0, 1.0};
//+
Point(15) = {x, -x, 0, 1.0};
//+
Point(16) = {-x, -x, 0, 1.0};
//+
Point(17) = {-x, x, 0, 1.0};
//+
Line(20) = {1, 14};
//+
Line(21) = {7, 12};
//+
Line(22) = {3, 17};
//+
Line(23) = {2, 16};
//+
Line(24) = {4, 15};
//+
Circle(25) = {13, 9, 16};
//+
Circle(26) = {16, 9, 11};
//+
Circle(27) = {11, 9, 17};
//+
Circle(28) = {17, 9, 12};
//+
Circle(29) = {12, 9, 14};
//+
Circle(30) = {14, 9, 10};
//+
Circle(31) = {10, 9, 15};
//+
Circle(32) = {15, 9, 13};
//+
Curve Loop(1) = {2, 3, 4, 1};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {-10, 11, 12, -3};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {8, 9, 10, -7};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {6, 7, -2, 5};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {-32, -24, -8, 18};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {-31, -17, -9, 24};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {17, -30, -20, -11};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {-12, 20, -29, -21};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {-4, 21, -28, -22};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {-19, -1, 22, -27};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {-23, -5, 19, -26};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {-25, -18, -6, 23};
//+
Plane Surface(12) = {12};
//+
Transfinite Curve {25, 32, 31, 30, 29, 28, 27, 26, 1, 5, 6, 8, 9, 11, 12, 4, 3, 7, 10, 2} = 41 Using Progression 1;
//+
Transfinite Curve {23, 24, 18, 17, 20, 21, 22, 19} = 21 Using Progression 1; //outer nrbc condition 
//+
Transfinite Surface {1} = {6, 9, 7, 3};
//+
Transfinite Surface {9} = {3, 7, 12, 17};
//+
Transfinite Surface {10} = {11, 6, 3, 17};
//+
Transfinite Surface {11} = {16, 2, 6, 11};
//+
Transfinite Surface {4} = {2, 8, 9, 6};
//+
Transfinite Surface {12} = {16, 13, 8, 2};
//+
Transfinite Surface {5} = {13, 15, 4, 8};
//+
Transfinite Surface {3} = {8, 4, 5, 9};
//+
Transfinite Surface {6} = {4, 15, 10, 5};
//+
Transfinite Surface {2} = {9, 5, 1, 7};
//+
Transfinite Surface {7} = {5, 10, 14, 1};
//+
Transfinite Surface {8} = {7, 1, 14, 12};
//+
Recombine Surface {1, 2, 3, 4, 5,6,7,8,9,10,11,12};
//+
Physical Curve("outer", 33) = {6, 8, 9, 11, 12, 4, 1, 5};
//+
Physical Curve("nrbc", 34) = {32, 31, 30, 29, 28, 27, 26, 25};
//+
Physical Curve("test", 35) = {10, 17};
//+
Physical Surface("fluid", 36) = {9, 8, 7, 6, 5, 12, 11, 10, 1, 2, 3, 4};

Mesh.Format = 1;
Mesh.MshFileVersion = 2.2 ;

