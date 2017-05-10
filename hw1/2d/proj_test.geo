Point(1) = {0.0    , 0.0   , 0.0, 1};
Point(2) = {2.0    , 0.0   , 0.0, 1};
Point(3) = {3.4142 ,-1.4142, 0.0, 1};
Point(4) = {1.4142 ,-1.4142, 0.0, 1};

Point(5) = {0.6364 ,-0.6364, 0.0, 1};
Point(6) = {1.1364 ,-0.6364, 0.0, 1};
Point(7) = {1.4192 ,-0.9192, 0.0, 1};
Point(8) = {0.9192 ,-0.9192, 0.0, 1};

Point(9)  = {0.5    , 0.0   , 0.0, 1};
Point(10) = {1.9192 ,-1.4192, 0.0, 1};

Point(11) = {2.6364 ,-0.6364, 0.0, 1};
Point(12) = {2.9192 ,-0.9192, 0.0, 1};

// Create lines for the points
Line(4) = {4 ,8 };
Line(5) = {8 ,7 };
Line(6) = {7 ,6 };
Line(7) = {6 ,5 };

Line(8) = {5 ,1 };

Line(9)  = {6 ,9 };
Line(10) = {10,7  };

Line(11) = {1, 9  };
Line(12) = {9, 2  };
Line(13) = {3, 10 };
Line(14) = {10, 4 };

Line(15) = {12, 7 };
Line(16) = { 6, 11};

Line(17) = { 3, 12};
Line(18) = {12, 11};
Line(19) = {11, 2 };

Transfinite Line { 5, 7, 11, 14} = 3 Using Progression 1;
Transfinite Line { 4, 10, 17   } = 3 Using Progression 1;

Transfinite Line {  8, 9, 19   } = 4 Using Progression 1;

Transfinite Line { 6, 18       } = 2 Using Progression 1;
Transfinite Line {12, 13, 15,16} =10 Using Progression 1;

Line Loop(5)  = {-14, 10, -5, -4};
Plane Surface(15) = {5};
Transfinite Surface {15};

Line Loop(6)  = {-7 , 9, -11, -8};
Plane Surface(16) = {6};
Transfinite Surface {16};

Line Loop(7)  = {-13, 17,  15, -10};
Plane Surface(17) = {7};
Transfinite Surface {17};

Line Loop(8)  = { 18,-16,  -6, -15};
Plane Surface(18) = {8};
Transfinite Surface {18};

Line Loop(9)  = { 19,-12,  -9,  16};
Plane Surface(19) = {9};
Transfinite Surface {19};

Physical Surface(20) = {16, 19, 18, 17, 15};
Physical Line(21) = {4, 5, 6, 7, 8};
