Point(1)  = {-1.0    , 0.0   , 0.0, 1};
Point(2)  = { 4.0    , 0.0   , 0.0, 1};

Point(3)  = {-1.0    ,-2.0   , 0.0, 1};
Point(4)  = { 4.0    ,-2.0   , 0.0, 1};

Point(5)  = { 1.0    ,-2.0   , 0.0, 1};
Point(6)  = { 2.0    ,-2.0   , 0.0, 1};

Point(7)  = { 1.0    ,-1.9   , 0.0, 1};
Point(8)  = { 2.0    ,-1.9   , 0.0, 1};

Point(9)  = { 1.0    , 0.0   , 0.0, 1};
Point(10) = { 2.0    , 0.0   , 0.0, 1};

Point(11) = {-1.0    ,-1.9   , 0.0, 1};
Point(12) = { 4.0    ,-1.9   , 0.0, 1};

// Create lines for the points
Line(1)  = {3 ,5 };
Line(2)  = {5 ,7 };
Line(3)  = {7 ,8 };
Line(4)  = {8 ,6 };
Line(5)  = {6 ,4 };
Line(6)  = {4 ,12};
Line(7)  = {12,2 };
Line(8)  = {2 ,10};
Line(9)  = {10,9 };
Line(10) = {9 ,1 };
Line(11) = {1 ,11};
Line(12) = {11,3 };

Line(13) = {7 ,9 };
Line(14) = {8 ,10};

Line(15) = {7 ,11};
Line(16) = {12,8 };

Transfinite Line { 1, 10, 15}    =30 Using Progression 1;
Transfinite Line { 5, 16,  8}    =30 Using Progression 1;
Transfinite Line { 3, 9     }    =18 Using Progression 1;
Transfinite Line {7, 13, 14, 11} =45 Using Progression 1;
Transfinite Line {2, 4, 12, 6  } = 4 Using Progression 1;


Line Loop(17) = {5, 6, 16, 4};
Plane Surface(18) = {17};
Transfinite Surface {18};

Line Loop(19) = {-16, -14, 8, 7};
Plane Surface(20) = {19};
Transfinite Surface {20};

Line Loop(21) = {3, 14, 9, -13};
Plane Surface(22) = {21};
Transfinite Surface {22};

Line Loop(23) = {-15, 11, 10, 13};
Plane Surface(24) = {23};
Transfinite Surface {24};

Line Loop(25) = {1, 2, 15, 12};
Plane Surface(26) = {25};
Transfinite Surface {26};

Physical Surface("fluid") = {18, 20, 22, 24, 26};
Physical Point("point") = {1};
Physical Line("inflow") = {11, 12};
Physical Line("outflow") = {7, 6};
