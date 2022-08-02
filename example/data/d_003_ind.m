% 3-bus test system with induction motor                       

Bus.con = [ ...
   1  100.00  1.00000  0.00000  1  1;
   2  100.00  1.00000 -0.09646  1  1;
   3  100.00  1.00000 -0.13938  1  1;
  ];

SW.con = [ ...
    2 100.0 100.00  1.06000  0.00000 99.00000 -99.00000 1.1 0.9  0 1 1 1;
  ];

PV.con = [ ...
    3 100.0 100.00  0.50000  1.03000  0.20000 -0.10000 1.1 0.9 1  1;
  ];

PQ.con = [ ...
   1 100.0 100.00  1.00000  0.20000 2.0 0.0 0  1;
   3 100.0 100.00  0.30000  0.10000 1.1 0.9 1  1;
  ];

Shunt.con = [ ...
  1 100.0 100.00 60  0.00000  0.19000  1;
  3 100.0 100.00 60  0.00000  0.04300  1;
  ];

Line.con = [ ...
   1    2  100.00 100.00 60 0   0.0000  0.01920  0.20000  0.05280  0.00000  0.00000 1.000    0.000    0.000  1;
   1    3  100.00 100.00 60 0   0.0000  0.04520  0.20000  0.04080  0.00000  0.00000 1.000    0.000    0.000  1;
   2    3  100.00 100.00 60 0   0.0000  0.05700  0.20000  0.03680  0.00000  0.00000 1.000    0.000    0.000  1;
  ];

Ltc.con=[ ...
];

Demand.con = [ ... 
 ];

Supply.con = [ ...  
];

Ind.con = [... 
  2  100  100  60  3  0  0.01  0.15  0.05  0.15  0.001  0.04  500  3  0.0  0.26  -0.13  0  0  1;
];



Bus.names = { ...
  'bus-1 '; 'bus-2 '; 'bus-3 '};

