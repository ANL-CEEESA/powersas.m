% 3-bus test system with synchronous generators (with controllers), const-impedance load and induction motors

Bus.con = [ ...
   1  100.00  1.00000  0.00000  1  1;
   2  100.00  1.00000 -0.09646  1  1;
   3  100.00  1.00000 -0.13938  1  1;
  ];

SW.con = [ ...
    2 100.0 100.00  1.00000  0.00000 99.00000 -99.00000 1.1 0.9  0 1 1 1;
  ];

PV.con = [ ...
    3 100.0 100.00  0.50000  1.00000  0.20000 -0.10000 1.1 0.9 1  1;
  ];

PQ.con = [ ...
   1 100.0 100.00  0.00000  0.00000 2.0 0.0 0  1;
  ];
  
Shunt.con = [ ...
  1 100.0 100.00 60  1.00000  0.20000  1;
  3 100.0 100.00 60  0.40000  0.10000  1;
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

Syn.con = [ ...      xl    ra    xd    xd'   xd"   Td'  Td"   xq     xq'    xq"   Tq'  Tq"   M       D  |-> optional
  2  100 100  60  3  0.05  0.01  0.50  0.30  0.23  7.4  0.03  0.69  0.446  0.4   0.3 0.033  2.296  2  0  0  1  1  0  0  0;
  3  100 100  60  3  0.00  0.01  0.35  0.20  0.13  6.1  0.04  0.50  0.360  0.13  0.3 0.099  3.080  2  0  0  1  1  0  0  0;
 ];


Ind.con = [...
  1  100  100  60  3  0  0.01  0.15  0.05  0.15  0.001  0.04  500  0.08  0.35  0  0  0  0  1; 
  2  100  100  60  3  0  0.01  0.15  0.05  0.15  0.001  0.04  500  2.25  0.33  0  0  0  0  1; 
];

Tg.con = [ ...
    1       2       0      0.05    10    -10   0.1      0.3;
    2       2       0      0.05    10    -10   0.1      0.3
    ];

Exc.con = [ ...
    2       3       5      -5      20     0.02   0.02    1.2    1.0   0.004   0.15  0   0   1;
    1       3       5      -5      20     0.02   0.02    1.2    1.0   0.004   0.15  0   0   1;
    ];



Bus.names = { ...
  'bus-1 '; 'bus-2 '; 'bus-3 '};

