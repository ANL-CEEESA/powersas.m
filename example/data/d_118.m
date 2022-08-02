% IEEE 118-bus test system
% From Power Systems Test Case Archive, University of Washington
% https://labs.ece.uw.edu/pstca/pf118/pg_tca118bus.htm
%
% Converted from IEEE CDF format to PSAT (MATLAB) format using PSAT (MATLAB)

Bus.con = [ ...
   1  100.00  0.95500  0.19150  1  1;
   2  100.00  0.97100  0.20092  1  1;
   3  100.00  0.96700  0.20693  1  1;
   4  100.00  0.99800  0.27182  1  1;
   5  100.00  1.00200  0.27958  1  1;
   6  100.00  0.99000  0.23197  1  1;
   7  100.00  0.98900  0.22422  1  1;
   8  100.00  1.01500  0.36722  1  1;
   9  100.00  1.04200  0.49382  1  1;
  10  100.00  1.05000  0.62614  1  1;
  11  100.00  0.98500  0.22698  1  1;
  12  100.00  0.99000  0.21796  1  1;
  13  100.00  0.96800  0.20296  1  1;
  14  100.00  0.98300  0.20544  1  1;
  15  100.00  0.97000  0.20026  1  1;
  16  100.00  0.98300  0.21270  1  1;
  17  100.00  0.99500  0.24426  1  1;
  18  100.00  0.97300  0.20560  1  1;
  19  100.00  0.96200  0.19747  1  1;
  20  100.00  0.95600  0.21277  1  1;
  21  100.00  0.95700  0.24047  1  1;
  22  100.00  0.96900  0.28503  1  1;
  23  100.00  0.99900  0.37085  1  1;
  24  100.00  0.99200  0.36849  1  1;
  25  100.00  1.05000  0.49182  1  1;
  26  100.00  1.01500  0.52290  1  1;
  27  100.00  0.96800  0.27234  1  1;
  28  100.00  0.96100  0.24222  1  1;
  29  100.00  0.96300  0.22489  1  1;
  30  100.00  0.98500  0.33219  1  1;
  31  100.00  0.96700  0.22691  1  1;
  32  100.00  0.96300  0.26285  1  1;
  33  100.00  0.97000  0.18942  1  1;
  34  100.00  0.98400  0.20090  1  1;
  35  100.00  0.98000  0.19295  1  1;
  36  100.00  0.98000  0.19295  1  1;
  37  100.00  0.99000  0.20885  1  1;
  38  100.00  0.96100  0.29857  1  1;
  39  100.00  0.97000  0.14969  1  1;
  40  100.00  0.97000  0.13082  1  1;
  41  100.00  0.96600  0.12307  1  1;
  42  100.00  0.98500  0.15102  1  1;
  43  100.00  0.97700  0.20001  1  1;
  44  100.00  0.98400  0.24335  1  1;
  45  100.00  0.98600  0.27527  1  1;
  46  100.00  1.00500  0.32419  1  1;
  47  100.00  1.01700  0.36301  1  1;
  48  100.00  1.02000  0.34938  1  1;
  49  100.00  1.02500  0.36689  1  1;
  50  100.00  1.00100  0.33130  1  1;
  51  100.00  0.96600  0.28561  1  1;
  52  100.00  0.95600  0.26896  1  1;
  53  100.00  0.94600  0.25196  1  1;
  54  100.00  0.95500  0.26787  1  1;
  55  100.00  0.95200  0.26281  1  1;
  56  100.00  0.95400  0.26606  1  1;
  57  100.00  0.97000  0.28709  1  1;
  58  100.00  0.95900  0.27213  1  1;
  59  100.00  0.98500  0.33943  1  1;
  60  100.00  0.99300  0.40544  1  1;
  61  100.00  0.99500  0.42099  1  1;
  62  100.00  0.99800  0.41022  1  1;
  63  100.00  0.96800  0.39841  1  1;
  64  100.00  0.98300  0.42923  1  1;
  65  100.00  1.00500  0.48379  1  1;
  66  100.00  1.05000  0.48098  1  1;
  67  100.00  1.01900  0.43492  1  1;
  68  100.00  1.00300  0.48166  1  1;
  69  100.00  1.03500  0.52360  1  1;
  70  100.00  0.98400  0.39474  1  1;
  71  100.00  0.98600  0.38757  1  1;
  72  100.00  0.98000  0.36840  1  1;
  73  100.00  0.99100  0.38389  1  1;
  74  100.00  0.95800  0.37818  1  1;
  75  100.00  0.96700  0.40020  1  1;
  76  100.00  0.94300  0.38045  1  1;
  77  100.00  1.00600  0.46688  1  1;
  78  100.00  1.00300  0.46157  1  1;
  79  100.00  1.00900  0.46679  1  1;
  80  100.00  1.04000  0.50597  1  1;
  81  100.00  0.99600  0.49121  1  1;
  82  100.00  0.98800  0.47597  1  1;
  83  100.00  0.98400  0.49679  1  1;
  84  100.00  0.97900  0.54105  1  1;
  85  100.00  0.98500  0.56819  1  1;
  86  100.00  0.98600  0.54430  1  1;
  87  100.00  1.01500  0.54882  1  1;
  88  100.00  0.98700  0.62291  1  1;
  89  100.00  1.00500  0.69373  1  1;
  90  100.00  0.98500  0.58186  1  1;
  91  100.00  0.98000  0.58207  1  1;
  92  100.00  0.99000  0.59132  1  1;
  93  100.00  0.98500  0.53842  1  1;
  94  100.00  0.98900  0.50060  1  1;
  95  100.00  0.98000  0.48361  1  1;
  96  100.00  0.99200  0.48070  1  1;
  97  100.00  1.01100  0.48721  1  1;
  98  100.00  1.02300  0.47880  1  1;
  99  100.00  1.01000  0.47239  1  1;
 100  100.00  1.01700  0.48970  1  1;
 101  100.00  0.99100  0.51742  1  1;
 102  100.00  0.98900  0.56488  1  1;
 103  100.00  1.01000  0.42441  1  1;
 104  100.00  0.97100  0.37956  1  1;
 105  100.00  0.96500  0.36029  1  1;
 106  100.00  0.96100  0.35575  1  1;
 107  100.00  0.95200  0.30686  1  1;
 108  100.00  0.96600  0.33934  1  1;
 109  100.00  0.96700  0.33144  1  1;
 110  100.00  0.97300  0.31667  1  1;
 111  100.00  0.98000  0.34538  1  1;
 112  100.00  0.97500  0.26257  1  1;
 113  100.00  0.99300  0.24421  1  1;
 114  100.00  0.96000  0.25702  1  1;
 115  100.00  0.96000  0.25688  1  1;
 116  100.00  1.00500  0.47407  1  1;
 117  100.00  0.97300  0.19106  1  1;
 118  100.00  0.94900  0.38294  1  1;
  ];

SW.con = [ ...
  69 100.0 100.00  1.13500  0.00000  3.00000 -3.00000 1.1 0.9  5.13863 1 1 1;
  ];

PV.con = [ ...
   1 100.0 100.00  0.00000  0.95500  0.15000 -0.05000 1.1 0.9 1  1;
   4 100.0 100.00 -0.09000  0.99800  3.00000 -3.00000 1.1 0.9 1  1;
   6 100.0 100.00  0.00000  0.99000  0.50000 -0.13000 1.1 0.9 1  1;
   8 100.0 100.00 -0.28000  1.01500  3.00000 -3.00000 1.1 0.9 1  1;
  10 100.0 100.00  4.50000  1.05000  2.00000 -1.47000 1.1 0.9 1  1;
  12 100.0 100.00  0.85000  0.99000  1.20000 -0.35000 1.1 0.9 1  1;
  15 100.0 100.00  0.00000  0.97000  0.30000 -0.10000 1.1 0.9 1  1;
  18 100.0 100.00  0.00000  0.97300  0.50000 -0.16000 1.1 0.9 1  1;
  19 100.0 100.00  0.00000  0.96200  0.24000 -0.08000 1.1 0.9 1  1;
  24 100.0 100.00 -0.13000  0.99200  3.00000 -3.00000 1.1 0.9 1  1;
  25 100.0 100.00  2.20000  1.05000  1.40000 -0.47000 1.1 0.9 1  1;
  26 100.0 100.00  3.14000  1.01500 10.00000 -10.00000 1.1 0.9 1  1;
  27 100.0 100.00 -0.09000  0.96800  3.00000 -3.00000 1.1 0.9 1  1;
  31 100.0 100.00  0.07000  0.96700  3.00000 -3.00000 1.1 0.9 1  1;
  32 100.0 100.00  0.00000  0.96300  0.42000 -0.14000 1.1 0.9 1  1;
  34 100.0 100.00  0.00000  0.98400  0.24000 -0.08000 1.1 0.9 1  1;
  36 100.0 100.00  0.00000  0.98000  0.24000 -0.08000 1.1 0.9 1  1;
  40 100.0 100.00 -0.46000  0.97000  3.00000 -3.00000 1.1 0.9 1  1;
  42 100.0 100.00 -0.59000  0.98500  3.00000 -3.00000 1.1 0.9 1  1;
  46 100.0 100.00  0.19000  1.00500  1.00000 -1.00000 1.1 0.9 1  1;
  49 100.0 100.00  2.04000  1.02500  2.10000 -0.85000 1.1 0.9 1  1;
  54 100.0 100.00  0.48000  0.95500  3.00000 -3.00000 1.1 0.9 1  1;
  55 100.0 100.00  0.00000  0.95200  0.23000 -0.08000 1.1 0.9 1  1;
  56 100.0 100.00  0.00000  0.95400  0.15000 -0.08000 1.1 0.9 1  1;
  59 100.0 100.00  1.55000  0.98500  1.80000 -0.60000 1.1 0.9 1  1;
  61 100.0 100.00  1.60000  0.99500  3.00000 -1.00000 1.1 0.9 1  1;
  62 100.0 100.00  0.00000  0.99800  0.20000 -0.20000 1.1 0.9 1  1;
  65 100.0 100.00  3.91000  1.00500  2.00000 -0.67000 1.1 0.9 1  1;
  66 100.0 100.00  3.92000  1.05000  2.00000 -0.67000 1.1 0.9 1  1;
  70 100.0 100.00  0.00000  0.98400  0.32000 -0.10000 1.1 0.9 1  1;
  72 100.0 100.00 -0.12000  0.98000  1.00000 -1.00000 1.1 0.9 1  1;
  73 100.0 100.00 -0.06000  0.99100  1.00000 -1.00000 1.1 0.9 1  1;
  74 100.0 100.00  0.00000  0.95800  0.09000 -0.06000 1.1 0.9 1  1;
  76 100.0 100.00  0.00000  0.94300  0.23000 -0.08000 1.1 0.9 1  1;
  77 100.0 100.00  0.00000  1.00600  0.70000 -0.20000 1.1 0.9 1  1;
  80 100.0 100.00  4.77000  1.04000  2.80000 -1.65000 1.1 0.9 1  1;
  85 100.0 100.00  0.00000  0.98500  0.23000 -0.08000 1.1 0.9 1  1;
  87 100.0 100.00  0.04000  1.01500 10.00000 -1.00000 1.1 0.9 1  1;
  89 100.0 100.00  6.07000  1.00500  3.00000 -2.10000 1.1 0.9 1  1;
  90 100.0 100.00 -0.85000  0.98500  3.00000 -3.00000 1.1 0.9 1  1;
  91 100.0 100.00 -0.10000  0.98000  1.00000 -1.00000 1.1 0.9 1  1;
  92 100.0 100.00  0.00000  0.99000  0.09000 -0.03000 1.1 0.9 1  1;
  99 100.0 100.00 -0.42000  1.01000  1.00000 -1.00000 1.1 0.9 1  1;
 100 100.0 100.00  2.52000  1.01700  1.55000 -0.50000 1.1 0.9 1  1;
 103 100.0 100.00  0.40000  1.01000  0.40000 -0.15000 1.1 0.9 1  1;
 104 100.0 100.00  0.00000  0.97100  0.23000 -0.08000 1.1 0.9 1  1;
 105 100.0 100.00  0.00000  0.96500  0.23000 -0.08000 1.1 0.9 1  1;
 107 100.0 100.00 -0.22000  0.95200  2.00000 -2.00000 1.1 0.9 1  1;
 110 100.0 100.00  0.00000  0.97300  0.23000 -0.08000 1.1 0.9 1  1;
 111 100.0 100.00  0.36000  0.98000 10.00000 -1.00000 1.1 0.9 1  1;
 112 100.0 100.00 -0.43000  0.97500 10.00000 -1.00000 1.1 0.9 1  1;
 113 100.0 100.00 -0.06000  0.99300  2.00000 -1.00000 1.1 0.9 1  1;
 116 100.0 100.00 -1.84000  1.00500 10.00000 -10.00000 1.1 0.9 1  1;
  ];

PQ.con = [ ...
   1 100.0 100.00  0.51000  0.27000 1.1 0.9 1  1;
   2 100.0 100.00  0.20000  0.09000 1.1 0.9 1  1;
   3 100.0 100.00  0.39000  0.10000 1.1 0.9 1  1;
   4 100.0 100.00  0.30000  0.12000 1.1 0.9 1  1;
   6 100.0 100.00  0.52000  0.22000 1.1 0.9 1  1;
   7 100.0 100.00  0.19000  0.02000 1.1 0.9 1  1;
  11 100.0 100.00  0.70000  0.23000 1.1 0.9 1  1;
  12 100.0 100.00  0.47000  0.10000 1.1 0.9 1  1;
  13 100.0 100.00  0.34000  0.16000 1.1 0.9 1  1;
  14 100.0 100.00  0.14000  0.01000 1.1 0.9 1  1;
  15 100.0 100.00  0.90000  0.30000 1.1 0.9 1  1;
  16 100.0 100.00  0.25000  0.10000 1.1 0.9 1  1;
  17 100.0 100.00  0.11000  0.03000 1.1 0.9 1  1;
  18 100.0 100.00  0.60000  0.34000 1.1 0.9 1  1;
  19 100.0 100.00  0.45000  0.25000 1.1 0.9 1  1;
  20 100.0 100.00  0.18000  0.03000 1.1 0.9 1  1;
  21 100.0 100.00  0.14000  0.08000 1.1 0.9 1  1;
  22 100.0 100.00  0.10000  0.05000 1.1 0.9 1  1;
  23 100.0 100.00  0.07000  0.03000 1.1 0.9 1  1;
  27 100.0 100.00  0.62000  0.13000 1.1 0.9 1  1;
  28 100.0 100.00  0.17000  0.07000 1.1 0.9 1  1;
  29 100.0 100.00  0.24000  0.04000 1.1 0.9 1  1;
  31 100.0 100.00  0.43000  0.27000 1.1 0.9 1  1;
  32 100.0 100.00  0.59000  0.23000 1.1 0.9 1  1;
  33 100.0 100.00  0.23000  0.09000 1.1 0.9 1  1;
  34 100.0 100.00  0.59000  0.26000 1.1 0.9 1  1;
  35 100.0 100.00  0.33000  0.09000 1.1 0.9 1  1;
  36 100.0 100.00  0.31000  0.17000 1.1 0.9 1  1;
  39 100.0 100.00  0.27000  0.11000 1.1 0.9 1  1;
  40 100.0 100.00  0.20000  0.23000 1.1 0.9 1  1;
  41 100.0 100.00  0.37000  0.10000 1.1 0.9 1  1;
  42 100.0 100.00  0.37000  0.23000 1.1 0.9 1  1;
  43 100.0 100.00  0.18000  0.07000 1.1 0.9 1  1;
  44 100.0 100.00  0.16000  0.08000 1.1 0.9 1  1;
  45 100.0 100.00  0.53000  0.22000 1.1 0.9 1  1;
  46 100.0 100.00  0.28000  0.10000 1.1 0.9 1  1;
  47 100.0 100.00  0.34000  0.00000 1.1 0.9 1  1;
  48 100.0 100.00  0.20000  0.11000 1.1 0.9 1  1;
  49 100.0 100.00  0.87000  0.30000 1.1 0.9 1  1;
  50 100.0 100.00  0.17000  0.04000 1.1 0.9 1  1;
  51 100.0 100.00  0.17000  0.08000 1.1 0.9 1  1;
  52 100.0 100.00  0.18000  0.05000 1.1 0.9 1  1;
  53 100.0 100.00  0.23000  0.11000 1.1 0.9 1  1;
  54 100.0 100.00  1.13000  0.32000 1.1 0.9 1  1;
  55 100.0 100.00  0.63000  0.22000 1.1 0.9 1  1;
  56 100.0 100.00  0.84000  0.18000 1.1 0.9 1  1;
  57 100.0 100.00  0.12000  0.03000 1.1 0.9 1  1;
  58 100.0 100.00  0.12000  0.03000 1.1 0.9 1  1;
  59 100.0 100.00  2.77000  1.13000 1.1 0.9 1  1;
  60 100.0 100.00  0.78000  0.03000 1.1 0.9 1  1;
  62 100.0 100.00  0.77000  0.14000 1.1 0.9 1  1;
  66 100.0 100.00  0.39000  0.18000 1.1 0.9 1  1;
  67 100.0 100.00  0.28000  0.07000 1.1 0.9 1  1;
  70 100.0 100.00  0.66000  0.20000 1.1 0.9 1  1;
  74 100.0 100.00  0.68000  0.27000 1.1 0.9 1  1;
  75 100.0 100.00  0.47000  0.11000 1.1 0.9 1  1;
  76 100.0 100.00  0.68000  0.36000 1.1 0.9 1  1;
  77 100.0 100.00  0.61000  0.28000 1.1 0.9 1  1;
  78 100.0 100.00  0.71000  0.26000 1.1 0.9 1  1;
  79 100.0 100.00  0.39000  0.32000 1.1 0.9 1  1;
  80 100.0 100.00  1.30000  0.26000 1.1 0.9 1  1;
  82 100.0 100.00  0.54000  0.27000 1.1 0.9 1  1;
  83 100.0 100.00  0.20000  0.10000 1.1 0.9 1  1;
  84 100.0 100.00  0.11000  0.07000 1.1 0.9 1  1;
  85 100.0 100.00  0.24000  0.15000 1.1 0.9 1  1;
  86 100.0 100.00  0.21000  0.10000 1.1 0.9 1  1;
  88 100.0 100.00  0.48000  0.10000 1.1 0.9 1  1;
  90 100.0 100.00  0.78000  0.42000 1.1 0.9 1  1;
  92 100.0 100.00  0.65000  0.10000 1.1 0.9 1  1;
  93 100.0 100.00  0.12000  0.07000 1.1 0.9 1  1;
  94 100.0 100.00  0.30000  0.16000 1.1 0.9 1  1;
  95 100.0 100.00  0.42000  0.31000 1.1 0.9 1  1;
  96 100.0 100.00  0.38000  0.15000 1.1 0.9 1  1;
  97 100.0 100.00  0.15000  0.09000 1.1 0.9 1  1;
  98 100.0 100.00  0.34000  0.08000 1.1 0.9 1  1;
 100 100.0 100.00  0.37000  0.18000 1.1 0.9 1  1;
 101 100.0 100.00  0.22000  0.15000 1.1 0.9 1  1;
 102 100.0 100.00  0.05000  0.03000 1.1 0.9 1  1;
 103 100.0 100.00  0.23000  0.16000 1.1 0.9 1  1;
 104 100.0 100.00  0.38000  0.25000 1.1 0.9 1  1;
 105 100.0 100.00  0.31000  0.26000 1.1 0.9 1  1;
 106 100.0 100.00  0.43000  0.16000 1.1 0.9 1  1;
 107 100.0 100.00  0.28000  0.12000 1.1 0.9 1  1;
 108 100.0 100.00  0.02000  0.01000 1.1 0.9 1  1;
 109 100.0 100.00  0.08000  0.03000 1.1 0.9 1  1;
 110 100.0 100.00  0.39000  0.30000 1.1 0.9 1  1;
 112 100.0 100.00  0.25000  0.13000 1.1 0.9 1  1;
 114 100.0 100.00  0.08000  0.03000 1.1 0.9 1  1;
 115 100.0 100.00  0.22000  0.07000 1.1 0.9 1  1;
 117 100.0 100.00  0.20000  0.08000 1.1 0.9 1  1;
 118 100.0 100.00  0.33000  0.15000 1.1 0.9 1  1;
  ];

Shunt.con = [ ...
   5 100.0 100.00 60  0.00000 -0.40000  1;
  34 100.0 100.00 60  0.00000  0.14000  1;
  37 100.0 100.00 60  0.00000 -0.25000  1;
  44 100.0 100.00 60  0.00000  0.10000  1;
  45 100.0 100.00 60  0.00000  0.10000  1;
  46 100.0 100.00 60  0.00000  0.10000  1;
  48 100.0 100.00 60  0.00000  0.15000  1;
  74 100.0 100.00 60  0.00000  0.12000  1;
  79 100.0 100.00 60  0.00000  0.20000  1;
  82 100.0 100.00 60  0.00000  0.20000  1;
  83 100.0 100.00 60  0.00000  0.10000  1;
 105 100.0 100.00 60  0.00000  0.20000  1;
 107 100.0 100.00 60  0.00000  0.06000  1;
 110 100.0 100.00 60  0.00000  0.06000  1;
  ];

Line.con = [ ...
   2    1  100.00 100.00 60 0   0.0000  0.03030  0.09990  0.02540  0.00000  0.00000 0    0.000    0.000  1;% 1
   3    1  100.00 100.00 60 0   0.0000  0.01290  0.04240  0.01082  0.00000  0.00000 0    0.000    0.000  1;% 2
   5    4  100.00 100.00 60 0   0.0000  0.00176  0.00798  0.00210  0.00000  0.00000 0    0.000    0.000  1;% 3
   5    3  100.00 100.00 60 0   0.0000  0.02410  0.10800  0.02840  0.00000  0.00000 0    0.000    0.000  1;% 4
   6    5  100.00 100.00 60 0   0.0000  0.01190  0.05400  0.01426  0.00000  0.00000 0    0.000    0.000  1;% 5
   7    6  100.00 100.00 60 0   0.0000  0.00459  0.02080  0.00550  0.00000  0.00000 0    0.000    0.000  1;% 6
   9    8  100.00 100.00 60 0   0.0000  0.00244  0.03050  1.16200  0.00000  0.00000 0    0.000    0.000  1;% 7
   8    5  100.00 100.00 60 0   0.0000  0.00000  0.02670  0.00000  0.98500  0.00000 0    0.000    0.000  1;% 8
  10    9  100.00 100.00 60 0   0.0000  0.00258  0.03220  1.23000  0.00000  0.00000 0    0.000    0.000  1;% 9
  11    4  100.00 100.00 60 0   0.0000  0.02090  0.06880  0.01748  0.00000  0.00000 0    0.000    0.000  1;% 10
  11    5  100.00 100.00 60 0   0.0000  0.02030  0.06820  0.01738  0.00000  0.00000 0    0.000    0.000  1;% 11
  12   11  100.00 100.00 60 0   0.0000  0.00595  0.01960  0.00502  0.00000  0.00000 0    0.000    0.000  1;% 12
  12    2  100.00 100.00 60 0   0.0000  0.01870  0.06160  0.01572  0.00000  0.00000 0    0.000    0.000  1;% 13
  12    3  100.00 100.00 60 0   0.0000  0.04840  0.16000  0.04060  0.00000  0.00000 0    0.000    0.000  1;% 14
  12    7  100.00 100.00 60 0   0.0000  0.00862  0.03400  0.00874  0.00000  0.00000 0    0.000    0.000  1;% 15
  13   11  100.00 100.00 60 0   0.0000  0.02225  0.07310  0.01876  0.00000  0.00000 0    0.000    0.000  1;% 16
  14   12  100.00 100.00 60 0   0.0000  0.02150  0.07070  0.01816  0.00000  0.00000 0    0.000    0.000  1;% 17
  15   13  100.00 100.00 60 0   0.0000  0.07440  0.24440  0.06268  0.00000  0.00000 0    0.000    0.000  1;% 18
  15   14  100.00 100.00 60 0   0.0000  0.05950  0.19500  0.05020  0.00000  0.00000 0    0.000    0.000  1;% 19
  16   12  100.00 100.00 60 0   0.0000  0.02120  0.08340  0.02140  0.00000  0.00000 0    0.000    0.000  1;% 20
  17   15  100.00 100.00 60 0   0.0000  0.01320  0.04370  0.04440  0.00000  0.00000 0    0.000    0.000  1;% 21
  17   16  100.00 100.00 60 0   0.0000  0.04540  0.18010  0.04660  0.00000  0.00000 0    0.000    0.000  1;% 22
  18   17  100.00 100.00 60 0   0.0000  0.01230  0.05050  0.01298  0.00000  0.00000 0    0.000    0.000  1;% 23
  19   18  100.00 100.00 60 0   0.0000  0.01119  0.04930  0.01142  0.00000  0.00000 0    0.000    0.000  1;% 24
  20   19  100.00 100.00 60 0   0.0000  0.02520  0.11700  0.02980  0.00000  0.00000 0    0.000    0.000  1;% 25
  19   15  100.00 100.00 60 0   0.0000  0.01200  0.03940  0.01010  0.00000  0.00000 0    0.000    0.000  1;% 26
  21   20  100.00 100.00 60 0   0.0000  0.01830  0.08490  0.02160  0.00000  0.00000 0    0.000    0.000  1;% 27
  22   21  100.00 100.00 60 0   0.0000  0.02090  0.09700  0.02460  0.00000  0.00000 0    0.000    0.000  1;% 28
  23   22  100.00 100.00 60 0   0.0000  0.03420  0.15900  0.04040  0.00000  0.00000 0    0.000    0.000  1;% 29
  24   23  100.00 100.00 60 0   0.0000  0.01350  0.04920  0.04980  0.00000  0.00000 0    0.000    0.000  1;% 30
  25   23  100.00 100.00 60 0   0.0000  0.01560  0.08000  0.08640  0.00000  0.00000 0    0.000    0.000  1;% 31
  26   25  100.00 100.00 60 0   0.0000  0.00000  0.03820  0.00000  0.96000  0.00000 0    0.000    0.000  1;% 32
  27   25  100.00 100.00 60 0   0.0000  0.03180  0.16300  0.17640  0.00000  0.00000 0    0.000    0.000  1;% 33
  28   27  100.00 100.00 60 0   0.0000  0.01913  0.08550  0.02160  0.00000  0.00000 0    0.000    0.000  1;% 34
  29   28  100.00 100.00 60 0   0.0000  0.02370  0.09430  0.02380  0.00000  0.00000 0    0.000    0.000  1;% 35
  30   17  100.00 100.00 60 0   0.0000  0.00000  0.03880  0.00000  0.96000  0.00000 0    0.000    0.000  1;% 36
  30    8  100.00 100.00 60 0   0.0000  0.00431  0.05040  0.51400  0.00000  0.00000 0    0.000    0.000  1;% 37
  30   26  100.00 100.00 60 0   0.0000  0.00799  0.08600  0.90800  0.00000  0.00000 0    0.000    0.000  1;% 38
  31   17  100.00 100.00 60 0   0.0000  0.04740  0.15630  0.03990  0.00000  0.00000 0    0.000    0.000  1;% 39
  31   29  100.00 100.00 60 0   0.0000  0.01080  0.03310  0.00830  0.00000  0.00000 0    0.000    0.000  1;% 40
  32   23  100.00 100.00 60 0   0.0000  0.03170  0.11530  0.11730  0.00000  0.00000 0    0.000    0.000  1;% 41
  32   31  100.00 100.00 60 0   0.0000  0.02980  0.09850  0.02510  0.00000  0.00000 0    0.000    0.000  1;% 42
  32   27  100.00 100.00 60 0   0.0000  0.02290  0.07550  0.01926  0.00000  0.00000 0    0.000    0.000  1;% 43
  33   15  100.00 100.00 60 0   0.0000  0.03800  0.12440  0.03194  0.00000  0.00000 0    0.000    0.000  1;% 44
  34   19  100.00 100.00 60 0   0.0000  0.07520  0.24700  0.06320  0.00000  0.00000 0    0.000    0.000  1;% 45
  36   35  100.00 100.00 60 0   0.0000  0.00224  0.01020  0.00268  0.00000  0.00000 0    0.000    0.000  1;% 46
  37   35  100.00 100.00 60 0   0.0000  0.01100  0.04970  0.01318  0.00000  0.00000 0    0.000    0.000  1;% 47
  37   33  100.00 100.00 60 0   0.0000  0.04150  0.14200  0.03660  0.00000  0.00000 0    0.000    0.000  1;% 48
  36   34  100.00 100.00 60 0   0.0000  0.00871  0.02680  0.00568  0.00000  0.00000 0    0.000    0.000  1;% 49
  37   34  100.00 100.00 60 0   0.0000  0.00256  0.00940  0.00984  0.00000  0.00000 0    0.000    0.000  1;% 50
  38   37  100.00 100.00 60 0   0.0000  0.00000  0.03750  0.00000  0.93500  0.00000 0    0.000    0.000  1;% 51
  39   37  100.00 100.00 60 0   0.0000  0.03210  0.10600  0.02700  0.00000  0.00000 0    0.000    0.000  1;% 52
  40   37  100.00 100.00 60 0   0.0000  0.05930  0.16800  0.04200  0.00000  0.00000 0    0.000    0.000  1;% 53
  38   30  100.00 100.00 60 0   0.0000  0.00464  0.05400  0.42200  0.00000  0.00000 0    0.000    0.000  1;% 54
  40   39  100.00 100.00 60 0   0.0000  0.01840  0.06050  0.01552  0.00000  0.00000 0    0.000    0.000  1;% 55
  41   40  100.00 100.00 60 0   0.0000  0.01450  0.04870  0.01222  0.00000  0.00000 0    0.000    0.000  1;% 56
  42   40  100.00 100.00 60 0   0.0000  0.05550  0.18300  0.04660  0.00000  0.00000 0    0.000    0.000  1;% 57
  42   41  100.00 100.00 60 0   0.0000  0.04100  0.13500  0.03440  0.00000  0.00000 0    0.000    0.000  1;% 58
  44   43  100.00 100.00 60 0   0.0000  0.06080  0.24540  0.06068  0.00000  0.00000 0    0.000    0.000  1;% 59
  43   34  100.00 100.00 60 0   0.0000  0.04130  0.16810  0.04226  0.00000  0.00000 0    0.000    0.000  1;% 60
  45   44  100.00 100.00 60 0   0.0000  0.02240  0.09010  0.02240  0.00000  0.00000 0    0.000    0.000  1;% 61
  46   45  100.00 100.00 60 0   0.0000  0.04000  0.13560  0.03320  0.00000  0.00000 0    0.000    0.000  1;% 62
  47   46  100.00 100.00 60 0   0.0000  0.03800  0.12700  0.03160  0.00000  0.00000 0    0.000    0.000  1;% 63
  48   46  100.00 100.00 60 0   0.0000  0.06010  0.18900  0.04720  0.00000  0.00000 0    0.000    0.000  1;% 64
  49   47  100.00 100.00 60 0   0.0000  0.01910  0.06250  0.01604  0.00000  0.00000 0    0.000    0.000  1;% 65
  49   42  100.00 100.00 60 0   0.0000  0.07150  0.32300  0.08600  0.00000  0.00000 0    0.000    0.000  1;% 66
  49   42  100.00 100.00 60 0   0.0000  0.07150  0.32300  0.08600  0.00000  0.00000 0    0.000    0.000  1;% 67
  49   45  100.00 100.00 60 0   0.0000  0.06840  0.18600  0.04440  0.00000  0.00000 0    0.000    0.000  1;% 68
  49   48  100.00 100.00 60 0   0.0000  0.01790  0.05050  0.01258  0.00000  0.00000 0    0.000    0.000  1;% 69
  50   49  100.00 100.00 60 0   0.0000  0.02670  0.07520  0.01874  0.00000  0.00000 0    0.000    0.000  1;% 70
  51   49  100.00 100.00 60 0   0.0000  0.04860  0.13700  0.03420  0.00000  0.00000 0    0.000    0.000  1;% 71
  52   51  100.00 100.00 60 0   0.0000  0.02030  0.05880  0.01396  0.00000  0.00000 0    0.000    0.000  1;% 72
  53   52  100.00 100.00 60 0   0.0000  0.04050  0.16350  0.04058  0.00000  0.00000 0    0.000    0.000  1;% 73
  54   53  100.00 100.00 60 0   0.0000  0.02630  0.12200  0.03100  0.00000  0.00000 0    0.000    0.000  1;% 74
  54   49  100.00 100.00 60 0   0.0000  0.07300  0.28900  0.07380  0.00000  0.00000 0    0.000    0.000  1;% 75
  54   49  100.00 100.00 60 0   0.0000  0.08690  0.29100  0.07300  0.00000  0.00000 0    0.000    0.000  1;% 76
  55   54  100.00 100.00 60 0   0.0000  0.01690  0.07070  0.02020  0.00000  0.00000 0    0.000    0.000  1;% 77
  56   54  100.00 100.00 60 0   0.0000  0.00275  0.00955  0.00732  0.00000  0.00000 0    0.000    0.000  1;% 78
  56   55  100.00 100.00 60 0   0.0000  0.00488  0.01510  0.00374  0.00000  0.00000 0    0.000    0.000  1;% 79
  57   56  100.00 100.00 60 0   0.0000  0.03430  0.09660  0.02420  0.00000  0.00000 0    0.000    0.000  1;% 80
  57   50  100.00 100.00 60 0   0.0000  0.04740  0.13400  0.03320  0.00000  0.00000 0    0.000    0.000  1;% 81
  58   56  100.00 100.00 60 0   0.0000  0.03430  0.09660  0.02420  0.00000  0.00000 0    0.000    0.000  1;% 82
  58   51  100.00 100.00 60 0   0.0000  0.02550  0.07190  0.01788  0.00000  0.00000 0    0.000    0.000  1;% 83
  59   54  100.00 100.00 60 0   0.0000  0.05030  0.22930  0.05980  0.00000  0.00000 0    0.000    0.000  1;% 84
  59   56  100.00 100.00 60 0   0.0000  0.08250  0.25100  0.05690  0.00000  0.00000 0    0.000    0.000  1;% 85
  59   56  100.00 100.00 60 0   0.0000  0.08030  0.23900  0.05360  0.00000  0.00000 0    0.000    0.000  1;% 86
  59   55  100.00 100.00 60 0   0.0000  0.04739  0.21580  0.05646  0.00000  0.00000 0    0.000    0.000  1;% 87
  60   59  100.00 100.00 60 0   0.0000  0.03170  0.14500  0.03760  0.00000  0.00000 0    0.000    0.000  1;% 88
  61   59  100.00 100.00 60 0   0.0000  0.03280  0.15000  0.03880  0.00000  0.00000 0    0.000    0.000  1;% 89
  61   60  100.00 100.00 60 0   0.0000  0.00264  0.01350  0.01456  0.00000  0.00000 0    0.000    0.000  1;% 90
  62   60  100.00 100.00 60 0   0.0000  0.01230  0.05610  0.01468  0.00000  0.00000 0    0.000    0.000  1;% 91
  62   61  100.00 100.00 60 0   0.0000  0.00824  0.03760  0.00980  0.00000  0.00000 0    0.000    0.000  1;% 92
  63   59  100.00 100.00 60 0   0.0000  0.00000  0.03860  0.00000  0.96000  0.00000 0    0.000    0.000  1;% 93
  64   63  100.00 100.00 60 0   0.0000  0.00172  0.02000  0.21600  0.00000  0.00000 0    0.000    0.000  1;% 94
  64   61  100.00 100.00 60 0   0.0000  0.00000  0.02680  0.00000  0.98500  0.00000 0    0.000    0.000  1;% 95
  65   38  100.00 100.00 60 0   0.0000  0.00901  0.09860  1.04600  0.00000  0.00000 0    0.000    0.000  1;% 96
  65   64  100.00 100.00 60 0   0.0000  0.00269  0.03020  0.38000  0.00000  0.00000 0    0.000    0.000  1;% 97
  66   49  100.00 100.00 60 0   0.0000  0.01800  0.09190  0.02480  0.00000  0.00000 0    0.000    0.000  1;% 98
  66   49  100.00 100.00 60 0   0.0000  0.01800  0.09190  0.02480  0.00000  0.00000 0    0.000    0.000  1;% 99
  66   62  100.00 100.00 60 0   0.0000  0.04820  0.21800  0.05780  0.00000  0.00000 0    0.000    0.000  1;% 100
  67   62  100.00 100.00 60 0   0.0000  0.02580  0.11700  0.03100  0.00000  0.00000 0    0.000    0.000  1;% 101
  65   66  100.00 100.00 60 0   0.0000  0.00000  0.03700  0.00000  0.93500  0.00000 0    0.000    0.000  1;% 102
  67   66  100.00 100.00 60 0   0.0000  0.02240  0.10150  0.02682  0.00000  0.00000 0    0.000    0.000  1;% 103
  68   65  100.00 100.00 60 0   0.0000  0.00138  0.01600  0.63800  0.00000  0.00000 0    0.000    0.000  1;% 104
  69   47  100.00 100.00 60 0   0.0000  0.08440  0.27780  0.07092  0.00000  0.00000 0    0.000    0.000  1;% 105
  69   49  100.00 100.00 60 0   0.0000  0.09850  0.32400  0.08280  0.00000  0.00000 0    0.000    0.000  1;% 106
  68   69  100.00 100.00 60 0   0.0000  0.00000  0.03700  0.00000  0.93500  0.00000 0    0.000    0.000  1;% 107
  70   69  100.00 100.00 60 0   0.0000  0.03000  0.12700  0.12200  0.00000  0.00000 0    0.000    0.000  1;% 108
  70   24  100.00 100.00 60 0   0.0000  0.00221  0.41150  0.10198  0.00000  0.00000 0    0.000    0.000  1;% 109
  71   70  100.00 100.00 60 0   0.0000  0.00882  0.03550  0.00878  0.00000  0.00000 0    0.000    0.000  1;% 110
  72   24  100.00 100.00 60 0   0.0000  0.04880  0.19600  0.04880  0.00000  0.00000 0    0.000    0.000  1;% 111
  72   71  100.00 100.00 60 0   0.0000  0.04460  0.18000  0.04444  0.00000  0.00000 0    0.000    0.000  1;% 112
  73   71  100.00 100.00 60 0   0.0000  0.00866  0.04540  0.01178  0.00000  0.00000 0    0.000    0.000  1;% 113
  74   70  100.00 100.00 60 0   0.0000  0.04010  0.13230  0.03368  0.00000  0.00000 0    0.000    0.000  1;% 114
  75   70  100.00 100.00 60 0   0.0000  0.04280  0.14100  0.03600  0.00000  0.00000 0    0.000    0.000  1;% 115
  75   69  100.00 100.00 60 0   0.0000  0.04050  0.12200  0.12400  0.00000  0.00000 0    0.000    0.000  1;% 116
  75   74  100.00 100.00 60 0   0.0000  0.01230  0.04060  0.01034  0.00000  0.00000 0    0.000    0.000  1;% 117
  77   76  100.00 100.00 60 0   0.0000  0.04440  0.14800  0.03680  0.00000  0.00000 0    0.000    0.000  1;% 118
  77   69  100.00 100.00 60 0   0.0000  0.03090  0.10100  0.10380  0.00000  0.00000 0    0.000    0.000  1;% 119
  77   75  100.00 100.00 60 0   0.0000  0.06010  0.19990  0.04978  0.00000  0.00000 0    0.000    0.000  1;% 120
  78   77  100.00 100.00 60 0   0.0000  0.00376  0.01240  0.01264  0.00000  0.00000 0    0.000    0.000  1;% 121
  79   78  100.00 100.00 60 0   0.0000  0.00546  0.02440  0.00648  0.00000  0.00000 0    0.000    0.000  1;% 122
  80   77  100.00 100.00 60 0   0.0000  0.01700  0.04850  0.04720  0.00000  0.00000 0    0.000    0.000  1;% 123
  80   77  100.00 100.00 60 0   0.0000  0.02940  0.10500  0.02280  0.00000  0.00000 0    0.000    0.000  1;% 124
  80   79  100.00 100.00 60 0   0.0000  0.01560  0.07040  0.01870  0.00000  0.00000 0    0.000    0.000  1;% 125
  81   68  100.00 100.00 60 0   0.0000  0.00175  0.02020  0.80800  0.00000  0.00000 0    0.000    0.000  1;% 126
  81   80  100.00 100.00 60 0   0.0000  0.00000  0.03700  0.00000  0.93500  0.00000 0    0.000    0.000  1;% 127
  82   77  100.00 100.00 60 0   0.0000  0.02980  0.08530  0.08174  0.00000  0.00000 0    0.000    0.000  1;% 128
  83   82  100.00 100.00 60 0   0.0000  0.01120  0.03665  0.03796  0.00000  0.00000 0    0.000    0.000  1;% 129
  84   83  100.00 100.00 60 0   0.0000  0.06250  0.13200  0.02580  0.00000  0.00000 0    0.000    0.000  1;% 130
  85   83  100.00 100.00 60 0   0.0000  0.04300  0.14800  0.03480  0.00000  0.00000 0    0.000    0.000  1;% 131
  85   84  100.00 100.00 60 0   0.0000  0.03020  0.06410  0.01234  0.00000  0.00000 0    0.000    0.000  1;% 132
  86   85  100.00 100.00 60 0   0.0000  0.03500  0.12300  0.02760  0.00000  0.00000 0    0.000    0.000  1;% 133
  87   86  100.00 100.00 60 0   0.0000  0.02828  0.20740  0.04450  0.00000  0.00000 0    0.000    0.000  1;% 134
  88   85  100.00 100.00 60 0   0.0000  0.02000  0.10200  0.02760  0.00000  0.00000 0    0.000    0.000  1;% 135
  89   85  100.00 100.00 60 0   0.0000  0.02390  0.17300  0.04700  0.00000  0.00000 0    0.000    0.000  1;% 136
  89   88  100.00 100.00 60 0   0.0000  0.01390  0.07120  0.01934  0.00000  0.00000 0    0.000    0.000  1;% 137
  90   89  100.00 100.00 60 0   0.0000  0.05180  0.18800  0.05280  0.00000  0.00000 0    0.000    0.000  1;% 138
  90   89  100.00 100.00 60 0   0.0000  0.02380  0.09970  0.10600  0.00000  0.00000 0    0.000    0.000  1;% 139
  91   90  100.00 100.00 60 0   0.0000  0.02540  0.08360  0.02140  0.00000  0.00000 0    0.000    0.000  1;% 140
  92   89  100.00 100.00 60 0   0.0000  0.00990  0.05050  0.05480  0.00000  0.00000 0    0.000    0.000  1;% 141
  92   89  100.00 100.00 60 0   0.0000  0.03930  0.15810  0.04140  0.00000  0.00000 0    0.000    0.000  1;% 142
  92   91  100.00 100.00 60 0   0.0000  0.03870  0.12720  0.03268  0.00000  0.00000 0    0.000    0.000  1;% 143
  93   92  100.00 100.00 60 0   0.0000  0.02580  0.08480  0.02180  0.00000  0.00000 0    0.000    0.000  1;% 144
  94   92  100.00 100.00 60 0   0.0000  0.04810  0.15800  0.04060  0.00000  0.00000 0    0.000    0.000  1;% 145
  94   93  100.00 100.00 60 0   0.0000  0.02230  0.07320  0.01876  0.00000  0.00000 0    0.000    0.000  1;% 146
  95   94  100.00 100.00 60 0   0.0000  0.01320  0.04340  0.01110  0.00000  0.00000 0    0.000    0.000  1;% 147
  96   80  100.00 100.00 60 0   0.0000  0.03560  0.18200  0.04940  0.00000  0.00000 0    0.000    0.000  1;% 148
  96   82  100.00 100.00 60 0   0.0000  0.01620  0.05300  0.05440  0.00000  0.00000 0    0.000    0.000  1;% 149
  96   94  100.00 100.00 60 0   0.0000  0.02690  0.08690  0.02300  0.00000  0.00000 0    0.000    0.000  1;% 150
  97   80  100.00 100.00 60 0   0.0000  0.01830  0.09340  0.02540  0.00000  0.00000 0    0.000    0.000  1;% 151
  98   80  100.00 100.00 60 0   0.0000  0.02380  0.10800  0.02860  0.00000  0.00000 0    0.000    0.000  1;% 152
  99   80  100.00 100.00 60 0   0.0000  0.04540  0.20600  0.05460  0.00000  0.00000 0    0.000    0.000  1;% 153
 100   92  100.00 100.00 60 0   0.0000  0.06480  0.29500  0.04720  0.00000  0.00000 0    0.000    0.000  1;% 154
 100   94  100.00 100.00 60 0   0.0000  0.01780  0.05800  0.06040  0.00000  0.00000 0    0.000    0.000  1;% 155
  96   95  100.00 100.00 60 0   0.0000  0.01710  0.05470  0.01474  0.00000  0.00000 0    0.000    0.000  1;% 156
  97   96  100.00 100.00 60 0   0.0000  0.01730  0.08850  0.02400  0.00000  0.00000 0    0.000    0.000  1;% 157
 100   98  100.00 100.00 60 0   0.0000  0.03970  0.17900  0.04760  0.00000  0.00000 0    0.000    0.000  1;% 158
 100   99  100.00 100.00 60 0   0.0000  0.01800  0.08130  0.02160  0.00000  0.00000 0    0.000    0.000  1;% 159
 101  100  100.00 100.00 60 0   0.0000  0.02770  0.12620  0.03280  0.00000  0.00000 0    0.000    0.000  1;% 160
 102   92  100.00 100.00 60 0   0.0000  0.01230  0.05590  0.01464  0.00000  0.00000 0    0.000    0.000  1;% 161
 102  101  100.00 100.00 60 0   0.0000  0.02460  0.11200  0.02940  0.00000  0.00000 0    0.000    0.000  1;% 162
 103  100  100.00 100.00 60 0   0.0000  0.01600  0.05250  0.05360  0.00000  0.00000 0    0.000    0.000  1;% 163
 104  100  100.00 100.00 60 0   0.0000  0.04510  0.20400  0.05410  0.00000  0.00000 0    0.000    0.000  1;% 164
 104  103  100.00 100.00 60 0   0.0000  0.04660  0.15840  0.04070  0.00000  0.00000 0    0.000    0.000  1;% 165
 105  103  100.00 100.00 60 0   0.0000  0.05350  0.16250  0.04080  0.00000  0.00000 0    0.000    0.000  1;% 166
 106  100  100.00 100.00 60 0   0.0000  0.06050  0.22900  0.06200  0.00000  0.00000 0    0.000    0.000  1;% 167
 105  104  100.00 100.00 60 0   0.0000  0.00994  0.03780  0.00986  0.00000  0.00000 0    0.000    0.000  1;% 168
 106  105  100.00 100.00 60 0   0.0000  0.01400  0.05470  0.01434  0.00000  0.00000 0    0.000    0.000  1;% 169
 107  105  100.00 100.00 60 0   0.0000  0.05300  0.18300  0.04720  0.00000  0.00000 0    0.000    0.000  1;% 170
 108  105  100.00 100.00 60 0   0.0000  0.02610  0.07030  0.01844  0.00000  0.00000 0    0.000    0.000  1;% 171
 107  106  100.00 100.00 60 0   0.0000  0.05300  0.18300  0.04720  0.00000  0.00000 0    0.000    0.000  1;% 172
 109  108  100.00 100.00 60 0   0.0000  0.01050  0.02880  0.00760  0.00000  0.00000 0    0.000    0.000  1;% 173
 110  103  100.00 100.00 60 0   0.0000  0.03906  0.18130  0.04610  0.00000  0.00000 0    0.000    0.000  1;% 174
 110  109  100.00 100.00 60 0   0.0000  0.02780  0.07620  0.02020  0.00000  0.00000 0    0.000    0.000  1;% 175
 111  110  100.00 100.00 60 0   0.0000  0.02200  0.07550  0.02000  0.00000  0.00000 0    0.000    0.000  1;% 176
 112  110  100.00 100.00 60 0   0.0000  0.02470  0.06400  0.06200  0.00000  0.00000 0    0.000    0.000  1;% 177
 113   17  100.00 100.00 60 0   0.0000  0.00913  0.03010  0.00768  0.00000  0.00000 0    0.000    0.000  1;% 178
 113   32  100.00 100.00 60 0   0.0000  0.06150  0.20300  0.05180  0.00000  0.00000 0    0.000    0.000  1;% 179
 114   32  100.00 100.00 60 0   0.0000  0.01350  0.06120  0.01628  0.00000  0.00000 0    0.000    0.000  1;% 180
 115   27  100.00 100.00 60 0   0.0000  0.01640  0.07410  0.01972  0.00000  0.00000 0    0.000    0.000  1;% 181
 115  114  100.00 100.00 60 0   0.0000  0.00230  0.01040  0.00276  0.00000  0.00000 0    0.000    0.000  1;% 182
 116   68  100.00 100.00 60 0   0.0000  0.00034  0.00405  0.16400  0.00000  0.00000 0    0.000    0.000  1;% 183
 117   12  100.00 100.00 60 0   0.0000  0.03290  0.14000  0.03580  0.00000  0.00000 0    0.000    0.000  1;% 184
 118   75  100.00 100.00 60 0   0.0000  0.01450  0.04810  0.01198  0.00000  0.00000 0    0.000    0.000  1;% 185
 118   76  100.00 100.00 60 0   0.0000  0.01640  0.05440  0.01356  0.00000  0.00000 0    0.000    0.000  1;% 186
  ];

Bus.names = { ...
  'bus-1   100'; 'bus-2   100'; 'bus-3   100'; 'bus-4   100'; 'bus-5   100';
  'bus-6   100'; 'bus-7   100'; 'bus-8   100'; 'bus-9   100'; 'bus-10  100';
  'bus-11  100'; 'bus-12  100'; 'bus-13  100'; 'bus-14  100'; 'bus-15  100';
  'bus-16  100'; 'bus-17  100'; 'bus-18  100'; 'bus-19  100'; 'bus-20  100';
  'bus-21  100'; 'bus-22  100'; 'bus-23  100'; 'bus-24  100'; 'bus-25  100';
  'bus-26  100'; 'bus-27  100'; 'bus-28  100'; 'bus-29  100'; 'bus-30  100';
  'bus-31  100'; 'bus-32  100'; 'bus-33  100'; 'bus-34  100'; 'bus-35  100';
  'bus-36  100'; 'bus-37  100'; 'bus-38  100'; 'bus-39  100'; 'bus-40  100';
  'bus-41  100'; 'bus-42  100'; 'bus-43  100'; 'bus-44  100'; 'bus-45  100';
  'bus-46  100'; 'bus-47  100'; 'bus-48  100'; 'bus-49  100'; 'bus-50  100';
  'bus-51  100'; 'bus-52  100'; 'bus-53  100'; 'bus-54  100'; 'bus-55  100';
  'bus-56  100'; 'bus-57  100'; 'bus-58  100'; 'bus-59  100'; 'bus-60  100';
  'bus-61  100'; 'bus-62  100'; 'bus-63  100'; 'bus-64  100'; 'bus-65  100';
  'bus-66  100'; 'bus-67  100'; 'bus-68  100'; 'bus-69  100'; 'bus-70  100';
  'bus-71  100'; 'bus-72  100'; 'bus-73  100'; 'bus-74  100'; 'bus-75  100';
  'bus-76  100'; 'bus-77  100'; 'bus-78  100'; 'bus-79  100'; 'bus-80  100';
  'bus-81  100'; 'bus-82  100'; 'bus-83  100'; 'bus-84  100'; 'bus-85  100';
  'bus-86  100'; 'bus-87  100'; 'bus-88  100'; 'bus-89  100'; 'bus-90  100';
  'bus-91  100'; 'bus-92  100'; 'bus-93  100'; 'bus-94  100'; 'bus-95  100';
  'bus-96  100'; 'bus-97  100'; 'bus-98  100'; 'bus-99  100'; 'bus-100 100';
  'bus-101 100'; 'bus-102 100'; 'bus-103 100'; 'bus-104 100'; 'bus-105 100';
  'bus-106 100'; 'bus-107 100'; 'bus-108 100'; 'bus-109 100'; 'bus-110 100';
  'bus-111 100'; 'bus-112 100'; 'bus-113 100'; 'bus-114 100'; 'bus-115 100';
  'bus-116 100'; 'bus-117 100'; 'bus-118 100'};
