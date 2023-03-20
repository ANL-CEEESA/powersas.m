% 03/16/23 File data originated from Matpower data file
% 

Bus.con = [ ...
      1      345    1.031  -0.2141    1    1;
      2      345    1.002  -0.2502    1    1;
      3      345   0.9535  -0.3343    1    1;
      4      345    0.931  -0.3363    1    2;
      5      345   0.9174  -0.3049    1    2;
      6      345   0.9186  -0.2919    1    2;
      7      345   0.8626  -0.3252    1    2;
      8      345     0.88  -0.3322    1    2;
      9      345   0.9788  -0.2449    1    2;
     10      345   0.9562  -0.2565    1    2;
     11      345   0.9429  -0.2703    1    2;
     12      138   0.9449  -0.3003    1    2;
     13      345   0.9523  -0.2722    1    2;
     14      345   0.9465  -0.3029    2    3;
     15      345   0.9565   -0.314    2    3;
     16      345   0.9767  -0.2878    2    3;
     17      345   0.9749  -0.3076    1    1;
     18      345   0.9652  -0.3262    1    1;
     19      345     1.03  -0.2038    2    4;
     20      138   0.9798  -0.2297    2    4;
     21      345    0.993  -0.2428    2    3;
     22      345    1.029  -0.1615    2    3;
     23      345    1.023  -0.1651    2    4;
     24      345   0.9871  -0.2857    2    3;
     25      345    1.021  -0.2356    1    1;
     26      345     1.01  -0.2661    1    1;
     27      345   0.9876  -0.3071    1    1;
     28      345     1.02  -0.2014    2    3;
     29      345    1.024  -0.1509    2    3;
     30     34.5    1.005  -0.2041    1    1;
     31     34.5   0.9064  -0.1583    1    2;
     32       21   0.9752  -0.1068    1    2;
     33       21   0.9972   -0.112    2    4;
     34     15.5    1.012  -0.1387    2    4;
     35     15.5    1.049 -0.07315    2    3;
     36     12.5    1.063 -0.02552    2    4;
     37     12.5    1.028  -0.1139    1    1;
     38     34.5    1.026 -0.02585    2    3;
     39      345     1.03  -0.1913    1    2;
   ];

SW.con = [ ...
     39     1199      345     1.03  -0.1913   0.4794  -0.1445     1.05     0.95   0.6012 1 1  1;
   ];

PV.con = [ ...
  30      275     34.5   0.9091    1.048    0.557  -0.2112     1.05     0.95 1  1;
  31      836     34.5   0.6853     1.04   0.5141  -0.1467     1.05     0.95 1  1;
  32    843.7       21   0.7704   0.9831   0.5295  -0.2136     1.05     0.95 1  1;
  33     1175       21    0.538   0.9972   0.4666   -0.182     1.05     0.95 1  1;
  34     1080     15.5   0.4703    1.012   0.5661  -0.1741     1.05     0.95 1  1;
  35     1086     15.5   0.5987    1.049   0.5469  -0.2164     1.05     0.95 1  1;
  36     1025     12.5   0.5462    1.063   0.5544   -0.243     1.05     0.95 1  1;
  37    970.2     12.5   0.5566    1.028   0.4571  -0.2228     1.05     0.95 1  1;
  38     1684     34.5   0.4928    1.026   0.4957  -0.2119     1.05     0.95 1  1;
  10      900      345   0.1111    1.004   0.3287  -0.3287     1.05     0.95 1  1;
  20      900      138   0.1667    1.004   0.3287  -0.3287     1.05     0.95 1  1;
   2      900      345   0.2222    1.004   0.3287  -0.3287     1.05     0.95 1  1;
  25      900      345   0.2778    1.004   0.3287  -0.3287     1.05     0.95 1  1;
   ];

PQ.con = [ ...
   3 100      345        6      2.5     1.05     0.95 0 1;
   4 100      345      4.5     1.84     1.05     0.95 0 1;
   7 100      345    2.338      8.4     1.05     0.95 0 1;
   8 100      345     5.22     1.76     1.05     0.95 0 1;
  12 100      138      1.2      0.3     1.05     0.95 0 1;
  15 100      345      3.2     1.53     1.05     0.95 0 1;
  16 100      345    3.294     3.23     1.05     0.95 0 1;
  18 100      345     1.58      0.3     1.05     0.95 0 1;
  20 100      138      6.8     1.03     1.05     0.95 0 1;
  21 100      345     2.74     1.15     1.05     0.95 0 1;
  23 100      345    2.475    0.846     1.05     0.95 0 1;
  24 100      345    3.086   -0.922     1.05     0.95 0 1;
  25 100      345     2.24    0.472     1.05     0.95 0 1;
  26 100      345     1.39     0.17     1.05     0.95 0 1;
  27 100      345     2.81    0.755     1.05     0.95 0 1;
  28 100      345     2.06    0.276     1.05     0.95 0 1;
  29 100      345    2.835    1.269     1.05     0.95 0 1;
  31 100     34.5      0.8      0.4     1.05     0.95 0 1;
  39 100      345        4      2.5     1.05     0.95 0 1;
   1 100      345        0        0     1.05     0.95 0 1;
   5 100      345        0        0     1.05     0.95 0 1;
   6 100      345        0        0     1.05     0.95 0 1;
   9 100      345        0        0     1.05     0.95 0 1;
  11 100      345        0        0     1.05     0.95 0 1;
  13 100      345        0        0     1.05     0.95 0 1;
  14 100      345        0        0     1.05     0.95 0 1;
  17 100      345        0        0     1.05     0.95 0 1;
  19 100      345        0        0     1.05     0.95 0 1;
  22 100      345        0        0     1.05     0.95 0 1;
   ];

Shunt.con = [ ...
   4 100      345 60        0        1 1;
   5 100      345 60        0        2 1;
   ];

Line.con = [ ...
   1    2 100      345 60 0         0   0.0035   0.0411   0.6987        0        0     1.75     1.92     1.92  1;
   1   39 100      345 60 0         0    0.001    0.025     0.75        0        0     1.25     1.37     1.37  1;
   2    3 100      345 60 0         0   0.0013   0.0151   0.2572        0        0      6.5     7.15     7.15  1;
   2   25 100      345 60 0         0    0.007   0.0086    0.146        0        0      2.6     2.86     2.86  1;
   3    4 100      345 60 0         0   0.0013   0.0213   0.2214        0        0     1.25     1.37     1.37  1;
   3   18 100      345 60 0         0   0.0011   0.0133   0.2138        0        0     1.25     1.37     1.37  1;
   4    5 100      345 60 0         0   0.0008   0.0128   0.1342        0        0        3      3.3      3.3  1;
   4   14 100      345 60 0         0   0.0008   0.0129   0.1382        0        0        3      3.3      3.3  1;
   5    6 100      345 60 0         0   0.0002   0.0026   0.0434        0        0     4.62     5.08     5.08  1;
   5    8 100      345 60 0         0   0.0008   0.0112   0.1476        0        0     4.62     5.08     5.08  1;
   6    7 100      345 60 0         0   0.0006   0.0092    0.113        0        0        7      7.7      7.7  1;
   6   11 100      345 60 0         0   0.0007   0.0082   0.1389        0        0     4.66     5.12     5.12  1;
   7    8 100      345 60 0         0   0.0004   0.0046    0.078        0        0     4.62     5.08     5.08  1;
   8    9 100      345 60 0         0   0.0023   0.0363   0.3804        0        0     4.66     5.12     5.12  1;
   9   39 100      345 60 0         0    0.001    0.025      1.2        0        0        4      4.4      4.4  1;
  10   11 100      345 60 0         0   0.0004   0.0043   0.0729        0        0     4.62     5.08     5.08  1;
  10   13 100      345 60 0         0   0.0004   0.0043   0.0729        0        0     4.29     4.71     4.71  1;
  13   14 100      345 60 0         0   0.0009   0.0101   0.1723        0        0      3.8     4.18     4.18  1;
  14   15 100      345 60 0         0   0.0018   0.0217    0.366        0        0     1.25     1.37     1.37  1;
  15   16 100      345 60 0         0   0.0009   0.0094    0.171        0        0     4.29     4.71     4.71  1;
  16   17 100      345 60 0         0   0.0007   0.0089   0.1342        0        0      2.6     2.86     2.86  1;
  16   19 100      345 60 0         0   0.0016   0.0195    0.304        0        0     5.52     6.07     6.07  1;
  16   21 100      345 60 0         0   0.0008   0.0135   0.2548        0        0     4.29     4.71     4.71  1;
  16   24 100      345 60 0         0   0.0003   0.0059    0.068        0        0      2.3     2.53     2.53  1;
  17   18 100      345 60 0         0   0.0007   0.0082   0.1319        0        0        3      3.3      3.3  1;
  17   27 100      345 60 0         0   0.0013   0.0173   0.3216        0        0     1.25     1.37     1.37  1;
  21   22 100      345 60 0         0   0.0008    0.014   0.2565        0        0      6.5     7.15     7.15  1;
  22   23 100      345 60 0         0   0.0006   0.0096   0.1846        0        0     1.25     1.37     1.37  1;
  23   24 100      345 60 0         0   0.0022    0.035    0.361        0        0     4.66     5.12     5.12  1;
  25   26 100      345 60 0         0   0.0032   0.0323    0.513        0        0     1.25     1.37     1.37  1;
  26   27 100      345 60 0         0   0.0014   0.0147   0.2396        0        0      3.8     4.18     4.18  1;
  26   28 100      345 60 0         0   0.0043   0.0474   0.7802        0        0     1.75     1.92     1.92  1;
  26   29 100      345 60 0         0   0.0057   0.0625    1.029        0        0      2.3     2.53     2.53  1;
  28   29 100      345 60 0         0   0.0014   0.0151    0.249        0        0        4      4.4      4.4  1;
   2   30 100      345 60 0        10        0   0.0181        0    1.025        0      3.8     4.18     4.18  1;
  31    6 100     34.5 60 0       0.1        0    0.025        0      0.9        0        7      7.7      7.7  1;
  10   32 100      345 60 0     16.43        0     0.02        0     1.07        0     8.16     8.97     8.97  1;
  12   11 100      138 60 0       0.4   0.0016   0.0435        0    1.006        0     1.25     1.37     1.37  1;
  12   13 100      138 60 0       0.4   0.0016   0.0435        0    1.006        0     1.25     1.37     1.37  1;
  19   20 100      345 60 0       2.5   0.0007   0.0138        0     1.06        0      2.3     2.53     2.53  1;
  19   33 100      345 60 0     16.43   0.0007   0.0142        0     1.07        0        7      7.7      7.7  1;
  20   34 100      138 60 0     8.903   0.0009    0.018        0    1.009        0      5.6     6.16     6.16  1;
  22   35 100      345 60 0     22.26        0   0.0143        0    1.025        0     7.15     7.86     7.86  1;
  23   36 100      345 60 0      27.6   0.0005   0.0272        0        1        0      5.6     6.16     6.16  1;
  25   37 100      345 60 0      27.6   0.0006   0.0232        0    1.025        0     5.52     6.07     6.07  1;
  29   38 100      345 60 0        10   0.0008   0.0156        0    1.025        0     8.37      9.2      9.2  1;
   ];

Bus.names = {...
      'BUS1'; 'BUS2'; 'LOAD3'; 'LOAD4'; 'BUS5'; 
      'BUS6'; 'LOAD7'; 'LOAD8'; 'BUS9'; 'BUS10'; 
      'BUS11'; 'LOAD12'; 'BUS13'; 'BUS14'; 'LOAD15'; 
      'LOAD16'; 'BUS17'; 'LOAD18'; 'BUS19'; 'LOAD20'; 
      'LOAD21'; 'BUS22'; 'LOAD23'; 'LOAD24'; 'LOAD25'; 
      'LOAD26'; 'LOAD27'; 'LOAD28'; 'LOAD29'; 'GEN30'; 
      'GEN31'; 'GEN32'; 'GEN33'; 'GEN34'; 'GEN35'; 
      'GEN36'; 'GEN37'; 'GEN38'; 'GEN39'};

