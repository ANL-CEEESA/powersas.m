% Event definition for 4-bus system
% Author: Rui Yao <ruiyao@ieee.org>
%
% Copyright (C) 2021, UChicago Argonne, LLC. All rights reserved.
%
% OPEN SOURCE LICENSE
% 
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% 
% 
% ******************************************************************************************************
% DISCLAIMER
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
% WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY 
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
% OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% ***************************************************************************************************
%
% 

nlvl = 14;
taylorN =  4;
alphaTol =1.0000e-04;
diffTol =1.0000e-05;
diffTolMax =1.0000e-02;
method = 0;
diffTolCtrl =5.0000e-05;

% evtId, startTime, endTime, type, num, method, dt
% Type:
%   0: Black start / Full start (steady state)
%   1: Add line
%   2: Add static load
%   3: Add Motor load
%   4: Add syn gen
%  50: Dyn simulation
%   6: Fault
%   7: Cut line
%   8: Cut static load
%   9: Cut Motor load
%  10: Cut syn gen
%  99: End of simulation
% Method: X.Y
% Method X: Differential eq. 0:HE, 1:Modified Euler, 2:RK-4(not active).
% Method Y: Algebraic eq. 0:HE, 1:N-R.
eventList=[...
   1    0.0000  0.0000    0    1 0.0  0.0100
   6   60.0000  0.0000    2    1 0.0  0.0100
   7   90.0000  0.0000    2    2 0.0  0.0100
   8  120.0000  0.0000    2    3 0.0  0.0100
   9  150.0000  0.0000    2    4 0.0  0.0100
  10  180.0000  0.0000    2    1 0.0  0.0100
  11  210.0000  0.0000    2    2 0.0  0.0100
  12  240.0000  0.0000    2    3 0.0  0.0100
  13  270.0000  0.0000    2    4 0.0  0.0100
  18  500.0000  0.0000   99    0 0.0  0.0100
];

% Blackstart data
bsSyn=[...
];

bsBus=[...
];

bsInd=[...
];

% Plain start data
Efstd= 1.2000;

% Line event data
evtLine=[...
];

evtLineSpec=[...
];

% Static load event data
evtZip=[...
   1    1    1    1
   2    1    2    2
   3    1    3    3
   4    1    4    4
%    5    1    8    8
%    6    1    9   12
%    7    1   13   16
%    8    1   17   20
%    9    1   21   24
];

evtZipSpec=[...
 ];

evtZipSpec2=[...
   3 100.0000 100.0000 60.0000  0.0648  0.0648  0.0648  0.0359  0.0359  0.0359    0    1
   2 100.0000 100.0000 60.0000  0.0648  0.0648  0.0648  0.0359  0.0359  0.0359    0    1
   3 100.0000 100.0000 60.0000 -0.0648 -0.0648 -0.0648 -0.0359 -0.0359 -0.0359    0    1
   2 100.0000 100.0000 60.0000 -0.0648 -0.0648 -0.0648 -0.0359 -0.0359 -0.0359    0    1
];

% Motor load event data
evtInd=[...
];

evtIndSpec=[...
];

% Syncronous generator event data
evtSyn=[...
];

evtSynSpec=[...
];

% Fault event data
evtFault=[...
 ];

evtFaultSpec=[...
 ];

% Dynamic simulation event data
% evtId 2pqSt pqEnd 4pvSt pvEnd 6indSt indEnd 8zipSt zipEnd 10shSt shEnd
% 12zipRampSt zipRampEnd 14TmechSt TmechEnd 16PmSt PmEnd 18EfSt EfEnd 
% 20VrefSt VrefEnd 22Eq1St Eq1End
%       2         4         6         8        10        12        14        16        18        20        22
evtDyn=[...
];

evtDynPQ=[...
 ];

evtDynPV=[...
 ];

evtDynInd=[...
 ];

evtDynZip=[...
 ];

evtDynSh=[...
 ];

evtDynZipRamp=[...
 ];

evtDynTmech=[...
];

evtDynPm=[...
];

evtDynEf=[...
 ];

evtDynVref=[...
 ];

evtDynEq1=[...
 ];

