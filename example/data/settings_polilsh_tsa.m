% Event definition for modified Polish 2383-bus system transient stability analysis
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
% Simulation data
nlvl = 15;
taylorN =  4;
alphaTol =1.0000e-04;
diffTol =1.0000e-05;
diffTolMax =1.0000e-02;
method = 0;
diffTolCtrl =5.0000e-05;

% COLUMNS:
% 1 evtId, 2 startTime, 3 endTime, 4 type, 5 num, 6 method, 7 dt
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
% Method: X.YZ
% Method X: Differential eq. 
%           0:HE, 1:Modified Euler, 2:RK-4, 3: Trapzoidal.
% Method Y: Algebraic eq. 
%           0:HE, 1:N-R.
% Method Z: Step control 
%           0:fixed step, 1:adaptive step.
eventList=[...
   1    0.0000   0.0000    0    1 0.0  0.0000
%   2    0.2000  0.0000    1    1 0.0  0.0000
   1.1  0.5000  0.0000    6    1 0.0  0.0000
   1.2  0.7500  0.0000    6    2 0.0  0.0000
   1.3  1.5000  0.0000    6    3 0.0  0.0000
   1.4  1.9500  0.0000    6    4 0.0  0.0000
%    2   0.0000  200.0000   50    1 0.0  0.0000
   3   10.00    0.0000   99    0 0.0  0.0000
];

% Blackstart data
bsSyn=[...
 ];

bsBus=[...
 ];

bsInd=[...
 ];

% Plain start data
Efstd=1.3;

% Line event data
evtLine=[...
   1    1    1
 ];

evtLineSpec=[...
  11    0    0  0.0000  0.0000
 ];

% Static load event data
evtZip=[...
 ];

evtZipSpec=[...
 ];

evtZipSpec2=[...
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
    1   1   2
    2   3   4
    3   5   5
    4   6   6
 ];

evtFaultSpec=[...
     114,  0.00, 0, 0.02,  0;
      74,  0.00, 0, 0.02,  0;
     114,  0.00, 0, 0.02,  1;
      74,  0.00, 0, 0.02,  1;
    1674,  0.00, 0, 0.1,  0;
    1674,  0.00, 0, 0.1,  1;
 ];

% Dynamic simulation event data
evtDyn=[...
   1    0    0    0    0    0    0    1 12    0    0    0    0    0    0    0    0    0    0    0    0    0    0 
 ];

evtDynPQ=[...
 ];

evtDynPV=[...
 ];

evtDynInd=[...
 ];

evtDynZip=[...
 126  2.0
 188  2.0
 189  2.0
 190  2.0
 191  2.0
 193  2.0
 195  2.0
 197  2.0
 199  2.0
 200  2.0
 201  2.0
 203  2.0
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

