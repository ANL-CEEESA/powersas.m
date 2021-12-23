% Event definition for 39-bus system
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
nlvl = 15;
taylorN =  4;
alphaTol =1.0000e-03;
diffTol =1.0000e-06;
diffTolMax =1.0000e-02;
method=0;
diffTolCtrl=2e-6;

% evtId, startTime, endTime, type, num
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
% Method X: Differential eq. 
%           0:HE, 1:Modified Euler, 2:RK-4, 3: TRAP.
% Method Y: Algebraic eq. 0:HE, 1:N-R.
eventList=[...
   1    0.0000  0.0000    0    1 
   1.1  1.5000  0.0000    6    1 
   1.2  1.6700  0.0000    6    2 
%    1.3  1.5000  4.0000   50    1 
   2    5.000  0.0000   99    0 
];

% Blackstart data
bsSyn=[...
 ];

bsBus=[...
 ];

bsInd=[...
 ];

% Plain start data
Efstd=[...
    1.1922
    1.2157
    1.2140
    1.2213
    1.2137
    1.2444
    1.2607
    1.2887
    1.3101
    1.3264 
];


% Line event data
evtLine=[...
 ];

evtLineSpec=[...
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
    1   1   1
    2   2   2
 ];

evtFaultSpec=[...
    14, 0.00, 0, 0.013,  0;
    14, 0.00, 0, 0.013,  1;
 ];

% Dynamic simulation event data
evtDyn=[...
   1    0    0    0    0    0    0    0    0    0    0    0    0    1    2    1    2    0    0    0    0    0    0 
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
   1  0.3459
   2  -0.3459
 ];

evtDynPm=[...
   1  0.3459
   2  -0.3459
 ];

evtDynEf=[...
 ];

evtDynVref=[...
 ];

evtDynEq1=[...
 ];

