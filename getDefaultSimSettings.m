function simSettings=getDefaultSimSettings()
    % Simulation data
% simSettings.nlvl = 15;
% simSettings.taylorN =  4;
% simSettings.alphaTol =1.0000e-04;
% simSettings.diffTol =1.0000e-06;
% simSettings.diffTolMax =1.0000e-03;
% simSettings.method = 0;
% simSettings.diffTolCtrl =5.0000e-05;

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
simSettings.eventList=[...
   1    0.0000   0.0000    0    1 0.0  0.0000
];

% Blackstart data
simSettings.bsSyn=[...
 ];

simSettings.bsBus=[...
 ];

simSettings.bsInd=[...
 ];

% Plain start data
simSettings.Efstd=1.2;

% Line event data
simSettings.evtLine=[];

simSettings.evtLineSpec=[];

% Static load event data
simSettings.evtZip=[];

simSettings.evtZipSpec=[];

simSettings.evtZipSpec2=[];

% Motor load event data
simSettings.evtInd=[];

simSettings.evtIndSpec=[];

% Syncronous generator event data
simSettings.evtSyn=[];

simSettings.evtSynSpec=[...
 ];

% Fault event data
simSettings.evtFault=[];

simSettings.evtFaultSpec=[];

% Dynamic simulation event data
simSettings.evtDyn=[];

simSettings.evtDynPQ=[...
 ];

simSettings.evtDynPV=[...
 ];

simSettings.evtDynInd=[...
 ];

simSettings.evtDynZip=[];

simSettings.evtDynSh=[];

simSettings.evtDynZipRamp=[];

simSettings.evtDynTmech=[];

simSettings.evtDynPm=[];

simSettings.evtDynEf=[];

simSettings.evtDynVref=[];

simSettings.evtDynEq1=[];

end