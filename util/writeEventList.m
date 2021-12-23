function writeEventList(outFileName,...
    nlvl,taylorN,alphaTol,diffTol,diffTolMax,method,diffTolCtrl,eventList,...
    bsSyn,bsBus,bsInd,Efstd,evtLine,evtLineSpec,evtZip,evtZipSpec,evtZipSpec2,evtInd,evtIndSpec,evtSyn,evtSynSpec,evtFault,evtFaultSpec,...
    evtDyn,evtDynPQ,evtDynPV,evtDynInd,evtDynZip,evtDynSh,evtDynZipRamp,evtDynTmech,evtDynPm,evtDynEf,evtDynVref,evtDynEq1,intervalEndTimes)
% function for writing event definitions
%
% FUNCTION writeRestorationPlan
%
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

EvtType=getSimEventType();
MethodType=getSimMethodType();
SimSet=getSimSettingDefault();

fid = fopen(outFileName, 'w');

fprintf(fid,'%% Simulation data\r\n');
fprintf(fid,'nlvl =%3d;\r\n',nlvl);
fprintf(fid,'taylorN =%3d;\r\n',taylorN);
fprintf(fid,'alphaTol =%.4e;\r\n',alphaTol);
fprintf(fid,'diffTol =%.4e;\r\n',diffTol);
fprintf(fid,'diffTolMax =%.4e;\r\n',diffTolMax);
fprintf(fid,'method =%2d;\r\n',method);
fprintf(fid,'diffTolCtrl =%.4e;\r\n',diffTolCtrl);
fprintf(fid,'\r\n');

fprintf(fid,'%% evtId, startTime, endTime, type, num, method, dt\r\n');
fprintf(fid,'%% Type:\r\n');
fprintf(fid,'%% %3d: Black start / Full start (steady state)\r\n',EvtType.START);
fprintf(fid,'%% %3d: Add line\r\n',EvtType.ADD_LINE);
fprintf(fid,'%% %3d: Add static load\r\n',EvtType.ADD_LOAD);
fprintf(fid,'%% %3d: Add Motor load\r\n',EvtType.ADD_IND);
fprintf(fid,'%% %3d: Add syn gen\r\n',EvtType.ADD_SYN);
fprintf(fid,'%% %3d: Dyn simulation\r\n',EvtType.DYN_SIM);
fprintf(fid,'%% %3d: Fault\r\n',EvtType.FAULT);
fprintf(fid,'%% %3d: Cut line\r\n',EvtType.CUT_LINE);
fprintf(fid,'%% %3d: Cut static load\r\n',EvtType.CUT_LOAD);
fprintf(fid,'%% %3d: Cut Motor load\r\n',EvtType.CUT_IND);
fprintf(fid,'%% %3d: Cut syn gen\r\n',EvtType.CUT_SYN);
fprintf(fid,'%% %3d: End of simulation\r\n',EvtType.END);
fprintf(fid,'%% Method: X.Y\r\n');
fprintf(fid,'%% Method X: Differential eq. 0:HE, 1:Modified Euler, 2:RK-4(not active).\r\n');
fprintf(fid,'%% Method Y: Algebraic eq. 0:HE, 1:N-R.\r\n');

fprintf(fid,'eventList=[...\r\n');
fprintf(fid,'%6d %9.3f %9.3f %6d %6d %3.1f %9.3f\r\n',eventList');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% Blackstart data\r\n');
fprintf(fid,'%% synId, Ef, Pm, pShare\r\n');
fprintf(fid,'bsSyn=[...\r\n');
fprintf(fid,'%6d %9.6f %9.6f %9.6f\r\n',bsSyn');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% busId, Pz, Pi, Ps, Qz, Qi, Qse\r\n');
fprintf(fid,'bsBus=[...\r\n');
fprintf(fid,'%6d %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f\r\n',bsBus');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% motorId, Tm\r\n');
fprintf(fid,'bsInd=[...\r\n');
fprintf(fid,'%6d %9.6f \r\n',bsInd');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% Plain start data\r\n');
fprintf(fid,'Efstd=%9.6f;\r\n',Efstd);
fprintf(fid,'\r\n');

fprintf(fid,'%% Line event data\r\n');
fprintf(fid,'%% evtId   start   end\r\n');
fprintf(fid,'evtLine=[...\r\n');
fprintf(fid,'%6d %6d %6d\r\n',evtLine');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% Line events spec\r\n');
fprintf(fid,'%% lineId evtType(0-add, 1-trip, 2-fault, 3-resume) pos r-fault x-fault\r\n');
fprintf(fid,'evtLineSpec=[...\r\n');
fprintf(fid,'%6d %6d %6d %9.6f %9.6f %% %6d\r\n',[evtLineSpec,(1:size(evtLineSpec,1))']');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% Static load event data\r\n');
fprintf(fid,'%% evtId   start   end\r\n');
fprintf(fid,'evtZip=[...\r\n');
fprintf(fid,'%6d %6d %6d %6d\r\n',evtZip');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% zip event specs\r\n');
fprintf(fid,'%% zipId type(0-add, 1-trip)\r\n');
fprintf(fid,'evtZipSpec=[...\r\n');
fprintf(fid,'%6d %6d %% %6d\r\n',[evtZipSpec,(1:size(evtZipSpec,1))']');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% zip event specs 2\r\n');
fprintf(fid,'%% Direct add\r\n');
fprintf(fid,'evtZipSpec2=[...\r\n');
fprintf(fid,'%%  Bus          S         V         f         g        Ip         P         b        Iq         Q\r\n');
fprintf(fid,'%6d %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %6d %6d %% %6d\r\n',[evtZipSpec2,(1:size(evtZipSpec2,1))']');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% Motor load event data\r\n');
fprintf(fid,'%% evtId   start   end\r\n');
fprintf(fid,'evtInd=[...\r\n');
fprintf(fid,'%6d %6d %6d\r\n',evtInd');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% ind event specs\r\n');
fprintf(fid,'%% indId type(0-add, 1-adjust, 2-trip) Tm s0\r\n');
fprintf(fid,'evtIndSpec=[...\r\n');
fprintf(fid,'%6d %6d %9.6f %9.6f %% %6d\r\n',[evtIndSpec,(1:size(evtIndSpec,1))']');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% Syncronous generator event data\r\n');
fprintf(fid,'%% evtId   start   end\r\n');
fprintf(fid,'evtSyn=[...\r\n');
fprintf(fid,'%6d %6d %6d\r\n',evtSyn');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% syn event specs\r\n');
fprintf(fid,'%% synId type(0-add, 1-adjust, 2-trip) d0 Pm0 Ef0\r\n');
fprintf(fid,'evtSynSpec=[...\r\n');
fprintf(fid,'%6d %6d %9.6f %9.6f %9.6f %% %6d\r\n',[evtSynSpec,(1:size(evtSynSpec,1))']');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% Fault event data\r\n');
fprintf(fid,'evtFault=[...\r\n');
fprintf(fid,'%6d %6d %6d\r\n',evtFault');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'evtFaultSpec=[...\r\n');
fprintf(fid,'%6d %9.6f %9.6f %9.6f %6d %% %6d\r\n',[evtFaultSpec,(1:size(evtFaultSpec,1))']');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% Dynamic simulation event data\r\n');
fprintf(fid,'%% evtId 2pqSt pqEnd 4pvSt pvEnd 6indSt indEnd 8zipSt zipEnd 10shSt shEnd\r\n');
fprintf(fid,'%% 12zipRampSt zipRampEnd 14TmechSt TmechEnd 16PmSt PmEnd 18EfSt EfEnd \r\n');
fprintf(fid,'%% 20VrefSt VrefEnd 22Eq1St Eq1End\r\n');
fprintf(fid,'%%           2             4             6             8            10            12            14            16            18            20            22\r\n');
fprintf(fid,'evtDyn=[...\r\n');
fprintf(fid,'%6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \r\n',evtDyn');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% busId pIncr qIncr\r\n');
fprintf(fid,'evtDynPQ=[...\r\n');
fprintf(fid,'%6d %9.6f %9.6f %% %6d\r\n',[evtDynPQ,(1:size(evtDynPQ,1))']');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% busId pIncr\r\n');
fprintf(fid,'evtDynPV=[...\r\n');
fprintf(fid,'%6d %9.6f %% %6d\r\n',[evtDynPV,(1:size(evtDynPV,1))']');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% indId TmRamp\r\n');
fprintf(fid,'evtDynInd=[...\r\n');
fprintf(fid,'%6d %9.6f %% %6d\r\n',[evtDynInd,(1:size(evtDynInd,1))']');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% busId zipRamp\r\n');
fprintf(fid,'evtDynZip=[...\r\n');
fprintf(fid,'%6d %9.6f %% %6d\r\n',[evtDynZip,(1:size(evtDynZip,1))']');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% busId shRamp\r\n');
fprintf(fid,'evtDynSh=[...\r\n');
fprintf(fid,'%6d %9.6f %% %6d\r\n',[evtDynSh,(1:size(evtDynSh,1))']');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%%  Bus          S         V         f         g        Ip         P         b        Iq         Q\r\n');
fprintf(fid,'evtDynZipRamp=[...\r\n');
fprintf(fid,'%6d %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %6d %6d %% %6d\r\n',[evtDynZipRamp,(1:size(evtDynZipRamp,1))']');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% synId TmechRamp\r\n');
fprintf(fid,'evtDynTmech=[...\r\n');
fprintf(fid,'%6d %9.6f %% %6d\r\n',[evtDynTmech,(1:size(evtDynTmech,1))']');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% synId PmRamp\r\n');
fprintf(fid,'evtDynPm=[...\r\n');
fprintf(fid,'%6d %9.6f %% %6d\r\n',[evtDynPm,(1:size(evtDynPm,1))']');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% synId Ramp\r\n');
fprintf(fid,'evtDynEf=[...\r\n');
fprintf(fid,'%6d %9.6f %% %6d\r\n',[evtDynEf,(1:size(evtDynEf,1))']');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% synId Ramp\r\n');
fprintf(fid,'evtDynVref=[...\r\n');
fprintf(fid,'%6d %9.6f %% %6d\r\n',[evtDynVref,(1:size(evtDynVref,1))']');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'%% synId Ramp\r\n');
fprintf(fid,'evtDynEq1=[...\r\n');
fprintf(fid,'%6d %9.6f %% %6d\r\n',[evtDynEq1,(1:size(evtDynEq1,1))']');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'intervalEndTimes=[...\r\n');
fprintf(fid,'%9.4f %% %6d\r\n',[intervalEndTimes,(1:size(intervalEndTimes,1))']');
fprintf(fid,'];\r\n');
fprintf(fid,'\r\n');

fclose(fid);
end