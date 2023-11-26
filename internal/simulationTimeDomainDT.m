function [stateCurve,t,finalAlpha,alphaList,diffList,exitFlag]=simulationTimeDomainHem(SimData,SysData,SysPara,x0)
% Numerical integration approach for dynamic simulation
%
%%FUNCTION simulationTimeDomainNI
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
% INPUT
%   SimData - Simulation parameters
%   SysData - System data for simulation
%   SysPara - Parameters representing the events happening in the system
%   x0 - Initial system state
%
% OUTPUT
%   stateCurve - A list of states in the order of time
%   t - A list of time points (starting with 0)
%	finalAlpha - Final value of alpha
%	alphaList - List of alpha (from multi-stage HE)
%   exitFlag - 
%       0  - Success and normally exit
%       -1 - Fail to finish (due to computation errors or failures)
%       1  - System in absolute steady state
%       2  - Generator transients faded away
%       3  - Suggest using error reduction mode
%
% MODES (controlled by SimData.method)
%   SimData.method has format X.YZ
%   X - Method for solving differential equations (DE)
%       0 - HE
%       1 - Modified Euler
%       2 - RK4
%       3 - Trapezoidal
%
%   Y - Method for solving algebraic equations (AE)
%       0 - HE
%       1 - NR
%
%   Z - Variable step (only applies to X>1)
%       0 - Fixed step
%       1 - Adaptive step
% 
[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
[nState,idxs]...
    =getIndexDyn(SysData);
% [V0,Q0,s0,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,tgovmech0]=unfoldX(x0,SysData);
% [~,~,sNew,dNew,~,eq1New,eq2New,ed1New,ed2New,psidNew,psiqNew,~,EfNew,~,~,~,~,~,~,~]=unfoldX(xNew,SysData);
[maxTime,segTime,dt,nlvl,taylorN,alphaTol,diffTol,diffTolMax,method,varOpt]=unfoldSimData(SimData);
[pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,~,~,~,~,Tmech1,Varref1,Ef1,Pm1,Eq11]=unfoldSysPara(SysPara);

nbus=size(bus,1);
nline=size(line,1);
nInd=size(ind,1);
nZip=size(zip,1);
nSyn=size(syn,1);
nTg=size(tg,1);
nExc=size(exc,1);

if isempty(pqIncr);pqIncr=zeros(size(pq,1),2);end
if isempty(pvIncr);pvIncr=zeros(size(pv,1),1);end
if isempty(Rind0);Rind0=ones(nInd,1);end
if isempty(Rind1);Rind1=zeros(nInd,1);end
if isempty(Reind0);Reind0=ones(nInd,1);end
if isempty(Reind1);Reind1=zeros(nInd,1);end
if isempty(Rzip0);Rzip0=ones(nZip,1);end
if isempty(Rzip1);Rzip1=zeros(nZip,1);end
if isempty(Ytr0)||isempty(Ytr1)||isempty(Ysh0)||isempty(Ysh1)
    [~,Ytr0x,Ysh0x,~,~,~,~]=getYMatrix(nbus,line);
    if isempty(Ytr0);Ytr0=Ytr0x;end
    if isempty(Ytr1);Ytr1=0*Ytr0;end
    if isempty(Ysh0);Ysh0=Ysh0x;end
    if isempty(Ysh1);Ysh1=0*Ysh0;end
end
if isempty(Tmech1);Tmech1=zeros(nTg,1);end
if isempty(Varref1);Varref1=zeros(nExc,1);end
if isempty(Ef1);Ef1=zeros(nSyn,1);end
if isempty(Pm1);Pm1=zeros(nSyn,1);end
if isempty(Eq11);Eq11=zeros(nSyn,1);end

busType=zeros(nbus,1);
if isempty(pv)
    pv=zeros(0,6);
end
if isempty(pq)
    pq=zeros(0,6);
end
if isempty(shunt)
    shunt=zeros(0,7);
end
if isempty(sw)
    sw=zeros(0,13);
end
busType(pv(:,1))=1;
busType(sw(:,1))=2;

isw=find(busType==2);
ipv=find(busType~=0);
ipq=find(busType==0);
npq=size(ipq,1);
npv=size(ipv,1);

% yShunt=zeros(nbus,1);
% yShunt(shunt(:,1))=shunt(:,5)+1j*shunt(:,6);
% Ysh0=Ysh0+yShunt;
if isempty(VspSq2)
    V0=x0(idxs.vIdx);
    VspSq2=[abs(V0).*abs(V0),zeros(nbus,1)];
end

segAlpha=min([maxTime,segTime]);

SysParax=foldSysPara(pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,[],[],[],[],Tmech1,Varref1,Ef1,Pm1,Eq11);
SimDatax=foldSimData(maxTime,segAlpha,dt,nlvl,taylorN,alphaTol,diffTol,diffTolMax,method,varOpt);

[t,stateCurve,finalAlpha,alphaList,diffList,exitFlag]=dtmMachinePFmultiStageDyn(SimDatax,SysData,SysParax,x0);

end