function [stateNew,finalAlpha,alphaList,diff,heData]=...
    solveAlgebraicHemDistFactorTest(SimData,SysData,SysPara,x0,xNew)
% [TEMPORARY] HE for solving the algebraic equations using distribution factors
%
% FUNCTION restorationAlgebraicHemDistFactorTest
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
% INPUT
%   SimData - Simulation parameters
%   SysData - System data for simulation
%   SysPara - Parameters representing the events happening in the system
%   x0 - Initial system state
%   xNew - Updated system state (algebraic variables to be solved)
%
% OUTPUT
%   stateNew - Updated state
%   finalAlpha - The ending length of this segment of simulation
%   alphaList - A list of alpha
%   diff - Equation mismatch
%	heData - HE coeffcients recorded
%

%
[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
[nState,idxs]...
    =getIndexDyn(SysData);
% [V0,Q0,s0,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,tgovmech0]=unfoldX(x0,SysData);
% [~,~,sNew,dNew,~,eq1New,eq2New,ed1New,ed2New,psidNew,psiqNew,~,EfNew,~,~,~,~,~,~,~]=unfoldX(xNew,SysData);
[maxTime,segTime,dt,nlvl,taylorN,alphaTol,diffTol,diffTolMax,method]=unfoldSimData(SimData);
[pqIncr,pvIncr,Rind0,Rind1,~,~,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,MatGV0,MatGV1,MatGRhs0,MatGRhs1]=unfoldSysPara(SysPara);
%
nbus=size(bus,1);
nline=size(line,1);
nInd=size(ind,1);
nZip=size(zip,1);
nSyn=size(syn,1);

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

yShunt=zeros(nbus,1);
yShunt(shunt(:,1))=shunt(:,5)+1j*shunt(:,6);
Ysh0=Ysh0+yShunt;

% VspSq2=zeros(nbus,1);
% pqIncrS=zeros(size(pq,1),2);pvIncrS=zeros(size(pv,1),1);
% Rzip0S=ones(nZip,1);Rzip1S=zeros(nZip,1);
% Rind0S=ones(nInd,1);Rind1S=zeros(nInd,1);
Reind0=ones(nInd,1);Reind1=zeros(nInd,1);

addLog('[START]Instant switching calc.','INFO');

stateNew=xNew;
SysPara=foldSysPara(pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,MatGV0,MatGV1,MatGRhs0,MatGRhs1);
try
    [Vx,Wx,Qx,fx,tempState,finalAlpha,alphaList,diff]=hemMachinePFmultiStageAlgebraic(SimData,SysData,SysPara,x0,xNew);
    stateNew(idxs.vIdx)=tempState(idxs.vIdx);
    stateNew(idxs.qIdx)=tempState(idxs.qIdx);
    stateNew(idxs.fIdx)=tempState(idxs.fIdx);
    heData={Vx,Wx,Qx,fx};
catch ME
    stateNew=[];
    finalAlpha=0;
    alphaList=[];
    diff=0;
    heData={[],[],[]};
end
addLog(['[INFO]Instant switching calc, final Alpha=',num2str(finalAlpha)],'INFO');
end