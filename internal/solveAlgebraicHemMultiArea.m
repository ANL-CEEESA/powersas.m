function [stateNew,finalAlpha,alphaList,diff]=...
    solveAlgebraicHemMultiArea(SimData,SysData,SysPara,SysDataList,SysParaList,links,grossMaps,sysMaps,x0,xNew,SysDataNew)
% Compute HE coefficients for solving algebraic equations in multiple areas (parallel computation enabled)
%
% FUNCTION hemMachinePFSalientcontinueAlgebraicMultiArea
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
%	SysDataList - The list of sub system data
%	SysParaList - The list of sub system parameters
%	links - the links connecting the sub systems
%	grossMaps - Mapping from agglomerated system to subsystem
%	sysMaps - Mapping from subsystem to agglomerated system
%   x0 - Initial system state
%   xNew - New system state (algebraic variables need to be updated)
%	SysDataNew - New system data
%	
% OUTPUT
%	stateNew - Solved state
%   finalAlpha - The ending length of this segment of simulation
%   alphaList - Record of alphas
%   diff - A list of errors
%
%
[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
if nargin<11
    SysDataNew=[];
end
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

nSys=length(SysDataList);

for iSys=1:nSys
    [busx,swx,pvx,pqx,shuntx,linex,indx,zipx,synx,excx,tgx,agcx,cacx,clusterx]=unfoldSysData(SysDataList{iSys});
    [pqIncrx,pvIncrx,Rind0x,Rind1x,~,~,Rzip0x,Rzip1x,Ytr0x,Ytr1x,Ysh0x,Ysh1x,VspSq2x,MatGV0x,MatGV1x,MatGRhs0x,MatGRhs1x]=unfoldSysPara(SysParaList{iSys});
    
    nbusx=size(busx,1);
    nlinex=size(linex,1);
    nIndx=size(indx,1);
    nZipx=size(zipx,1);
    nSynx=size(synx,1);

    busTypex=zeros(nbusx,1);
    if isempty(pvx)
        pvx=zeros(0,6);
    end
    if isempty(pqx)
        pqx=zeros(0,6);
    end
    if isempty(shuntx)
        shuntx=zeros(0,7);
    end
    if isempty(swx)
        swx=zeros(0,13);
    end
    busTypex(pvx(:,1))=1;
    busTypex(swx(:,1))=2;

    iswx=find(busTypex==2);
    ipvx=find(busTypex~=0);
    ipqx=find(busTypex==0);
    npqx=size(ipqx,1);
    npvx=size(ipvx,1);

    yShuntx=zeros(nbusx,1);
    yShuntx(shuntx(:,1))=shuntx(:,5)+1j*shuntx(:,6);
    Ysh0x=Ysh0x+yShuntx;

%     VspSq2x=zeros(nbusx,1);
    % pqIncrS=zeros(size(pq,1),2);pvIncrS=zeros(size(pv,1),1);
    % Rzip0S=ones(nZip,1);Rzip1S=zeros(nZip,1);
    % Rind0S=ones(nInd,1);Rind1S=zeros(nInd,1);
    Reind0x=ones(nIndx,1);Reind1x=zeros(nIndx,1);
    
    SysParaList{iSys}=foldSysPara(pqIncrx,pvIncrx,Rind0x,Rind1x,Reind0x,Reind1x,Rzip0x,Rzip1x,Ytr0x,Ytr1x,Ysh0x,Ysh1x,VspSq2x,MatGV0x,MatGV1x,MatGRhs0x,MatGRhs1x);
end

% try
    [~,~,~,VSol,QSol,finalAlpha,alphaList,diff]=hemMachinePFmultiStageAlgebraicMultiArea(SimData,SysData,SysPara,SysDataList,SysParaList,links,grossMaps,sysMaps,x0,xNew,SysDataNew); % Note parameters
    stateNew(idxs.vIdx)=VSol;
    stateNew(idxs.qIdx)=QSol;
% catch ME
%     stateNew=[];
%     finalAlpha=0;
%     alphaList=[];
%     diff=0;
% end
addLog(['[INFO]Instant switching calc, final Alpha=',num2str(finalAlpha)],'INFO');
end