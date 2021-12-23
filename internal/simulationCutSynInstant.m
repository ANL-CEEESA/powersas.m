function [stateUpd,finalAlpha,alphaList,diff,SysDataUpd,SysParaUpd,maps]=...
    simulationCutSynInstant(busTag,SysData,SysDataBase,SysPara,SimData,maps,...
    x0,synCutIdx)
% Simulate cut synchronous generator
%
% FUNCTION simulationCutSynInstant
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
%   busTag - Tags for determining the partial system that will be kept 
%   SysData - System data for simulation
%   SysDataBase - Base system data for simulation
%   SysPara - Parameters representing the events happening in the system
%   SimData - Simulation parameters
%	*Map - The mapping of current data to base data
%	inv*Map - Inverse mapping, base data to current data
%   x0 - Initial system state
%   synCutIdx - Index of synchronous generator
%
% OUTPUT
%   stateNew - Solved new state
%   finalAlpha - The ending alpha
%   alphaList - A list of alphas
%   diff - Equation mismatches
%   SysDataUpd - Updated system data 
%	*Map - The mapping of current data to base data (updated)
%	inv*Map - Inverse mapping, base data to current data (updated)
%
[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
[busBase,swBase,pvBase,pqBase,shuntBase,lineBase,indBase,zipBase,synBase,excBase,tgBase,agcBase,cacBase,clusterBase]=unfoldSysData(SysDataBase);
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,fault]=unfoldSysPara(SysPara);

busMap=maps.busMap;invBusMap=maps.invBusMap;
indMap=maps.indMap;invIndMap=maps.invIndMap;
synMap=maps.synMap;invSynMap=maps.invSynMap;
lineMap=maps.lineMap;invLineMap=maps.invLineMap;

[~,~,~,~,~,~,~,~,method]=unfoldSimData(SimData);
AEmethod=round(10*mod(method,1));

synToTg=zeros(size(synBase,1),1);
synToTg(tgBase(:,1))=1:size(tgBase,1);
synToExc=zeros(size(synBase,1),1);
synToExc(excBase(:,1))=1:size(excBase,1);

pqIncrS=zeros(size(pq,1),2);pvIncrS=zeros(size(pv,1),1);
nbus=size(bus,1);
nZip=size(zip,1);
nInd=size(ind,1);
nSyn=size(syn,1);
Rzip0S=ones(nZip,1);Rzip1S=zeros(nZip,1);
Rind0S=ones(nInd,1);Rind1S=zeros(nInd,1);

[~,Ytr0,Ysh0,~,~,~,~]=getYMatrix(nbus,line,fault);

[V,Q,s,d,w,eq1,eq2,ed1,ed2,psid,psiq,Pm,Ef,Vavrm,Vavrr,Vavrf,Vavrref,tgovg,tgovm,tgovmech,f,dpg,qplt,vg]=...
    unfoldX(x0,SysData);

[MatGV0,~,MatGRhs0,~]=getLinearInterpolatorSyn(...
    syn,[],d,d,0*Ef,ed1,ed2,0*Ef,ed1,ed2,Ef,eq1,eq2,Ef,eq1,eq2,psid,psiq,psid,psiq);

Vgtemp=V(syn(:,1));
Ig=MatGRhs0(:,1)+1j*MatGRhs0(:,2)-...
    ((MatGV0(:,1).*real(Vgtemp)+MatGV0(:,2).*imag(Vgtemp))+1j*(MatGV0(:,3).*real(Vgtemp)+MatGV0(:,4).*imag(Vgtemp)));

Yg=-Ig./Vgtemp;

newSynIdx=synMap(synCutIdx);
activeNewSynCutIdx=newSynIdx(newSynIdx~=0);
synTag=zeros(size(syn,1),1);
synTag(activeNewSynCutIdx)=1;

activeSynRemainIdx=find(synTag==0);
invActiveSynRemainIdx=zeros(nSyn,1);
invActiveSynRemainIdx(activeSynRemainIdx)=1:length(activeSynRemainIdx);
YgCut=accumarray(syn(activeNewSynCutIdx,1),Yg(activeNewSynCutIdx),[nbus,1]);

syn=syn(activeSynRemainIdx,:);
d=d(activeSynRemainIdx);
w=w(activeSynRemainIdx);
eq1=eq1(activeSynRemainIdx);
eq2=eq2(activeSynRemainIdx);
ed1=ed1(activeSynRemainIdx);
ed2=ed2(activeSynRemainIdx);
psid=psid(activeSynRemainIdx);
psiq=psiq(activeSynRemainIdx);
Pm=Pm(activeSynRemainIdx);
Ef=Ef(activeSynRemainIdx);

excTag=synTag(exc(:,1));
tgTag=synTag(tg(:,1));

Vavrm=Vavrm(excTag==0);
Vavrr=Vavrr(excTag==0);
Vavrf=Vavrf(excTag==0);
Vavrref=Vavrref(excTag==0);
exc=exc(excTag==0,:);

tgovg=tgovg(tgTag==0);
tgovm=tgovm(tgTag==0);
tgovmech=tgovmech(tgTag==0);
tg=tg(tgTag==0,:);

synMap(invSynMap)=invActiveSynRemainIdx;
invSynMap=invSynMap(activeSynRemainIdx);

SysDataUpd=foldSysData(bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster);

x0=foldX(SysDataUpd,V,Q,s,d,w,eq1,eq2,ed1,ed2,psid,psiq,Pm,Ef,Vavrm,Vavrr,Vavrf,Vavrref,tgovg,tgovm,tgovmech,f,dpg,qplt,vg);

Ysh0=Ysh0+YgCut;
Ysh1=-YgCut;
SysParaSim=foldSysPara(pqIncrS,pvIncrS,Rind0S,Rind1S,[],[],Rzip0S,Rzip1S,Ytr0,0*Ytr0,Ysh0,Ysh1,[],[],[],[],[]);

if ~isempty(busTag)&&~isempty(find(busTag==0,1))
    % Purge system with busTag
    [SysDataUpdFtr,busFtr,swFtr,pvFtr,pqFtr,shuntFtr,lineFtr,indFtr,zipFtr,synFtr,excFtr,tgFtr,...
        busMap,indMap,synMap,lineMap,invBusMap,invIndMap,invSynMap,invLineMap]=...
        filterSysData(busTag,SysDataUpd,busMap,indMap,synMap,lineMap,invBusMap,invIndMap,invSynMap,invLineMap);
    SysParaSim=filterSysPara(SysParaSim,busFtr,swFtr,pvFtr,pqFtr,shuntFtr,lineFtr,indFtr,zipFtr,synFtr,excFtr,tgFtr);
    SysPara=filterSysPara(SysPara,busFtr,swFtr,pvFtr,pqFtr,shuntFtr,lineFtr,indFtr,zipFtr,synFtr,excFtr,tgFtr);
    x0=filterX(SysDataUpd,SysDataUpdFtr,x0,busFtr,swFtr,pvFtr,pqFtr,shuntFtr,lineFtr,indFtr,zipFtr,synFtr,excFtr,tgFtr);
    SysDataUpd=SysDataUpdFtr;
end

if AEmethod==0
    [stateNew,finalAlpha,alphaList,diff]=solveAlgebraicHem(SimData,SysDataUpd,SysParaSim,x0,x0);
else
    [stateNew,flag,diff,loop]=solveAlgebraicNR(SimData,SysDataUpd,SysParaSim,x0,x0);
    if flag==0
        finalAlpha=1;
    else
        finalAlpha=0;
    end
    alphaList=linspace(0,finalAlpha,loop);
end

stateUpd=stateNew;
SysParaUpd=SysPara;

maps.busMap=busMap;maps.invBusMap=invBusMap;
maps.indMap=indMap;maps.invIndMap=invIndMap;
maps.synMap=synMap;maps.invSynMap=invSynMap;
maps.lineMap=lineMap;maps.invLineMap=invLineMap;

end