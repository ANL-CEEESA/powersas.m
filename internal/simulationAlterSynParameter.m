function [stateUpd,finalAlpha,alphaList,diff,SysDataUpd,SysParaUpd,maps]=...
    simulationAlterSynParameter(busTag,SysData,SysDataBase,SysPara,SimData,maps,...
    x0,synNew)
% Simulate suddenly change synchornous generator parameter 
%
% FUNCTION simulationAlterSynParameter
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
%   synNew - Data of updated synchronous generator
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
[busBase,swBase,pvBase,pqBase,shuntBase,lineBase,indBase,zipBase,synBase,excBase,tgBase]=unfoldSysData(SysDataBase);
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,fault]=unfoldSysPara(SysPara);

busMap=maps.busMap;invBusMap=maps.invBusMap;
indMap=maps.indMap;invIndMap=maps.invIndMap;
synMap=maps.synMap;invSynMap=maps.invSynMap;
lineMap=maps.lineMap;invLineMap=maps.invLineMap;

[~,~,~,~,~,~,~,~,method]=unfoldSimData(SimData);
AEmethod=round(10*mod(method,1));

% synToTg=zeros(size(synBase,1),1);
% synToTg(tgBase(:,1))=1:size(tgBase,1);
% synToExc=zeros(size(synBase,1),1);
% synToExc(excBase(:,1))=1:size(excBase,1);

pqIncrS=zeros(size(pq,1),2);pvIncrS=zeros(size(pv,1),1);
nbus=size(bus,1);
nZip=size(zip,1);
nInd=size(ind,1);
Rzip0S=ones(nZip,1);Rzip1S=zeros(nZip,1);
Rind0S=ones(nInd,1);Rind1S=zeros(nInd,1);

[~,Ytr0,Ysh0,~,~,~,~]=getYMatrix(nbus,line,fault);

% synAdd=synBase(synAddIdx,:);
% synAdd(:,1)=busMap(synAdd(:,1));
% nSynAdd=size(synAdd,1);
% syn=[syn;synAdd];
% addSynMap=zeros(size(synMap));
% addSynMap(synAddIdx)=(length(invSynMap)+1):(length(invSynMap)+length(synAddIdx));
% synMap(synAddIdx)=(length(invSynMap)+1):(length(invSynMap)+length(synAddIdx));
% invSynMap=[invSynMap;synAddIdx];

[V,Q,s,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,tgovmech0,f0,qplt0,vg0]=...
    unfoldX(x0,SysData);

% synIdx=synAdd(:,1);
% wgb=synAdd(:,4);
% model=synAdd(:,5);
% Xgl=synAdd(:,6);
% Rga=synAdd(:,7);
% Xgd=synAdd(:,8);
% Xgd1=synAdd(:,9);
% Xgd2=synAdd(:,10);
% Tgd1=synAdd(:,11);
% Tgd2=synAdd(:,12);
% Xgq=synAdd(:,13);
% Xgq1=synAdd(:,14);
% Xgq2=synAdd(:,15);
% Tgq1=synAdd(:,16);
% Tgq2=synAdd(:,17);
% Mg=synAdd(:,18);
% Dg=synAdd(:,19);
% TgAA=synAdd(:,24);
% gammad=Tgd2./Tgd1.*Xgd2./Xgd1.*(Xgd-Xgd1);
% gammaq=Tgq2./Tgq1.*Xgq2./Xgq1.*(Xgq-Xgq1);
% 
% VbusAdd=V(synAdd(:,1));
% CbusAdd=real(VbusAdd);DbusAdd=imag(VbusAdd);
% cdAdd=cos(d0Add);sdAdd=sin(d0Add);
% VdBusAdd=sdAdd.*CbusAdd-cdAdd.*DbusAdd;VqBusAdd=cdAdd.*CbusAdd+sdAdd.*DbusAdd;

% d0Tmp=atan2(DbusAdd,CbusAdd);
% d0Add(isnan(d0Add))=d0Tmp(isnan(d0Add));
% VqBusTmp=sqrt(CbusAdd.*CbusAdd+DbusAdd.*DbusAdd);VdBusTmp=0*VqBusTmp;
% Ef0(isnan(Ef0))=VqBusTmp(isnan(Ef0));
% 
% ed0Add=VdBusTmp;ed10Add=VdBusTmp;ed20Add=VdBusTmp;
% edNewAdd=zeros(nSynAdd,1);ed1NewAdd=zeros(nSynAdd,1);ed2NewAdd=zeros(nSynAdd,1);
% eq0Add=VqBusTmp;eq10Add=VqBusTmp;eq20Add=VqBusTmp;
% eqNewAdd=Ef0;eq1NewAdd=Ef0.*(1-TgAA./Tgd1);eq2NewAdd=Ef0;
% eq1NewAdd(model==4|model==3|model==2)=Ef0(model==4|model==3|model==2);
% psiq0Add=VdBusTmp;psid0Add=VqBusTmp;
% psiqNewAdd=zeros(nSynAdd,1);psidNewAdd=Ef0;
% w0Add=zeros(nSynAdd,1);wNewAdd=zeros(nSynAdd,1);
%
% PmNew=[Pm;Pm0];
% EfNew=[Ef;Ef0];
% 
% d0=[d;d0Tmp];dNew=[d;d0Add];
% w0=[w;w0Add];wNew=[w;wNewAdd];
% ed0=[0*Ef;ed0Add];edNew=[0*Ef;edNewAdd];
% ed10=[ed1;ed10Add];ed1New=[ed1;ed1NewAdd];
% ed20=[ed2;ed20Add];ed2New=[ed2;ed2NewAdd];
% eq0=[Ef;eq0Add];eqNew=[Ef;eqNewAdd];
% eq10=[eq1;eq10Add];eq1New=[eq1;eq1NewAdd];
% eq20=[eq2;eq20Add];eq2New=[eq2;eq2NewAdd];
% psid0=[psid;psid0Add];psidNew=[psid;psidNewAdd];
% psiq0=[psiq;psiq0Add];psiqNew=[psiq;psiqNewAdd];

[MatGV0,MatGV1,MatGRhs0,MatGRhs1]=getLinearInterpolatorSyn(...
    syn,[],d0,d0,0*Ef0,ed10,ed20,0*Ef0,ed10,ed20,Ef0,eq10,eq20,Ef0,eq10,eq20,psid0,psiq0,psid0,psiq0);

SysParaSim=foldSysPara(pqIncrS,pvIncrS,Rind0S,Rind1S,[],[],Rzip0S,Rzip1S,Ytr0,0*Ytr0,Ysh0,0*Ysh0,[],MatGV0,MatGV1,MatGRhs0,MatGRhs1);
SysDataTmp=foldSysData(bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster);
SysDataTmpNew=foldSysData(bus,sw,pv,pq,shunt,line,ind,zip,synNew,exc,tg,agc,cac,cluster);
x1=foldX(SysDataTmp,V,Q,s,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,...
    Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,tgovmech0,f0,dpg0,qplt0,vg0);

if ~isempty(busTag)&&~isempty(find(busTag==0,1))
    % Purge system with busTag
    [SysDataTmpNewFtr]=...
        filterSysData(busTag,SysDataTmp,busMap,indMap,synMap,lineMap,invBusMap,invIndMap,invSynMap,invLineMap);
    [SysDataTmpFtr,busFtr,swFtr,pvFtr,pqFtr,shuntFtr,lineFtr,indFtr,zipFtr,synFtr,excFtr,tgFtr,...
        busMap,indMap,synMap,lineMap,invBusMap,invIndMap,invSynMap,invLineMap]=...
        filterSysData(busTag,SysDataTmp,busMap,indMap,synMap,lineMap,invBusMap,invIndMap,invSynMap,invLineMap);
    SysParaSim=filterSysPara(SysParaSim,busFtr,swFtr,pvFtr,pqFtr,shuntFtr,lineFtr,indFtr,zipFtr,synFtr,excFtr,tgFtr);
    SysPara=filterSysPara(SysPara,busFtr,swFtr,pvFtr,pqFtr,shuntFtr,lineFtr,indFtr,zipFtr,synFtr,excFtr,tgFtr);
    x1=filterX(SysDataTmp,SysDataTmpFtr,x1,busFtr,swFtr,pvFtr,pqFtr,shuntFtr,lineFtr,indFtr,zipFtr,synFtr,excFtr,tgFtr);
    SysDataTmp=SysDataTmpFtr;
    SysDataTmpNew=SysDataTmpNewFtr;
end

if AEmethod==0
    [stateNew,finalAlpha,alphaList,diff]=solveAlgebraicHem(SimData,SysDataTmp,SysParaSim,x1,x1,SysDataTmpNew);
else
    [stateNew,flag,diff,loop]=solveAlgebraicNR(SimData,SysDataTmp,SysParaSim,x1,x1,SysDataTmpNew);
    if flag==0
        finalAlpha=1;
    else
        finalAlpha=0;
    end
    alphaList=linspace(0,finalAlpha,loop);
end

[Vupd,Qupd]=unfoldX(stateNew,SysDataTmpNew);

SysDataUpd=SysDataTmpNew;
SysParaUpd=SysPara;

stateUpd=foldX(SysDataUpd,Vupd,Qupd,s,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,...
    Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,tgovmech0,f0,dpg0,qplt0,vg0);

maps.busMap=busMap;maps.invBusMap=invBusMap;
maps.indMap=indMap;maps.invIndMap=invIndMap;
maps.synMap=synMap;maps.invSynMap=invSynMap;
maps.lineMap=lineMap;maps.invLineMap=invLineMap;

end