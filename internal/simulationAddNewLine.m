function [stateUpd,finalAlpha,alphaList,diff,SysDataUpd,maps]=...
    simulationAddNewLine(SysData,SysDataBase,SysPara,SimData,maps,x0,lineAddIdx)
% Simulate add line (connects new islands or buses)
%
% FUNCTION simulationAddNewLine
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
%   SysData - System data for simulation
%   SysDataBase - Base system data for simulation
%   SysPara - Parameters representing the events happening in the system
%   SimData - Simulation parameters
%	*Map - The mapping of current data to base data
%	inv*Map - Inverse mapping, base data to current data
%   x0 - Initial system state
%   lineAddIdx - Index of added line
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
[nState,idxs]...
    =getIndexDyn(SysData); 
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,fault]=unfoldSysPara(SysPara);

busMap=maps.busMap;invBusMap=maps.invBusMap;
indMap=maps.indMap;invIndMap=maps.invIndMap;
synMap=maps.synMap;invSynMap=maps.invSynMap;
lineMap=maps.lineMap;invLineMap=maps.invLineMap;

[~,~,~,~,~,~,~,~,method]=unfoldSimData(SimData);
AEmethod=round(10*mod(method,1));

pqIncrS=zeros(size(pq,1),2);pvIncrS=zeros(size(pv,1),1);
nbus=size(bus,1);
nZip=size(zip,1);
nInd=size(ind,1);
Rzip0S=ones(nZip,1);Rzip1S=zeros(nZip,1);
Rind0S=ones(nInd,1);Rind1S=zeros(nInd,1);

lineAdd=lineBase(lineAddIdx,:);
[connectBus,pos]=max([busMap(lineAdd(:,1)),busMap(lineAdd(:,2))],[],2);
[~,Ytr0,Ysh0,~,~,~,~]=getYMatrix(nbus,line,fault);

r=lineAdd(:,8);
x=lineAdd(:,9);
z=r+1j*x;
b=lineAdd(:,10);
chrg1 = 1j*0.5*lineAdd(:,16).*b;
chrg2 = 1j*0.5*lineAdd(:,16).*b;
y = lineAdd(:,16)./z;
lineAdd(lineAdd(:,11)==0,11)=1;
ts = lineAdd(:,11).*exp(1j*lineAdd(:,12)*pi/180);
ts2= ts.*conj(ts);
ytrfr=y./conj(ts);
ytrto=y./ts;
yshfr=(y+chrg1)./ts2-ytrfr;
yshto=y+chrg2-ytrto;

nbusAdd=length(connectBus);
yShAddSeg=zeros(nbusAdd,3);
yShAddSeg(pos==1,1)=yshfr(pos==1);yShAddSeg(pos==1,2)=ytrfr(pos==1);yShAddSeg(pos==1,3)=yshto(pos==1);
yShAddSeg(pos==2,1)=yshto(pos==2);yShAddSeg(pos==2,2)=ytrto(pos==2);yShAddSeg(pos==2,3)=yshfr(pos==2);

yShAdd=yShAddSeg(:,1)+1./(1./yShAddSeg(:,2)+1./yShAddSeg(:,3));
Ysh1=accumarray(connectBus,yShAdd,[nbus,1]);

SysPara=foldSysPara(pqIncrS,pvIncrS,Rind0S,Rind1S,[],[],Rzip0S,Rzip1S,Ytr0,0*Ytr0,Ysh0,Ysh1,[],[],[],[],[]);
xNew=x0;

if AEmethod==0
    [stateNew,finalAlpha,alphaList,diff]=solveAlgebraicHem(SimData,SysData,SysPara,x0,xNew);
else
    [stateNew,flag,diff,loop]=solveAlgebraicNR(SimData,SysData,SysPara,x0,xNew);
    if flag==0
        finalAlpha=1;
    else
        finalAlpha=0;
    end
    alphaList=linspace(0,finalAlpha,loop);
end

[V,Q,s,d,w,eq1,eq2,ed1,ed2,psid,psiq,Pm,Ef,Vavrm,Vavrr,Vavrf,Vavrref,tgovg,tgovm,tgovmech,f,dpg,qplt,vg]=unfoldX(stateNew,SysData);

newVolt=stateNew(idxs.vIdx(connectBus)).*yShAddSeg(:,2)./(yShAddSeg(:,2)+yShAddSeg(:,3));
newBusNum=zeros(nbusAdd,1);
newBusNum(pos==1)=lineAdd(pos==1,2);newBusNum(pos==2)=lineAdd(pos==2,1);
busMapAppend=zeros(size(busBase,1),1);
busMapAppend(newBusNum)=(nbus+1):(nbus+nbusAdd);
busMap(newBusNum)=(nbus+1):(nbus+nbusAdd);
invBusMap=[invBusMap;newBusNum];
busAdd=busBase(newBusNum,:);busAdd(:,1)=busMapAppend(busAdd(:,1));bus=[bus;busAdd];
swAdd=swBase(busMapAppend(swBase(:,1))~=0,:);swAdd(:,1)=busMapAppend(swAdd(:,1));sw=[sw;swAdd];
pvAdd=pvBase(busMapAppend(pvBase(:,1))~=0,:);pvAdd(:,1)=busMapAppend(pvAdd(:,1));pv=[pv;pvAdd];
% pq=[pq;pqBase(busMapAppend(pqBase(:,1))~=0,:)];
% shunt=[shunt;shuntBase(busMapAppend(shuntBase(:,1))~=0,:)];
nline=size(line,1);
lineAdd(:,1:2)=[busMap(lineAdd(:,1)),busMap(lineAdd(:,2))];
line=[line;lineAdd];
lineMap(lineAddIdx)=(nline+1):(nline+size(lineAdd,1));
invLineMap=[invLineMap;lineAddIdx];

busTag=zeros(max([invBusMap;agcBase(:,1)]),1);
busTag(invBusMap)=1;
agc=agcBase(busTag(agcBase(:,1))~=0,:);
agc(:,1)=busMap(agc(:,1));

SysDataUpd=foldSysData(bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster);

newV=[V;newVolt];
newQ=[Q;zeros(size(newVolt))];
newf=[f;f(connectBus)];
newDpg=[dpg;zeros(size(newVolt))];
stateUpd=foldX(SysDataUpd,newV,newQ,s,d,w,eq1,eq2,ed1,ed2,psid,psiq,Pm,Ef,Vavrm,Vavrr,Vavrf,Vavrref,tgovg,tgovm,tgovmech,newf,newDpg,qplt,vg);

maps.busMap=busMap;maps.invBusMap=invBusMap;
maps.indMap=indMap;maps.invIndMap=invIndMap;
maps.synMap=synMap;maps.invSynMap=invSynMap;
maps.lineMap=lineMap;maps.invLineMap=invLineMap;

end