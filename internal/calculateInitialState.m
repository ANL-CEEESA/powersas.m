function [SysDataUpd,x0,finalAlpha,alphaList,diff,t,stateCurve]=calculateInitialState(SysData,SimData,SysPara)
% Calculate the initial steady-state solution of a system
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
%
% OUTPUT
%	SysDataUpd - Updated system states
%   x0 - Solved system states
%	finalAlpha - Final value of alpha
%	alphaList - List of alpha (from multi-stage HE)
%	diff - mismatch vector of the equations
%
[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
[maxTime,segTime,dt,nlvl,taylorN,alphaTol,diffTol,diffTolMax,method]=unfoldSimData(SimData);
[pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,...
    MatGV0,MatGV1,MatGRhs0,MatGRhs1,Tmech1,Varref1,Ef1,Pm1,Eq11]=unfoldSysPara(SysPara);
ltc=[];sup=[];dem=[];busName=[];
% pv=[]; % with synchronous generators, the PV buses are not useful
[nIslands,islands]=searchIslands(bus(:,1),line(:,[1,2]));
if nIslands>1
    addLog('There are more than 1 islands in the system. Only the largest island will be analyzed.','INFO');
    maxIsland=findLargestIsland(nIslands,islands,bus,[],pq,sw);
    [bus,~,pq,sw,line,ltc,sup,dem,shunt,ind,zip,syn,busName,newToOld,oldToNew,pvInd,pqInd,swInd,lineInd,ltcInd,supInd,demInd,shuntInd,indInd,zipInd,synInd] ...
        =forgeLargestIsland(maxIsland,nIslands,islands,bus,[],pq,sw,line,ltc,sup,dem,shunt,ind,zip,syn,busName);
end

nbus=size(bus,1);
nline=size(line,1);
nIslands=1;
islands=ones(nbus,1);
if ~isempty(sw)
    refs=sw(1,1);
else
    refs=zeros(0,1);
end

synTag=zeros(nbus,1);
if ~isempty(syn);synTag(syn(:,1))=1;end
sw=sw(synTag(sw(:,1))==0,:);

% [bus,~,pq,sw,line,ltc,sup,dem,shunt,busName,ind,zip,syn,pm]=regulatePFData(bus,pv,pq,sw,line,ltc,sup,dem,shunt,busName,ind,zip,syn);
pv=pv(synTag(pv(:,1))==0,:);

if isempty(syn)&&isempty(sw)
    SysDataUpd=foldSysData(bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster);
    x0=[];
    finalAlpha=0;
    alphaList=[];
    diff=[];
    t=[];
    stateCurve=[];
    addLog('Steady-state analysis needs at Swing bus or synchronous machine(s) in the system.','ERROR');
    return;
end

[bus,pv,pq,sw,line,shunt,ind,zip,syn,exc,tg,agc,cac,cluster,pm,oldToNew,newToOld]=regulateSystemData(bus,pv,pq,sw,line,shunt,ind,zip,syn,exc,tg,agc,cac,cluster);
if ~isempty(Pm1)
    pm=Pm1;
end
nbus=size(bus,1);
nline=size(line,1);
% nlvl=10;
pq0=pq;
pv0=pv;
p0=pq(:,4);
q0=pq(:,5);
shunt0=shunt;
ind0=ind;
zip0=zip;
if ~isempty(pv);pvIncr=pv(:,4);else pvIncr=zeros(0,1); end
if ~isempty(pq);pqIncr=pq(:,[4,5]);else pqIncr=zeros(0,2); end
Rind0=zeros(size(ind,1),1);
Rind1=ones(size(ind,1),1);
Reind0=zeros(size(ind,1),1);
Reind1=ones(size(ind,1),1);
Rzip0=zeros(size(zip,1),1);
Rzip1=ones(size(zip,1),1);

busType=zeros(nbus,1);
if ~isempty(pv);busType(pv(:,1))=1;end
if ~isempty(sw);busType(sw(:,1))=2;end
if ~isempty(syn)
    nm=accumarray(syn(:,1),ones(size(syn,1),1),[nbus,1]);
    swGenBus=find(busType==2&nm~=0);
    if ~isempty(swGenBus)
        idxSwGen=find(syn(:,1)==swGenBus,1);
    else
        [~,idxSwGen]=min(pm);
    end
    pShare=zeros(size(syn,1),1);
    pShare(idxSwGen)=1/length(idxSwGen);
else
    pShare=[];
end
pShare=1/size(syn,1)*ones(size(syn,1),1);

% weightP(end)=0;
% pShare=(weightP)/sum(weightP);

V0=ones(nbus,1);
s0=zeros(size(ind,1),1);
Q0=zeros(nbus,1);
d0=zeros(size(syn,1),1);
Ef0=ones(size(syn,1),1);
if isempty(Ef1);Ef1=zeros(size(syn,1),1);end

[Y,Ytr,Ysh,ytrfr,ytrto,yshfr,yshto]=getYMatrix(nbus,line);
yShunt=zeros(nbus,1);
if ~isempty(shunt)
    yShunt(shunt(:,1))=shunt(:,5)+1j*shunt(:,6);
end
% if ~isempty(zip)%zipMode=0
%     yShunt=yShunt+accumarray(zip(:,1),(zip(:,5)+1j*zip(:,8)).*zip(:,12),[nbus,1]);
% end
Ysh=Ysh+yShunt;
Ysh0=Ysh*0;
Ysh1=Ysh;

PmX=zeros(size(syn,1),2);
PmX(:,2)=pm;

VspSq2=[ones(nbus,1),zeros(nbus,1)];
if ~isempty(pv);VspSq2(pv(:,1),2)=pv(:,5).*pv(:,5)-1;end
if ~isempty(sw);VspSq2(sw(:,1),2)=sw(:,4).*sw(:,4)-1;end
% taylorN=4;

pq(:,[4,5])=0;
if ~isempty(pv);pv(:,4)=0;end
% alphaTol=0.001;
% diffTol=1e-6;
% diffTolMax=1e-2;
startTag1=tic;

SysDatax=foldSysData(bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster);
[nState,idxs]...
    =getIndexDyn(SysDatax);
x=zeros(nState,1);
x(idxs.vIdx)=V0;
x(idxs.sIdx)=s0;
x(idxs.qIdx)=Q0;
x(idxs.deltaIdx)=d0;
x(idxs.efIdx)=Ef0;
x(idxs.pgIdx)=PmX(:,1);
SysParax=foldSysPara(pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr,[],Ysh0,Ysh1,VspSq2,[],[],[],[],[],[],Ef1,PmX(:,2),[],[]);
SysParax.nIslands=nIslands;
SysParax.islands=islands;
SysParax.refs=refs; 

addLog('###### STEADY-STATE ######','INFO');
[VSol,QSol,sSol,dSol,finalAlpha,alphaList,diff,t,stateCurve]=hemMachinePFmultiStage(SimData,SysDatax,SysParax,x,pShare);
addLog(['# STEADY-STATE # Comp. Time = ',num2str(toc(startTag1)),'(s).'],'INFO');

Ef=Ef0+1*Ef1;
pq=pq0;
pv=pv0;
zip=zip0;
ind=ind0;
shunt=shunt0;

SysDataUpd=foldSysData(bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster);

x=zeros(nState,1);
x(idxs.vIdx)=VSol;
x(idxs.qIdx)=QSol;
x(idxs.sIdx)=sSol;
x(idxs.deltaIdx)=dSol;
x(idxs.omegaIdx)=0.0;
x(idxs.efIdx)=Ef;
dx=0.0*x;
x0=perpareInitialState(SysDataUpd,x,dx);
end