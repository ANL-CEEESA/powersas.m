function [stateUpd,t,finalAlpha,alphaList,diffList,SysDataUpd,exitFlag,SimDataUpd]=generalDynSimulation(SimData,SysData,SysPara,x0,zipRamp,busMap)
% General interface for invoking dynamic simulation for a segment of time
%
% FUNCTION generalDynSimulation
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
%   zipRamp - Special data field for defining ramping data of ZIP load
%   busMap - Mapping of current system bus numbers to the base system
%   (useful for the cases where only partial system is studied)
%
% OUTPUT
%   stateUpd - A list of states in the order of time
%   t - A list of time points (starting with 0)
%   finalAlpha - The ending length of this segment of simulation
%   alphaList - Equivalent to t
%   diffList - A list of errors
%   SysDataUpd - Updated system data 
%   exitFlag - 
%       0  - Success and normally exit
%       -1 - Fail to finish (due to computation errors or failures)
%       1  - System in absolute steady state
%       2  - Generator transients faded away
%       3  - Suggest using error reduction mode
%   SimDataUpd - Updated simulation data
%
%%

if nargin<5
    zipRamp=[];
end
[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);

if nargin<6
    busMap=bus(:,1);
end

[pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,...
    MatGV0,MatGV1,MatGRhs0,MatGRhs1,Tmech1,Varref1,Ef1,Pm1,Eq11,fault]=unfoldSysPara(SysPara);

[~,~,~,~,~,~,~,~,method]=unfoldSimData(SimData);
DEmethod=floor(method);

nbus=size(bus,1);
nline=size(line,1);
nInd=size(ind,1);
nZip=size(zip,1);
nSyn=size(syn,1);
nTg=size(tg,1);
nExc=size(exc,1);
nAgc=size(agc,1);
nCac=size(cac,1);
nCluster=size(cluster,1);

% if nargin<7
%     mode=0;
% end

if isempty(Rzip0);Rzip0=ones(nZip,1);end
if isempty(Rzip1);Rzip1=zeros(nZip,1);end
if ~isempty(zipRamp)
    zipRamp(:,1)=busMap(zipRamp(:,1));
    zip=[zip;zipRamp];
    Rzip0=[Rzip0;zeros(size(zipRamp,1),1)];Rzip1=[Rzip1;ones(size(zipRamp,1),1)];
end

if isempty(pqIncr);pqIncr=zeros(size(pq,1),2);end
if isempty(pvIncr);pvIncr=zeros(size(pv,1),1);end
if isempty(Rind0);Rind0=ones(nInd,1);end
if isempty(Rind1);Rind1=zeros(nInd,1);end
if isempty(Reind0);Reind0=ones(nInd,1);end
if isempty(Reind1);Reind1=zeros(nInd,1);end
if isempty(Ytr0)||isempty(Ytr1)||isempty(Ysh0)||isempty(Ysh1)
    [~,Ytr0x,Ysh0x,~,~,~,~]=getYMatrix(nbus,line,fault);
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

SysPara=foldSysPara(pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,MatGV0,MatGV1,MatGRhs0,MatGRhs1,Tmech1,Varref1,Ef1,Pm1,Eq11);

SysDataTmp=foldSysData(bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster);
if DEmethod==0
    [stateUpd,t,finalAlpha,alphaList,diffList,exitFlag]=simulationTimeDomainHem(SimData,SysDataTmp,SysPara,x0);
    alphaList=cumsum(alphaList);
    SimDataUpd=SimData;
elseif DEmethod==4
    [stateUpd,t,finalAlpha,alphaList,diffList,exitFlag]=simulationTimeDomainDT(SimData,SysDataTmp,SysPara,x0);
    alphaList=cumsum(alphaList);
    SimDataUpd=SimData;
else
    [stateUpd,t,diffList,nxtDt,exitFlag]=simulationTimeDomainNI(SimData,SysDataTmp,SysPara,x0);
    finalAlpha=t(end);
    alphaList=t;
    SimDataUpd=SimData;
    SimDataUpd.dt=nxtDt;
end

if ~isempty(zipRamp)
    zipExt=zeros(max(bus(:,1)),12);
    hasLoad=zeros(size(zipExt,1),1);
    zipExt(bus(:,1),1)=bus(:,1);zipExt(bus(:,1),2)=100;zipExt(bus(:,1),3)=bus(:,2);zipExt(bus(:,1),4)=60;
    zipExt(:,5)=accumarray(zip(:,1),zip(:,5).*(Rzip0+finalAlpha*Rzip1),[size(zipExt,1),1]);hasLoad=hasLoad|(zipExt(:,5)~=0);
    zipExt(:,6)=accumarray(zip(:,1),zip(:,6).*(Rzip0+finalAlpha*Rzip1),[size(zipExt,1),1]);hasLoad=hasLoad|(zipExt(:,6)~=0);
    zipExt(:,7)=accumarray(zip(:,1),zip(:,7).*(Rzip0+finalAlpha*Rzip1),[size(zipExt,1),1]);hasLoad=hasLoad|(zipExt(:,7)~=0);
    zipExt(:,8)=accumarray(zip(:,1),zip(:,8).*(Rzip0+finalAlpha*Rzip1),[size(zipExt,1),1]);hasLoad=hasLoad|(zipExt(:,8)~=0);
    zipExt(:,9)=accumarray(zip(:,1),zip(:,9).*(Rzip0+finalAlpha*Rzip1),[size(zipExt,1),1]);hasLoad=hasLoad|(zipExt(:,9)~=0);
    zipExt(:,10)=accumarray(zip(:,1),zip(:,10).*(Rzip0+finalAlpha*Rzip1),[size(zipExt,1),1]);hasLoad=hasLoad|(zipExt(:,10)~=0);
    zipExt(:,12)=1;
    zip=zipExt(hasLoad~=0,:);
end
if ~isempty(pq)
    pq(:,[4,5])=pq(:,[4,5])+pqIncr(:,1:2)*finalAlpha;
    if size(pqIncr,2)>=4
        pq(:,[4,5])=pq(:,[4,5])+pqIncr(:,3:4)*finalAlpha*finalAlpha;
    end
end
if ~isempty(pv)
    pv(:,4)=pv(:,4)+pvIncr*finalAlpha;
end
if ~isempty(ind)
    ind(:,15:17)=ind(:,15:17).*repmat(Rind0+finalAlpha*Rind1,1,3);
end
if ~isempty(shunt)
    shunt(:,[5,6])=shunt(:,[5,6])+finalAlpha*[real(Ysh1(shunt(:,1))),imag(Ysh1(shunt(:,1)))];
end

SysDataUpd=foldSysData(bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster);

end