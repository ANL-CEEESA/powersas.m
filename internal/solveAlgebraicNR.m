function [stateNew,flag,diffRec,loop]=solveAlgebraicNR(SimData,SysData,SysPara,x0,xNew,SysDataNew)
% General Newton-Raphson computation for solving the algebraic equations
%
% FUNCTION restorationAlgebraicNR (will be renamed in a future version)
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
%	xNew - New system state (algebraic variables need to be solved)
%	SysDataNew - NOT USED
%
% OUTPUT
%	stateNew - Solved state (only algebraic variables changed)
%	flag
%		0 - returned normally
%		-1- Fail
%	diffRec - Final maximum difference
%	loop - number of loops to finish NR
%
% TODO % Check input argument SysDataNew if it is necessary

setting_func_restorationAlgebraicNR;

addLog('[#START#]Solving algebraic equations - NR.','INFO');

[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
if nargin>=6
    [~,~,~,~,~,~,~,~,synNew,~,~]=unfoldSysData(SysData);
else
    synNew=[];
end
[nState,idxs]...
    =getIndexDyn(SysData);
% [V0,Q0,s0,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,tgovmech0]=unfoldX(x0,SysData);
% [~,~,sNew,dNew,~,eq1New,eq2New,ed1New,ed2New,psidNew,psiqNew,~,EfNew,~,~,~,~,~,~,~]=unfoldX(xNew,SysData);
[maxTime,segTime,dt,nlvl,taylorN,alphaTol,diffTol,diffTolMax,method]=unfoldSimData(SimData);
AEmethod=round(10*mod(method,1));
if isempty(AEmethod)
    AEmethod=1;
end

[pqIncr,pvIncr,Rind0,Rind1,~,~,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,MatGV0,MatGV1,MatGRhs0,MatGRhs1]=unfoldSysPara(SysPara);
% The Y matrices may be with faults, so don't re-generate them!

[V0,Q0,s0,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,tgovmech0,f0,dpg0,qplt0,vg0]=unfoldX(x0,SysData);
eq0=Ef0;ed0=zeros(size(Ef0));
[~,~,sNew,dNew,~,eq1New,eq2New,ed1New,ed2New,psidNew,psiqNew,~,EfNew,~,~,~,~,~,~,~]=unfoldX(xNew,SysData);
eqNew=EfNew;edNew=zeros(size(EfNew));

[MatGV0,MatGV1,MatGRhs0,MatGRhs1]=getLinearInterpolatorSyn(syn,synNew,d0,dNew,ed0,ed10,ed20,edNew,ed1New,ed2New,eq0,eq10,eq20,eqNew,eq1New,eq2New,psid0,psiq0,psidNew,psiqNew);
nbus=size(bus,1);
if ~isempty(ind)
    [YshInd0,Yshind1]=getLinearInterpolatorInd(nbus,ind,s0,sNew);
    Ysh0=Ysh0+YshInd0;
    Ysh1=Ysh1+Yshind1;
end
[nIslands,islands,refs]=searchIslands(bus(:,1),line(:,1:2));

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

if isempty(VspSq2)
    V0=x0(idxs.vIdx);
    VspSq2=[abs(V0).*abs(V0),zeros(nbus,1)];
end

% Determine the frequency model of each island
freqTypeTag=zeros(nIslands,1);%0:sw,1:syn,2:steady-state f
freqKeptTag=zeros(nbus,1);
frefs=refs;
fswTag=zeros(nbus,1);
fsynTag=zeros(nbus,1);
fswTag(isw)=1;
fswTagxD=fswTag;
fsynTag(syn(:,1))=1;
D0=imag(V0);
for isl=1:nIslands
    if isempty(find(fswTag(islands==isl)==1, 1))
        if isempty(find(fsynTag(islands==isl)==1, 1))
            freqTypeTag(isl)=2;
            busesInIsland=find(islands==isl);
            [~,imin]=min(abs(D0(busesInIsland)));
            frefs(isl)=busesInIsland(imin(1));
            fswTagxD(frefs(isl))=1;
            freqKeptTag(busesInIsland)=1;
        else
            freqTypeTag(isl)=1;
        end
    end
end
freqKeptTagxRef=freqKeptTag;
freqKeptTagxRef(frefs)=0;
nFreqKept=sum(freqKeptTag);

if ~isempty(agc)
    agcExt=zeros(nbus,size(agc,2));
    agcExt(agc(:,1),:)=agc;
    fdk=agcExt(:,2)+agcExt(:,3); %1/R+D
else
    fdk=zeros(nbus,1);
end

yShunt=zeros(nbus,1);
yShunt(shunt(:,1))=shunt(:,5)+1j*shunt(:,6);
Ysh0=Ysh0+yShunt;
Ysh0=Ysh0+accumarray(zip(:,1),(Rzip0).*(zip(:,5)+1j*zip(:,8)),[nbus,1]);
Ysh1=Ysh1+accumarray(zip(:,1),(Rzip1).*(zip(:,5)+1j*zip(:,8)),[nbus,1]);

% VspSq2=zeros(nbus,1);
% pqIncrS=zeros(size(pq,1),2);pvIncrS=zeros(size(pv,1),1);
% Rzip0S=ones(nZip,1);Rzip1S=zeros(nZip,1);
% Rind0S=ones(nInd,1);Rind1S=zeros(nInd,1);
Reind0=ones(nInd,1);Reind1=zeros(nInd,1);

if ~isempty(Ytr1)
    YtrNew=Ytr0+Ytr1;
else
    YtrNew=Ytr0;
end
YshNew=Ysh0+Ysh1;
YNew=YtrNew+sparse(1:nbus,1:nbus,YshNew,nbus,nbus);

MatGVNew=MatGV0+MatGV1;
MatGRhsNew=MatGRhs0+MatGRhs1;

SNew=-accumarray(pq(:,1),pq(:,4)+1j*pq(:,5),[nbus,1])+accumarray(pv(:,1),pv(:,4),[nbus,1]);
SNew=SNew-accumarray(pq(:,1),pqIncr(:,1)+1j*pqIncr(:,2),[nbus,1])+accumarray(pv(:,1),pvIncr,[nbus,1]);
SNew=SNew-accumarray(zip(:,1),(Rzip0+Rzip1).*(zip(:,7)+1j*zip(:,10)),[nbus,1]);

INew=-accumarray(zip(:,1),(Rzip0+Rzip1).*(zip(:,6)-1j*zip(:,9)),[nbus,1]);
JNew=real(INew);
KNew=imag(INew);

GNew=real(YNew);BNew=imag(YNew);
J1=GNew+sparse(syn(:,1),syn(:,1),MatGVNew(:,1),nbus,nbus);
J2=-BNew+sparse(syn(:,1),syn(:,1),MatGVNew(:,2),nbus,nbus);
J3=BNew+sparse(syn(:,1),syn(:,1),MatGVNew(:,3),nbus,nbus);
J4=GNew+sparse(syn(:,1),syn(:,1),MatGVNew(:,4),nbus,nbus);
P2freq=sparse(1:nbus,1:nbus,-freqKeptTag.*fdk,nbus,nbus);
Freq2freq=sparse([1:nbus,1:nbus],[1:nbus,frefs(islands)'],[ones(1,nbus),-ones(1,nbus)],nbus,nbus);

Vtemp=V0;
ftemp=f0;
QTemp=Q0;

busTag=ones(nbus,1);
busTag(sw(:,1))=0;
jacRowTag=[busTag;busTag;freqKeptTagxRef];
jacColTag=[busTag;fswTagxD==0;freqKeptTag];

idxNonSw=find(busType~=2);
idxNonSwD=find(busType~=2&fswTagxD==1);

if AEmethod==1 % NR and damped NR
    for loop=1:maxLoop
        Ctemp=real(Vtemp);
        Dtemp=imag(Vtemp);
        Atemp=abs(Vtemp);
        Vgtemp=Vtemp(syn(:,1));
        Ig=accumarray(syn(:,1),MatGRhsNew(:,1)+1j*MatGRhsNew(:,2),[nbus,1])-...
            accumarray(syn(:,1),(MatGVNew(:,1).*real(Vgtemp)+MatGVNew(:,2).*imag(Vgtemp))+1j*(MatGVNew(:,3).*real(Vgtemp)+MatGVNew(:,4).*imag(Vgtemp)),[nbus,1]);
        diff=conj(SNew)+INew.*abs(Vtemp)+Ig.*conj(Vtemp)-(YNew*Vtemp).*conj(Vtemp);
        diffP=real(diff)+freqKeptTag.*(-fdk.*ftemp+dpg0);
        diffVQ=imag(diff);
        diffF=ftemp-ftemp(frefs(islands));
        if ~isempty(pv)
            diffVQ(pv(:,1))=abs(Vtemp(pv(:,1))).*abs(Vtemp(pv(:,1)))-sum(VspSq2(pv(:,1)),2);
        end
        tDiff=[diffP(idxNonSw);diffVQ(idxNonSw);diffF(freqKeptTagxRef==1)];
        if max(abs(tDiff))<diffTol
            QTemp=imag(diff);
            break;
        end
        dCtemp=sparse(1:nbus,1:nbus,Ctemp,nbus,nbus);
        dDtemp=sparse(1:nbus,1:nbus,Dtemp,nbus,nbus);
        Jac1=[dCtemp*J1+dDtemp*J3,dCtemp*J2+dDtemp*J4;...
            dCtemp*J3-dDtemp*J1,dCtemp*J4-dDtemp*J2];
        Jac2=[sparse(1:nbus,1:nbus,J1*Ctemp+J2*Dtemp,nbus,nbus),sparse(1:nbus,1:nbus,J3*Ctemp+J4*Dtemp,nbus,nbus);...
            sparse(1:nbus,1:nbus,J3*Ctemp+J4*Dtemp,nbus,nbus),-sparse(1:nbus,1:nbus,J1*Ctemp+J2*Dtemp,nbus,nbus)];
        Jac3=[sparse(1:nbus,1:nbus,JNew.*Ctemp./Atemp,nbus,nbus),sparse(1:nbus,1:nbus,JNew.*Dtemp./Atemp,nbus,nbus);...
            sparse(1:nbus,1:nbus,KNew.*Ctemp./Atemp,nbus,nbus),sparse(1:nbus,1:nbus,KNew.*Dtemp./Atemp,nbus,nbus)];
        Jac4=[sparse(syn(:,1),syn(:,1),MatGRhsNew(:,1),nbus,nbus),sparse(syn(:,1),syn(:,1),MatGRhsNew(:,2),nbus,nbus);...
            sparse(syn(:,1),syn(:,1),MatGRhsNew(:,2),nbus,nbus),-sparse(syn(:,1),syn(:,1),MatGRhsNew(:,1),nbus,nbus)];
        Jac=Jac3+Jac4-Jac1-Jac2;
        
        JacVc=sparse(1:nbus,1:nbus,2*Ctemp,nbus,nbus);
        JacVd=sparse(1:nbus,1:nbus,2*Dtemp,nbus,nbus);
        
        Jac(nbus+ipv,:)=[JacVc(ipv,:),JacVd(ipv,:)];
        
        JacExt=[Jac,[P2freq;sparse(nbus,nbus)];...
            sparse(nbus,2*nbus),Freq2freq];
        
        corr=-JacExt(jacRowTag==1,jacColTag==1)\tDiff*stepRatio;
        corrFull=zeros(3*nbus,1);
        corrFull(jacColTag==1)=corr;
        Vtemp=Vtemp+corrFull(1:nbus)+1j*corrFull((nbus+1):(2*nbus));
        ftemp=ftemp+corrFull((2*nbus+1):(3*nbus));
        QTemp=imag(diff);
        %     QTemp(idxNonSw)
    end
elseif AEmethod==2 % NRJ - jacobian adjustment once
    for loop=1:maxLoop
        VtempOrig=Vtemp;
        ftempOrig=ftemp;
        QTempOrig=QTemp;
        % 1st step
        Ctemp=real(Vtemp);
        Dtemp=imag(Vtemp);
        Atemp=abs(Vtemp);
        Vgtemp=Vtemp(syn(:,1));
        Ig=accumarray(syn(:,1),MatGRhsNew(:,1)+1j*MatGRhsNew(:,2),[nbus,1])-...
            accumarray(syn(:,1),(MatGVNew(:,1).*real(Vgtemp)+MatGVNew(:,2).*imag(Vgtemp))+1j*(MatGVNew(:,3).*real(Vgtemp)+MatGVNew(:,4).*imag(Vgtemp)),[nbus,1]);
        diff=conj(SNew)+INew.*abs(Vtemp)+Ig.*conj(Vtemp)-(YNew*Vtemp).*conj(Vtemp);
        diffP=real(diff)+freqKeptTag.*(-fdk.*ftemp+dpg0);
        diffVQ=imag(diff);
        diffF=ftemp-ftemp(frefs(islands));
        if ~isempty(pv)
            diffVQ(pv(:,1))=abs(Vtemp(pv(:,1))).*abs(Vtemp(pv(:,1)))-sum(VspSq2(pv(:,1)),2);
        end
        tDiff=[diffP(idxNonSw);diffVQ(idxNonSw);diffF(freqKeptTagxRef==1)];
        if max(abs(tDiff))<diffTol
            QTemp=imag(diff);
            break;
        end
        dCtemp=sparse(1:nbus,1:nbus,Ctemp,nbus,nbus);
        dDtemp=sparse(1:nbus,1:nbus,Dtemp,nbus,nbus);
        Jac1=[dCtemp*J1+dDtemp*J3,dCtemp*J2+dDtemp*J4;...
            dCtemp*J3-dDtemp*J1,dCtemp*J4-dDtemp*J2];
        Jac2=[sparse(1:nbus,1:nbus,J1*Ctemp+J2*Dtemp,nbus,nbus),sparse(1:nbus,1:nbus,J3*Ctemp+J4*Dtemp,nbus,nbus);...
            sparse(1:nbus,1:nbus,J3*Ctemp+J4*Dtemp,nbus,nbus),-sparse(1:nbus,1:nbus,J1*Ctemp+J2*Dtemp,nbus,nbus)];
        Jac3=[sparse(1:nbus,1:nbus,JNew.*Ctemp./Atemp,nbus,nbus),sparse(1:nbus,1:nbus,JNew.*Dtemp./Atemp,nbus,nbus);...
            sparse(1:nbus,1:nbus,KNew.*Ctemp./Atemp,nbus,nbus),sparse(1:nbus,1:nbus,KNew.*Dtemp./Atemp,nbus,nbus)];
        Jac4=[sparse(syn(:,1),syn(:,1),MatGRhsNew(:,1),nbus,nbus),sparse(syn(:,1),syn(:,1),MatGRhsNew(:,2),nbus,nbus);...
            sparse(syn(:,1),syn(:,1),MatGRhsNew(:,2),nbus,nbus),-sparse(syn(:,1),syn(:,1),MatGRhsNew(:,1),nbus,nbus)];
        Jac=Jac3+Jac4-Jac1-Jac2;
        
        JacVc=sparse(1:nbus,1:nbus,2*Ctemp,nbus,nbus);
        JacVd=sparse(1:nbus,1:nbus,2*Dtemp,nbus,nbus);
        
        Jac(nbus+ipv,:)=[JacVc(ipv,:),JacVd(ipv,:)];
        
        JacExt=[Jac,[P2freq;sparse(nbus,nbus)];...
            sparse(nbus,2*nbus),Freq2freq];
        
        corr=-JacExt(jacRowTag==1,jacColTag==1)\tDiff*0.5*stepRatio;
        corrFull=zeros(3*nbus,1);
        corrFull(jacColTag==1)=corr;
        Vtemp=Vtemp+corrFull(1:nbus)+1j*corrFull((nbus+1):(2*nbus));
        ftemp=ftemp+corrFull((2*nbus+1):(3*nbus));
        QTemp=imag(diff);
        % 2nd step
        Ctemp=real(Vtemp);
        Dtemp=imag(Vtemp);
        Atemp=abs(Vtemp);
        Vgtemp=Vtemp(syn(:,1));
        Ig=accumarray(syn(:,1),MatGRhsNew(:,1)+1j*MatGRhsNew(:,2),[nbus,1])-...
            accumarray(syn(:,1),(MatGVNew(:,1).*real(Vgtemp)+MatGVNew(:,2).*imag(Vgtemp))+1j*(MatGVNew(:,3).*real(Vgtemp)+MatGVNew(:,4).*imag(Vgtemp)),[nbus,1]);
        diff=conj(SNew)+INew.*abs(Vtemp)+Ig.*conj(Vtemp)-(YNew*Vtemp).*conj(Vtemp);
        diffP=real(diff)+freqKeptTag.*(-fdk.*ftemp+dpg0);
        diffVQ=imag(diff);
        diffF=ftemp-ftemp(frefs(islands));
        if ~isempty(pv)
            diffVQ(pv(:,1))=abs(Vtemp(pv(:,1))).*abs(Vtemp(pv(:,1)))-sum(VspSq2(pv(:,1)),2);
        end
        tDiffx=[diffP(idxNonSw);diffVQ(idxNonSw);diffF(freqKeptTagxRef==1)];
        if max(abs(tDiffx))<diffTol
            QTemp=imag(diff);
            break;
        end
        dCtemp=sparse(1:nbus,1:nbus,Ctemp,nbus,nbus);
        dDtemp=sparse(1:nbus,1:nbus,Dtemp,nbus,nbus);
        Jac1=[dCtemp*J1+dDtemp*J3,dCtemp*J2+dDtemp*J4;...
            dCtemp*J3-dDtemp*J1,dCtemp*J4-dDtemp*J2];
        Jac2=[sparse(1:nbus,1:nbus,J1*Ctemp+J2*Dtemp,nbus,nbus),sparse(1:nbus,1:nbus,J3*Ctemp+J4*Dtemp,nbus,nbus);...
            sparse(1:nbus,1:nbus,J3*Ctemp+J4*Dtemp,nbus,nbus),-sparse(1:nbus,1:nbus,J1*Ctemp+J2*Dtemp,nbus,nbus)];
        Jac3=[sparse(1:nbus,1:nbus,JNew.*Ctemp./Atemp,nbus,nbus),sparse(1:nbus,1:nbus,JNew.*Dtemp./Atemp,nbus,nbus);...
            sparse(1:nbus,1:nbus,KNew.*Ctemp./Atemp,nbus,nbus),sparse(1:nbus,1:nbus,KNew.*Dtemp./Atemp,nbus,nbus)];
        Jac4=[sparse(syn(:,1),syn(:,1),MatGRhsNew(:,1),nbus,nbus),sparse(syn(:,1),syn(:,1),MatGRhsNew(:,2),nbus,nbus);...
            sparse(syn(:,1),syn(:,1),MatGRhsNew(:,2),nbus,nbus),-sparse(syn(:,1),syn(:,1),MatGRhsNew(:,1),nbus,nbus)];
        Jac=Jac3+Jac4-Jac1-Jac2;
        
        JacVc=sparse(1:nbus,1:nbus,2*Ctemp,nbus,nbus);
        JacVd=sparse(1:nbus,1:nbus,2*Dtemp,nbus,nbus);
        
        Jac(nbus+ipv,:)=[JacVc(ipv,:),JacVd(ipv,:)];
        
        JacExt=[Jac,[P2freq;sparse(nbus,nbus)];...
            sparse(nbus,2*nbus),Freq2freq];
        
        corr=-JacExt(jacRowTag==1,jacColTag==1)\tDiff*stepRatio;
        corrFull=zeros(3*nbus,1);
        corrFull(jacColTag==1)=corr;
        
        Vtemp=VtempOrig;
        ftemp=ftempOrig;
        QTemp=QTempOrig;
        
        Vtemp=Vtemp+corrFull(1:nbus)+1j*corrFull((nbus+1):(2*nbus));
        ftemp=ftemp+corrFull((2*nbus+1):(3*nbus));
        QTemp=imag(diff);
    end
elseif AEmethod==3 % NRM - 
    for loop=1:maxLoop
        Ctemp=real(Vtemp);
        Dtemp=imag(Vtemp);
        Atemp=abs(Vtemp);
        Vgtemp=Vtemp(syn(:,1));
        Ig=accumarray(syn(:,1),MatGRhsNew(:,1)+1j*MatGRhsNew(:,2),[nbus,1])-...
            accumarray(syn(:,1),(MatGVNew(:,1).*real(Vgtemp)+MatGVNew(:,2).*imag(Vgtemp))+1j*(MatGVNew(:,3).*real(Vgtemp)+MatGVNew(:,4).*imag(Vgtemp)),[nbus,1]);
        diff=conj(SNew)+INew.*abs(Vtemp)+Ig.*conj(Vtemp)-(YNew*Vtemp).*conj(Vtemp);
        diffP=real(diff)+freqKeptTag.*(-fdk.*ftemp+dpg0);
        diffVQ=imag(diff);
        diffF=ftemp-ftemp(frefs(islands));
        if ~isempty(pv)
            diffVQ(pv(:,1))=abs(Vtemp(pv(:,1))).*abs(Vtemp(pv(:,1)))-sum(VspSq2(pv(:,1)),2);
        end
        tDiff=[diffP(idxNonSw);diffVQ(idxNonSw);diffF(freqKeptTagxRef==1)];
        lambda=1e-3*norm(tDiff)^2;
        if max(abs(tDiff))<diffTol
            QTemp=imag(diff);
            break;
        end
        dCtemp=sparse(1:nbus,1:nbus,Ctemp,nbus,nbus);
        dDtemp=sparse(1:nbus,1:nbus,Dtemp,nbus,nbus);
        Jac1=[dCtemp*J1+dDtemp*J3,dCtemp*J2+dDtemp*J4;...
            dCtemp*J3-dDtemp*J1,dCtemp*J4-dDtemp*J2];
        Jac2=[sparse(1:nbus,1:nbus,J1*Ctemp+J2*Dtemp,nbus,nbus),sparse(1:nbus,1:nbus,J3*Ctemp+J4*Dtemp,nbus,nbus);...
            sparse(1:nbus,1:nbus,J3*Ctemp+J4*Dtemp,nbus,nbus),-sparse(1:nbus,1:nbus,J1*Ctemp+J2*Dtemp,nbus,nbus)];
        Jac3=[sparse(1:nbus,1:nbus,JNew.*Ctemp./Atemp,nbus,nbus),sparse(1:nbus,1:nbus,JNew.*Dtemp./Atemp,nbus,nbus);...
            sparse(1:nbus,1:nbus,KNew.*Ctemp./Atemp,nbus,nbus),sparse(1:nbus,1:nbus,KNew.*Dtemp./Atemp,nbus,nbus)];
        Jac4=[sparse(syn(:,1),syn(:,1),MatGRhsNew(:,1),nbus,nbus),sparse(syn(:,1),syn(:,1),MatGRhsNew(:,2),nbus,nbus);...
            sparse(syn(:,1),syn(:,1),MatGRhsNew(:,2),nbus,nbus),-sparse(syn(:,1),syn(:,1),MatGRhsNew(:,1),nbus,nbus)];
        Jac=Jac3+Jac4-Jac1-Jac2;
        
        JacVc=sparse(1:nbus,1:nbus,2*Ctemp,nbus,nbus);
        JacVd=sparse(1:nbus,1:nbus,2*Dtemp,nbus,nbus);
        
        Jac(nbus+ipv,:)=[JacVc(ipv,:),JacVd(ipv,:)];
        
        JacExt=[Jac,[P2freq;sparse(nbus,nbus)];...
            sparse(nbus,2*nbus),Freq2freq];
        
        filteredJac=JacExt(jacRowTag==1,jacColTag==1);        
        corr=-(filteredJac'*filteredJac+lambda*speye(size(filteredJac)))\(filteredJac'*tDiff)*stepRatio;
        corrFull=zeros(3*nbus,1);
        corrFull(jacColTag==1)=corr;
        Vtemp=Vtemp+corrFull(1:nbus)+1j*corrFull((nbus+1):(2*nbus));
        ftemp=ftemp+corrFull((2*nbus+1):(3*nbus));
        QTemp=imag(diff);
        %     QTemp(idxNonSw)
    end
elseif AEmethod==4 % NRMJ - 
    for loop=1:maxLoop
        VtempOrig=Vtemp;
        ftempOrig=ftemp;
        QTempOrig=QTemp;
        % 1st step
        Ctemp=real(Vtemp);
        Dtemp=imag(Vtemp);
        Atemp=abs(Vtemp);
        Vgtemp=Vtemp(syn(:,1));
        Ig=accumarray(syn(:,1),MatGRhsNew(:,1)+1j*MatGRhsNew(:,2),[nbus,1])-...
            accumarray(syn(:,1),(MatGVNew(:,1).*real(Vgtemp)+MatGVNew(:,2).*imag(Vgtemp))+1j*(MatGVNew(:,3).*real(Vgtemp)+MatGVNew(:,4).*imag(Vgtemp)),[nbus,1]);
        diff=conj(SNew)+INew.*abs(Vtemp)+Ig.*conj(Vtemp)-(YNew*Vtemp).*conj(Vtemp);
        diffP=real(diff)+freqKeptTag.*(-fdk.*ftemp+dpg0);
        diffVQ=imag(diff);
        diffF=ftemp-ftemp(frefs(islands));
        if ~isempty(pv)
            diffVQ(pv(:,1))=abs(Vtemp(pv(:,1))).*abs(Vtemp(pv(:,1)))-sum(VspSq2(pv(:,1)),2);
        end
        tDiff=[diffP(idxNonSw);diffVQ(idxNonSw);diffF(freqKeptTagxRef==1)];
        lambda=1e-3*norm(tDiff)^2;
        if max(abs(tDiff))<diffTol
            QTemp=imag(diff);
            break;
        end
        dCtemp=sparse(1:nbus,1:nbus,Ctemp,nbus,nbus);
        dDtemp=sparse(1:nbus,1:nbus,Dtemp,nbus,nbus);
        Jac1=[dCtemp*J1+dDtemp*J3,dCtemp*J2+dDtemp*J4;...
            dCtemp*J3-dDtemp*J1,dCtemp*J4-dDtemp*J2];
        Jac2=[sparse(1:nbus,1:nbus,J1*Ctemp+J2*Dtemp,nbus,nbus),sparse(1:nbus,1:nbus,J3*Ctemp+J4*Dtemp,nbus,nbus);...
            sparse(1:nbus,1:nbus,J3*Ctemp+J4*Dtemp,nbus,nbus),-sparse(1:nbus,1:nbus,J1*Ctemp+J2*Dtemp,nbus,nbus)];
        Jac3=[sparse(1:nbus,1:nbus,JNew.*Ctemp./Atemp,nbus,nbus),sparse(1:nbus,1:nbus,JNew.*Dtemp./Atemp,nbus,nbus);...
            sparse(1:nbus,1:nbus,KNew.*Ctemp./Atemp,nbus,nbus),sparse(1:nbus,1:nbus,KNew.*Dtemp./Atemp,nbus,nbus)];
        Jac4=[sparse(syn(:,1),syn(:,1),MatGRhsNew(:,1),nbus,nbus),sparse(syn(:,1),syn(:,1),MatGRhsNew(:,2),nbus,nbus);...
            sparse(syn(:,1),syn(:,1),MatGRhsNew(:,2),nbus,nbus),-sparse(syn(:,1),syn(:,1),MatGRhsNew(:,1),nbus,nbus)];
        Jac=Jac3+Jac4-Jac1-Jac2;
        
        JacVc=sparse(1:nbus,1:nbus,2*Ctemp,nbus,nbus);
        JacVd=sparse(1:nbus,1:nbus,2*Dtemp,nbus,nbus);
        
        Jac(nbus+ipv,:)=[JacVc(ipv,:),JacVd(ipv,:)];
        
        JacExt=[Jac,[P2freq;sparse(nbus,nbus)];...
            sparse(nbus,2*nbus),Freq2freq];
        
        filteredJac=JacExt(jacRowTag==1,jacColTag==1);        
        corr=-(filteredJac'*filteredJac+lambda*speye(size(filteredJac)))\(filteredJac'*tDiff)*0.5*stepRatio;
        corrFull=zeros(3*nbus,1);
        corrFull(jacColTag==1)=corr;
        Vtemp=Vtemp+corrFull(1:nbus)+1j*corrFull((nbus+1):(2*nbus));
        ftemp=ftemp+corrFull((2*nbus+1):(3*nbus));
        QTemp=imag(diff);
        % 2nd step
        Ctemp=real(Vtemp);
        Dtemp=imag(Vtemp);
        Atemp=abs(Vtemp);
        Vgtemp=Vtemp(syn(:,1));
        Ig=accumarray(syn(:,1),MatGRhsNew(:,1)+1j*MatGRhsNew(:,2),[nbus,1])-...
            accumarray(syn(:,1),(MatGVNew(:,1).*real(Vgtemp)+MatGVNew(:,2).*imag(Vgtemp))+1j*(MatGVNew(:,3).*real(Vgtemp)+MatGVNew(:,4).*imag(Vgtemp)),[nbus,1]);
        diff=conj(SNew)+INew.*abs(Vtemp)+Ig.*conj(Vtemp)-(YNew*Vtemp).*conj(Vtemp);
        diffP=real(diff)+freqKeptTag.*(-fdk.*ftemp+dpg0);
        diffVQ=imag(diff);
        diffF=ftemp-ftemp(frefs(islands));
        if ~isempty(pv)
            diffVQ(pv(:,1))=abs(Vtemp(pv(:,1))).*abs(Vtemp(pv(:,1)))-sum(VspSq2(pv(:,1)),2);
        end
        tDiffx=[diffP(idxNonSw);diffVQ(idxNonSw);diffF(freqKeptTagxRef==1)];
        lambdax=1e-3*norm(tDiffx)^2;
        if max(abs(tDiffx))<diffTol
            QTemp=imag(diff);
            break;
        end
        dCtemp=sparse(1:nbus,1:nbus,Ctemp,nbus,nbus);
        dDtemp=sparse(1:nbus,1:nbus,Dtemp,nbus,nbus);
        Jac1=[dCtemp*J1+dDtemp*J3,dCtemp*J2+dDtemp*J4;...
            dCtemp*J3-dDtemp*J1,dCtemp*J4-dDtemp*J2];
        Jac2=[sparse(1:nbus,1:nbus,J1*Ctemp+J2*Dtemp,nbus,nbus),sparse(1:nbus,1:nbus,J3*Ctemp+J4*Dtemp,nbus,nbus);...
            sparse(1:nbus,1:nbus,J3*Ctemp+J4*Dtemp,nbus,nbus),-sparse(1:nbus,1:nbus,J1*Ctemp+J2*Dtemp,nbus,nbus)];
        Jac3=[sparse(1:nbus,1:nbus,JNew.*Ctemp./Atemp,nbus,nbus),sparse(1:nbus,1:nbus,JNew.*Dtemp./Atemp,nbus,nbus);...
            sparse(1:nbus,1:nbus,KNew.*Ctemp./Atemp,nbus,nbus),sparse(1:nbus,1:nbus,KNew.*Dtemp./Atemp,nbus,nbus)];
        Jac4=[sparse(syn(:,1),syn(:,1),MatGRhsNew(:,1),nbus,nbus),sparse(syn(:,1),syn(:,1),MatGRhsNew(:,2),nbus,nbus);...
            sparse(syn(:,1),syn(:,1),MatGRhsNew(:,2),nbus,nbus),-sparse(syn(:,1),syn(:,1),MatGRhsNew(:,1),nbus,nbus)];
        Jac=Jac3+Jac4-Jac1-Jac2;
        
        JacVc=sparse(1:nbus,1:nbus,2*Ctemp,nbus,nbus);
        JacVd=sparse(1:nbus,1:nbus,2*Dtemp,nbus,nbus);
        
        Jac(nbus+ipv,:)=[JacVc(ipv,:),JacVd(ipv,:)];
        
        JacExt=[Jac,[P2freq;sparse(nbus,nbus)];...
            sparse(nbus,2*nbus),Freq2freq];
        
        filteredJac=JacExt(jacRowTag==1,jacColTag==1);        
        corr=-(filteredJac'*filteredJac+lambda*speye(size(filteredJac)))\(filteredJac'*tDiff)*stepRatio;
        corrFull=zeros(3*nbus,1);
        corrFull(jacColTag==1)=corr;
        
        Vtemp=VtempOrig;
        ftemp=ftempOrig;
        QTemp=QTempOrig;
        
        Vtemp=Vtemp+corrFull(1:nbus)+1j*corrFull((nbus+1):(2*nbus));
        ftemp=ftemp+corrFull((2*nbus+1):(3*nbus));
        QTemp=imag(diff);
    end
end

if loop>=maxLoop&&max(abs(tDiff))>=diffTol
    flag=-1;
else
    flag=0;
end

diffRec=max(abs(tDiff));

stateNew=xNew;
stateNew(idxs.vIdx)=Vtemp;
stateNew(idxs.qIdx)=QTemp;
stateNew(idxs.fIdx)=ftemp;

addLog(['[#END#] Solving algebraic equations - NR finished. Exitflag = ',num2str(flag), ', loop = ',num2str(loop),'.'],'INFO');

end