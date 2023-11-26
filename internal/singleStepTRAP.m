function [xTemp,diff,exitflag]=singleStepTRAP(SimData,SysData,SysPara,xt,t,dt,iden,AEmethod)
% Single step of trapezoidal formulation
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
%   xt - Initial system state
%	t - Relative time compared with SysPara starting time
%	dt - time step length
%	iden - Aux identifier
%	AEmethod - Method for solving algebraic equations
%		0 - HE
%		1 - NR
%
% OUTPUT
%	xTemp - State at the new step
%	diff - Error vector
%	exitflag - 
%		0 - Normally exit
%		-1 - Fail
%

[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
[nState,idxs]...
    =getIndexDyn(SysData);
[pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,~,~,~,~,Tmech1,Varref1,Ef1,Pm1,Eq11]=unfoldSysPara(SysPara);
[maxTime,segTime,~,nlvl,taylorN,alphaTol,diffTol,diffTolMax,method]=unfoldSimData(SimData);

if isfield(SysPara,'nIslands')&&isfield(SysPara,'islands')&&isfield(SysPara,'refs')
    nIslands=SysPara.nIslands;islands=SysPara.islands;refs=SysPara.refs;
else
    [nIslands,islands,refs]=searchIslands(bus(:,1),line(:,1:2));
end

nbus=size(bus,1);

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

% Determine the frequency model of each island
freqTypeTag=zeros(nIslands,1);%0:sw,1:syn,2:steady-state f
freqKeptTag=zeros(nbus,1);% corresponds to freqType=2
frefs=refs;
fswTag=zeros(nbus,1);
fsynTag=zeros(nbus,1);
fswTag(isw)=1;
fswTagxD=fswTag;
fsynTag(syn(:,1))=1;
D0=imag(xt(idxs.vIdx));
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
fIngIdx=((2*nbus+1):(3*nbus))';

stateFilterTag=ones(nState,1);
stateFilterTag(idxs.vIdx)=0;
stateFilterTag(idxs.qIdx)=0;
stateFilterTag(idxs.pgIdx)=0;
stateFilterTag(idxs.efIdx)=0;
stateFilterTag(idxs.tgovmIdx)=0;
stateFilterTag(idxs.fIdx)=0;
% stateFilterTag(dpgIdx)=0;
% busFilterTag=ones(nbus,1);
% busFilterTag(sw(:,1))=0;
% jacBusFilterTag=[busFilterTag;busFilterTag];

busTag=ones(nbus,1);
busTag(sw(:,1))=0;
jggRowTag=[busTag;busTag;freqKeptTagxRef];
jggColTag=[busTag;fswTagxD==0;freqKeptTag];

idxNonSw=find(busType~=2);
idxNonSwD=find(busType~=2&fswTagxD==1);

nInd=size(ind,1);
if ~isempty(ind)
    indIdx=ind(:,1);
    Rm1=ind(:,7);
    Xm1=ind(:,8);
    Zm1=ind(:,7)+1j*ind(:,8);
    Rm2=ind(:,9);
    Xm2=ind(:,10);
    Zmm=1j*ind(:,13);
    Tm0=ind(:,15)+ind(:,16)+ind(:,17);
    Tm1=-ind(:,16)-2*ind(:,17);
    Tm2=ind(:,17);
    Hm=ind(:,14);
end

nSyn=size(syn,1);
if ~isempty(syn)
    synIdx=syn(:,1);
    wgb=syn(:,4);
    model=syn(:,5);
    Xgl=syn(:,6);
    Rga=syn(:,7);
    Xgd=syn(:,8);
    Xgd1=syn(:,9);
    Xgd2=syn(:,10);
    Tgd1=syn(:,11);
    Tgd2=syn(:,12);
    Xgq=syn(:,13);
    Xgq1=syn(:,14);
    Xgq2=syn(:,15);
    Tgq1=syn(:,16);
    Tgq2=syn(:,17);
    Mg=syn(:,18);
    Dg=syn(:,19);
    TgAA=syn(:,24);
    gammad=Tgd2./Tgd1.*Xgd2./Xgd1.*(Xgd-Xgd1);
    gammaq=Tgq2./Tgq1.*Xgq2./Xgq1.*(Xgq-Xgq1);
end

nExc=size(exc,1);
if ~isempty(exc)
    excIdx=exc(:,1);
    VavrMax=exc(:,3);
    VavrMin=exc(:,4);
    muavr0=exc(:,5);
    Tavr1=exc(:,7);
    Tavr2=exc(:,6);
    vavrf0=exc(:,8);
    Vavr0=exc(:,9);
    Tavre=exc(:,10);
    Tavrr=exc(:,11);
end

nTg=size(tg,1);
if ~isempty(tg)
    tgIdx=tg(:,1);
    wtgref=tg(:,3);
    Rtg=tg(:,4);
    Ttgmax=tg(:,5);
    Ttgmin=tg(:,6);
    Ttg2=tg(:,7);
    Ttg1=tg(:,8);
end

[V0,Q0,s0,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,tgovmech0,f0,dpg0,qplt0,vg0]=unfoldX(xt,SysData);
if ~isempty(ind)
    [YshInd0,Yshind1]=getLinearInterpolatorInd(nbus,ind,s0,s0);
    Ysh0x=Ysh0+YshInd0;
    Ysh1x=Ysh1+Yshind1;
else
    Ysh0x=zeros(nbus,1);
    Ysh1x=zeros(nbus,1);
end
Ytr=Ytr0+Ytr1*t;
Ysh=Ysh0x+Ysh1x*t;
YtrNew=Ytr0+Ytr1*(t+dt);
YshNew=Ysh0x+Ysh1x*(t+dt);

SNew=-accumarray(pq(:,1),pq(:,4)+1j*pq(:,5),[nbus,1])+accumarray(pv(:,1),pv(:,4),[nbus,1]);
SNew=SNew-accumarray(pq(:,1),pqIncr(:,1)+1j*pqIncr(:,2),[nbus,1])*(t+dt)+accumarray(pv(:,1),pvIncr,[nbus,1])*(t+dt);
if size(pqIncr,2)>=4;SNew=SNew-accumarray(pq(:,1),pqIncr(:,3)+1j*pqIncr(:,4),[nbus,1])*(t+dt)*(t+dt);end
SNew=SNew-accumarray(zip(:,1),(Rzip0+Rzip1*(t+dt)).*(zip(:,7)+1j*zip(:,10)),[nbus,1]);

INew=-accumarray(zip(:,1),(Rzip0+Rzip1*(t+dt)).*(zip(:,6)-1j*zip(:,9)),[nbus,1]);
JNew=real(INew);
KNew=imag(INew);

YNew=YtrNew+sparse(1:nbus,1:nbus,YshNew,nbus,nbus);
GNew=real(YNew);BNew=imag(YNew);

% Calculate diff and jac, do N-R
Vtemp=V0;
pqx=pq;pqx(:,[4,5])=pq(:,[4,5])+pqIncr(:,1:2)*t;if size(pqIncr,2)>=4;pqx(:,[4,5])=pqx(:,[4,5])+pqIncr(:,3:4)*t*t;end
pvx=pv;pvx(:,4)=pv(:,4)+pvIncr*t;
indx=ind;if ~isempty(ind);indx(:,15:17)=ind(:,15:17).*repmat(Rind0+Rind1*t,1,3);indx(:,13)=ind(:,13)./(Reind0+Reind1*t);end
zipx=zip;if ~isempty(zip);zipx(:,5:10)=zip(:,5:10).*repmat(Rzip0+Rzip1*t,1,6);end
SysDatax=foldSysData(bus,sw,pvx,pqx,shunt,line,indx,zipx,syn,exc,tg,agc,cac,cluster);

pqxx=pq;pqxx(:,[4,5])=pq(:,[4,5])+pqIncr(:,1:2)*(t+dt);if size(pqIncr,2)>=4;pqxx(:,[4,5])=pqxx(:,[4,5])+pqIncr(:,3:4)*(t+dt)*(t+dt);end
pvxx=pv;pvxx(:,4)=pv(:,4)+pvIncr*(t+dt);
indxx=ind;if ~isempty(ind);indxx(:,15:17)=ind(:,15:17).*repmat(Rind0+Rind1*(t+dt),1,3);indxx(:,13)=ind(:,13)./(Reind0+Reind1*(t+dt));end
zipxx=zip;if ~isempty(zip);zipxx(:,5:10)=zip(:,5:10).*repmat(Rzip0+Rzip1*(t+dt),1,6);end
SysDataxx=foldSysData(bus,sw,pvxx,pqxx,shunt,line,indxx,zipxx,syn,exc,tg,agc,cac,cluster);

SysData=foldSysData(bus,sw,pvx,pqx,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster);
dx1=integ(SysDatax,SysPara,xt,dt);
xTemp=xt;
xTemp=xTemp+dx1;

iterMax=10;
exitflag=-1;
for iter=1:iterMax
    % calculate df and dg
    dx2=integ(SysDataxx,SysPara,xTemp,dt);
    diffF=xt-xTemp+(dx1+dx2)/2;
    diffF(stateFilterTag==0)=0;
    
    Ctemp=real(Vtemp);
    Dtemp=imag(Vtemp);
    Atemp=abs(Vtemp);
    Vgtemp=Vtemp(syn(:,1));
    ftemp=xTemp(idxs.fIdx);
    
    eq0=Ef0;ed0=zeros(size(Ef0));
    [~,~,sNew,dNew,wNew,eq1New,eq2New,ed1New,ed2New,psidNew,psiqNew,EdNew,EfNew,VavrmNew,VavrrNew,VavrfNew,VavrrefNew,tgovgNew,tgovmNew,tgovmechNew,fNew,dpgNew,qpltNew,vgNew]=unfoldX(xTemp,SysData);
    eqNew=EfNew;edNew=zeros(size(EfNew));
    
    [MatGV0,MatGV1,MatGRhs0,MatGRhs1]=getLinearInterpolatorSyn(syn,[],d0,dNew,ed0,ed10,ed20,edNew,ed1New,ed2New,eq0,eq10,eq20,eqNew,eq1New,eq2New,psid0,psiq0,psidNew,psiqNew);
    
    MatGVNew=MatGV0+MatGV1;
    MatGRhsNew=MatGRhs0+MatGRhs1;
    J1=GNew+sparse(syn(:,1),syn(:,1),MatGVNew(:,1),nbus,nbus);
    J2=-BNew+sparse(syn(:,1),syn(:,1),MatGVNew(:,2),nbus,nbus);
    J3=BNew+sparse(syn(:,1),syn(:,1),MatGVNew(:,3),nbus,nbus);
    J4=GNew+sparse(syn(:,1),syn(:,1),MatGVNew(:,4),nbus,nbus);
    P2freq=sparse(1:nbus,1:nbus,-freqKeptTag.*fdk,nbus,nbus);
    Freq2freq=sparse([1:nbus,1:nbus],[1:nbus,frefs(islands)'],[ones(1,nbus),-ones(1,nbus)],nbus,nbus);
    
    Ig=accumarray(syn(:,1),MatGRhsNew(:,1)+1j*MatGRhsNew(:,2),[nbus,1])-...
        accumarray(syn(:,1),(MatGVNew(:,1).*real(Vgtemp)+MatGVNew(:,2).*imag(Vgtemp))+1j*(MatGVNew(:,3).*real(Vgtemp)+MatGVNew(:,4).*imag(Vgtemp)),[nbus,1]);
%     diffG=conj(SNew)+INew.*abs(Vtemp)+Ig.*conj(Vtemp)-(YNew*Vtemp).*conj(Vtemp)+freqKeptTag.*(dpgNew-fdk.*fNew);
%     diffG(busFilterTag==0)=0;
%     diffGf=zeros(nbus,1);
    
    diff=conj(SNew)+INew.*abs(Vtemp)+Ig.*conj(Vtemp)-(YNew*Vtemp).*conj(Vtemp);
    diffP=real(diff)+freqKeptTag.*(-fdk.*ftemp+dpg0);
    diffVQ=imag(diff);
    diffGf=ftemp-ftemp(frefs(islands));
    if ~isempty(pv)
        diffVQ(pv(:,1))=abs(Vtemp(pv(:,1))).*abs(Vtemp(pv(:,1)))-(VspSq2(pv(:,1),1)+VspSq2(pv(:,1),2)*(t+dt));
    end
    
    diffGx=[diffP;diffVQ;diffGf];
    diffGx(jggRowTag==0)=0;
    
    if iter>1&&max(abs([diffF;diffGx]))<diffTol
        exitflag=0;
        break;
    end
    
    % construct Jgg
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
    
    Jgg=[Jac,[P2freq;sparse(nbus,nbus)];...
        sparse(nbus,2*nbus),Freq2freq];
    
    % construct Jgf
    
    jgfi=zeros(0,1);
    jgfj=zeros(0,1);
    jgfv=zeros(0,1);
    jfgi=zeros(0,1);
    jfgj=zeros(0,1);
    jfgv=zeros(0,1);
    jffi=zeros(0,1);
    jffj=zeros(0,1);
    jffv=zeros(0,1);
    
    if ~isempty(ind)
        % s
        deltasx=1e-3;
        Zind0=Zm1+Zmm.*(Rm2+1j*Xm2.*sNew)./(Rm2+(Zmm+1j*Xm2).*sNew);
        Yind0=1./Zind0;
        
        ZindNew=Zm1+Zmm.*(Rm2+1j*Xm2.*(sNew+deltasx))./(Rm2+(Zmm+1j*Xm2).*(sNew+deltasx));
        YindNew=1./ZindNew;
        
        Yshind1x=YindNew-Yind0;
        jgfi=[jgfi;indIdx;indIdx+nbus];
        jgfj=[jgfj;idxs.sIdx;idxs.sIdx];
        jgfv=[jgfv;Vtemp(indIdx).*conj(Vtemp(indIdx)).*real(Yshind1x/deltasx);Vtemp(indIdx).*conj(Vtemp(indIdx)).*imag(Yshind1x/deltasx)];
    end
    
    if ~isempty(syn)
        Vgtemp=Vtemp(synIdx);
        Cgtemp=real(Vgtemp);
        Dgtemp=imag(Vgtemp);
        
        cosd0=cos(dNew);
        sind0=sin(dNew);
        deltadNew=1e-5*abs(dNew)+1e-6;
        cosdx=cos(dNew+deltadNew);
        sindx=sin(dNew+deltadNew);
        
        MatGV0x=zeros(nSyn,4);
        MatGV1x=zeros(nSyn,4);
        
        MatGRhs0x=zeros(nSyn,2);
        MatGRhs1x=zeros(nSyn,2);
        
        for i=1:nSyn
            Mg0=[sind0(i),-cosd0(i);cosd0(i),sind0(i)];
            Mgx=[sindx(i),-cosdx(i);cosdx(i),sindx(i)];
            if model(i)==6||model(i)==5
                Yg=[Rga(i),Xgq2(i);-Xgd2(i),Rga(i)]/(Rga(i)*Rga(i)+Xgd2(i)*Xgq2(i));
                Mxye=Mg0'*Yg;
                Mdqe=Yg;
                dIdq=-Yg*(Mgx-Mg0)*[Cgtemp(i);Dgtemp(i)]/deltadNew(i);
                dVdq=-Yg*Mg0;
                jgfi=[jgfi;synIdx(i);synIdx(i);synIdx(i)+nbus;synIdx(i)+nbus];
                jgfj=[jgfj;idxs.ed2Idx(i);idxs.eq2Idx(i);idxs.ed2Idx(i);idxs.eq2Idx(i)];
                jgfv=[jgfv;Cgtemp(i)*Mxye(1,1)+Dgtemp(i)*Mxye(2,1);Cgtemp(i)*Mxye(1,2)+Dgtemp(i)*Mxye(2,2);Cgtemp(i)*Mxye(2,1)-Dgtemp(i)*Mxye(1,1);Cgtemp(i)*Mxye(2,2)-Dgtemp(i)*Mxye(1,2)];
                
                jffi=[jffi;idxs.eq2Idx(i);idxs.eq2Idx(i);idxs.eq2Idx(i)];
                jffj=[jffj;idxs.ed2Idx(i);idxs.eq2Idx(i);idxs.deltaIdx(i)];
                jffv=[jffv;-(Xgd1(i)-Xgd2(i)+gammad(i))*Mdqe(1,1)/Tgd2(i);-(Xgd1(i)-Xgd2(i)+gammad(i))*Mdqe(1,2)/Tgd2(i);-(Xgd1(i)-Xgd2(i)+gammad(i))*dIdq(1)/Tgd2(i)];
                jffi=[jffi;idxs.ed2Idx(i);idxs.ed2Idx(i);idxs.ed2Idx(i)];
                jffj=[jffj;idxs.ed2Idx(i);idxs.eq2Idx(i);idxs.deltaIdx(i)];
                jffv=[jffv;(Xgq1(i)-Xgq2(i)+gammaq(i))*Mdqe(2,1)/Tgq2(i);(Xgq1(i)-Xgq2(i)+gammaq(i))*Mdqe(2,2)/Tgq2(i);(Xgq1(i)-Xgq2(i)+gammaq(i))*dIdq(2)/Tgq2(i)];
                jffi=[jffi;idxs.eq1Idx(i);idxs.eq1Idx(i);idxs.eq1Idx(i)];
                jffj=[jffj;idxs.ed2Idx(i);idxs.eq2Idx(i);idxs.deltaIdx(i)];
                jffv=[jffv;-(Xgd(i)-Xgd1(i)-gammad(i))*Mdqe(1,1)/Tgd1(i);-(Xgd(i)-Xgd1(i)-gammad(i))*Mdqe(1,2)/Tgd1(i);-(Xgd(i)-Xgd1(i)-gammad(i))*dIdq(1)/Tgd1(i)];
                if model(i)==6
                    jffi=[jffi;idxs.ed1Idx(i);idxs.ed1Idx(i);idxs.ed1Idx(i)];
                    jffj=[jffj;idxs.ed2Idx(i);idxs.eq2Idx(i);idxs.deltaIdx(i)];
                    jffv=[jffv;(Xgq(i)-Xgq1(i)-gammaq(i))*Mdqe(2,1)/Tgq1(i);(Xgq(i)-Xgq1(i)-gammaq(i))*Mdqe(2,2)/Tgq1(i);(Xgq(i)-Xgq1(i)-gammaq(i))*dIdq(2)/Tgq1(i)];
                end
                
                jfgi=[jfgi;idxs.eq2Idx(i);idxs.eq2Idx(i)];
                jfgj=[jfgj;synIdx(i);synIdx(i)+nbus];
                jfgv=[jfgv;-(Xgd1(i)-Xgd2(i)+gammad(i))*dVdq(1,1)/Tgd2(i);-(Xgd1(i)-Xgd2(i)+gammad(i))*dVdq(1,2)/Tgd2(i)];
                jfgi=[jfgi;idxs.ed2Idx(i);idxs.ed2Idx(i)];
                jfgj=[jfgj;synIdx(i);synIdx(i)+nbus];
                jfgv=[jfgv;(Xgq1(i)-Xgq2(i)+gammaq(i))*dVdq(2,1)/Tgq2(i);(Xgq1(i)-Xgq2(i)+gammaq(i))*dVdq(2,2)/Tgq2(i)];
                jfgi=[jfgi;idxs.eq1Idx(i);idxs.eq1Idx(i)];
                jfgj=[jfgj;synIdx(i);synIdx(i)+nbus];
                jfgv=[jfgv;-(Xgd(i)-Xgd1(i)-gammad(i))*dVdq(1,1)/Tgd1(i);-(Xgd(i)-Xgd1(i)-gammad(i))*dVdq(1,2)/Tgd1(i)];
                if model(i)==6
                    jfgi=[jfgi;idxs.ed1Idx(i);idxs.ed1Idx(i)];
                    jfgj=[jfgj;synIdx(i);synIdx(i)+nbus];
                    jfgv=[jfgv;(Xgq(i)-Xgq1(i)-gammaq(i))*Mdqe(2,1)/Tgq1(i);(Xgq(i)-Xgq1(i)-gammaq(i))*Mdqe(2,2)/Tgq1(i)];
                end
                
                MatGRhs0x(i,:)=(Mg0'*Yg*[ed2New(i);eq2New(i)]).';
                MatGRhs1x(i,:)=(Mgx'*Yg*[ed2New(i);eq2New(i)]).'-MatGRhs0x(i,:);
            elseif model(i)==4
                Yg=[Rga(i),Xgq1(i);-Xgd1(i),Rga(i)]/(Rga(i)*Rga(i)+Xgd1(i)*Xgq1(i));
                Mxye=Mg0'*Yg;
                Mdqe=Yg;
                dVdq=-Yg*Mg0;
                dIdq=-Yg*(Mgx-Mg0)*[Cgtemp(i);Dgtemp(i)]/deltadNew(i);
                jgfi=[jgfi;synIdx(i);synIdx(i);synIdx(i)+nbus;synIdx(i)+nbus];
                jgfj=[jgfj;idxs.ed1Idx(i);idxs.eq1Idx(i);idxs.ed1Idx(i);idxs.eq1Idx(i)];
                jgfv=[jgfv;Cgtemp(i)*Mxye(1,1)+Dgtemp(i)*Mxye(2,1);Cgtemp(i)*Mxye(1,2)+Dgtemp(i)*Mxye(2,2);Cgtemp(i)*Mxye(2,1)-Dgtemp(i)*Mxye(1,1);Cgtemp(i)*Mxye(2,2)-Dgtemp(i)*Mxye(1,2)];
                
                jffi=[jffi;idxs.eq1Idx(i);idxs.eq1Idx(i);idxs.eq1Idx(i)];
                jffj=[jffj;idxs.ed1Idx(i);idxs.eq1Idx(i);idxs.deltaIdx(i)];
                jffv=[jffv;-(Xgd(i)-Xgd1(i))*Mdqe(1,1)/Tgd1(i);-(Xgd(i)-Xgd1(i))*Mdqe(1,2)/Tgd1(i);-(Xgd(i)-Xgd1(i))*dIdq(1)/Tgd1(i)];
                jffi=[jffi;idxs.ed1Idx(i);idxs.ed1Idx(i);idxs.ed1Idx(i)];
                jffj=[jffj;idxs.ed1Idx(i);idxs.eq1Idx(i);idxs.deltaIdx(i)];
                jffv=[jffv;(Xgq(i)-Xgq1(i))*Mdqe(2,1)/Tgq1(i);(Xgq(i)-Xgq1(i))*Mdqe(2,2)/Tgq1(i);(Xgq(i)-Xgq1(i))*dIdq(2)/Tgq1(i)];
                
                jfgi=[jfgi;idxs.eq1Idx(i);idxs.eq1Idx(i)];
                jfgj=[jfgj;synIdx(i);synIdx(i)+nbus];
                jfgv=[jfgv;-(Xgd(i)-Xgd1(i))*dVdq(1,1)/Tgd1(i);-(Xgd(i)-Xgd1(i))*dVdq(1,2)/Tgd1(i)];
                jfgi=[jfgi;idxs.ed1Idx(i);idxs.ed1Idx(i)];
                jfgj=[jfgj;synIdx(i);synIdx(i)+nbus];
                jfgv=[jfgv;(Xgq(i)-Xgq1(i))*Mdqe(2,1)/Tgq1(i);(Xgq(i)-Xgq1(i))*Mdqe(2,2)/Tgq1(i)];
                
                MatGRhs0x(i,:)=(Mg0'*Yg*[ed1New(i);eq1New(i)]).';
                MatGRhs1x(i,:)=(Mgx'*Yg*[ed1New(i);eq1New(i)]).'-MatGRhs0x(i,:);
            elseif model(i)==3||model(i)==2
                Yg=[Rga(i),Xgq(i);-Xgd1(i),Rga(i)]/(Rga(i)*Rga(i)+Xgd1(i)*Xgq(i));
                Mxye=Mg0'*Yg;
                Mdqe=Yg;
                dVdq=-Yg*Mg0;
                dIdq=-Yg*(Mgx-Mg0)*[Cgtemp(i);Dgtemp(i)]/deltadNew(i);
                if model(i)==3
                    jgfi=[jgfi;synIdx(i);synIdx(i)+nbus];
                    jgfj=[jgfj;idxs.eq1Idx(i);idxs.eq1Idx(i)];
                    jgfv=[jgfv;Cgtemp(i)*Mxye(1,2)+Dgtemp(i)*Mxye(2,2);Cgtemp(i)*Mxye(2,2)-Dgtemp(i)*Mxye(1,2)];
                    jffi=[jffi;idxs.eq1Idx(i);idxs.eq1Idx(i)];
                    jffj=[jffj;idxs.eq1Idx(i);idxs.deltaIdx(i)];
                    jffv=[jffv;-(Xgd(i)-Xgd1(i))*Mdqe(1,1)/Tgd1(i);-(Xgd(i)-Xgd1(i))*dIdq(1)/Tgd1(i)];
                    jfgi=[jfgi;idxs.eq1Idx(i);idxs.eq1Idx(i)];
                    jfgj=[jfgj;synIdx(i);synIdx(i)+nbus];
                    jfgv=[jfgv;-(Xgd(i)-Xgd1(i))*dVdq(1,1)/Tgd1(i);-(Xgd(i)-Xgd1(i))*dVdq(1,2)/Tgd1(i)];
                end
                
                MatGRhs0x(i,:)=(Mg0'*Yg*[edNew(i);eq1New(i)]).';
                MatGRhs1x(i,:)=(Mgx'*Yg*[edNew(i);eq1New(i)]).'-MatGRhs0x(i,:);
            elseif model(i)==8
                Yg=zeros(2);
                Mxye=Mg0'*diag([1/Xgd2(i);1/Xgq2(i)]);
                Mdqe=diag([1/Xgd2(i);1/Xgq2(i)]);
                dVdq=-Yg*Mg0;
                dIdq=-Yg*(Mgx-Mg0)*[Cgtemp(i);Dgtemp(i)]/deltadNew(i);
                jgfi=[jgfi;synIdx(i);synIdx(i);synIdx(i)+nbus;synIdx(i)+nbus];
                jgfj=[jgfj;idxs.ed2Idx(i);idxs.eq2Idx(i);idxs.ed2Idx(i);idxs.eq2Idx(i)];
                jgfv=[jgfv;-(Cgtemp(i)*Mxye(1,1)+Dgtemp(i)*Mxye(2,1));Cgtemp(i)*Mxye(1,2)+Dgtemp(i)*Mxye(2,2);-(Cgtemp(i)*Mxye(2,1)-Dgtemp(i)*Mxye(1,1));Cgtemp(i)*Mxye(2,2)-Dgtemp(i)*Mxye(1,2)];
                jgfi=[jgfi;synIdx(i);synIdx(i);synIdx(i)+nbus;synIdx(i)+nbus];
                jgfj=[jgfj;idxs.psiqIdx(i);idxs.psidIdx(i);idxs.psiqIdx(i);idxs.psidIdx(i)];
                jgfv=[jgfv;-(Cgtemp(i)*Mxye(1,1)+Dgtemp(i)*Mxye(2,1));-(Cgtemp(i)*Mxye(1,2)+Dgtemp(i)*Mxye(2,2));-(Cgtemp(i)*Mxye(2,1)-Dgtemp(i)*Mxye(1,1));-(Cgtemp(i)*Mxye(2,2)-Dgtemp(i)*Mxye(1,2))];
                
                jffi=[jffi;idxs.eq2Idx(i);idxs.eq2Idx(i);idxs.eq2Idx(i);idxs.eq2Idx(i);idxs.eq2Idx(i)];
                jffj=[jffj;idxs.ed2Idx(i);idxs.eq2Idx(i);idxs.psiqIdx(i);idxs.psidIdx(i);idxs.deltaIdx(i)];
                jffv=[jffv;-(Xgd1(i)-Xgd2(i)+gammad(i))*Mdqe(1,1)/Tgd2(i);-(Xgd1(i)-Xgd2(i)+gammad(i))*Mdqe(1,2)/Tgd2(i);...
                    -(Xgd1(i)-Xgd2(i)+gammad(i))*Mdqe(1,1)/Tgd2(i);(Xgd1(i)-Xgd2(i)+gammad(i))*Mdqe(1,2)/Tgd2(i);-(Xgd1(i)-Xgd2(i)+gammad(i))*dIdq(1)/Tgd2(i)];
                jffi=[jffi;idxs.ed2Idx(i);idxs.ed2Idx(i);idxs.ed2Idx(i);idxs.ed2Idx(i);idxs.ed2Idx(i)];
                jffj=[jffj;idxs.ed2Idx(i);idxs.eq2Idx(i);idxs.psiqIdx(i);idxs.psidIdx(i);idxs.deltaIdx(i)];
                jffv=[jffv;(Xgq1(i)-Xgq2(i)+gammaq(i))*Mdqe(2,1)/Tgq2(i);(Xgq1(i)-Xgq2(i)+gammaq(i))*Mdqe(2,2)/Tgq2(i);...
                    (Xgq1(i)-Xgq2(i)+gammaq(i))*Mdqe(2,1)/Tgq2(i);-(Xgq1(i)-Xgq2(i)+gammaq(i))*Mdqe(2,2)/Tgq2(i);(Xgq1(i)-Xgq2(i)+gammaq(i))*dIdq(2)/Tgq2(i)];
                jffi=[jffi;idxs.eq1Idx(i);idxs.eq1Idx(i);idxs.eq1Idx(i);idxs.eq1Idx(i);idxs.eq1Idx(i)];
                jffj=[jffj;idxs.ed2Idx(i);idxs.eq2Idx(i);idxs.psiqIdx(i);idxs.psidIdx(i);idxs.deltaIdx(i)];
                jffv=[jffv;-(Xgd(i)-Xgd1(i)-gammad(i))*Mdqe(1,1)/Tgd1(i);-(Xgd(i)-Xgd1(i)-gammad(i))*Mdqe(1,2)/Tgd1(i);...
                    -(Xgd(i)-Xgd1(i)-gammad(i))*Mdqe(1,1)/Tgd1(i);(Xgd(i)-Xgd1(i)-gammad(i))*Mdqe(1,2)/Tgd1(i);-(Xgd(i)-Xgd1(i)-gammad(i))*dIdq(1)/Tgd1(i)];
                jffi=[jffi;idxs.ed1Idx(i);idxs.ed1Idx(i);idxs.ed1Idx(i);idxs.ed1Idx(i);idxs.ed1Idx(i)];
                jffj=[jffj;idxs.ed2Idx(i);idxs.eq2Idx(i);idxs.psiqIdx(i);idxs.psidIdx(i);idxs.deltaIdx(i)];
                jffv=[jffv;(Xgq(i)-Xgq1(i)-gammaq(i))*Mdqe(2,1)/Tgq1(i);(Xgq(i)-Xgq1(i)-gammaq(i))*Mdqe(2,2)/Tgq1(i);...
                    (Xgq(i)-Xgq1(i)-gammaq(i))*Mdqe(2,1)/Tgq1(i);-(Xgq(i)-Xgq1(i)-gammaq(i))*Mdqe(2,2)/Tgq1(i);(Xgq(i)-Xgq1(i)-gammaq(i))*dIdq(2)/Tgq1(i)];
                
                MatGRhs0x(i,:)=(Mg0'*[(eq2New(i)-psidNew(i))/Xgd2(i);(-ed2New(i)-psiqNew(i))/Xgq2(i)]).';
                MatGRhs1x(i,:)=(Mgx'*[(eq2New(i)-psidNew(i))/Xgd2(i);(-ed2New(i)-psiqNew(i))/Xgq2(i)]).'-MatGRhs0x(i,:);
            end
            MatVg0x=Mg0'*Yg*Mg0;
            MatVgNewx=Mgx'*Yg*Mgx;
            MatGV0x(i,:)=[MatVg0x(1,1),MatVg0x(1,2),MatVg0x(2,1),MatVg0x(2,2)];
            MatGV1x(i,:)=[MatVgNewx(1,1),MatVgNewx(1,2),MatVgNewx(2,1),MatVgNewx(2,2)]-MatGV0x(i,:);
            
            jgfi=[jgfi;synIdx(i);synIdx(i)+nbus];
            jgfj=[jgfj;idxs.deltaIdx(i);idxs.deltaIdx(i)];
            jgfv=[jgfv;...
                (Cgtemp(i)*(MatGRhs1x(i,1)-MatGV1x(i,1)*Cgtemp(i)-MatGV1x(i,2)*Dgtemp(i))+Dgtemp(i)*(MatGRhs1x(i,2)-MatGV1x(i,3)*Cgtemp(i)-MatGV1x(i,4)*Dgtemp(i)))/deltadNew(i);...
                (-Dgtemp(i)*(MatGRhs1x(i,1)-MatGV1x(i,1)*Cgtemp(i)-MatGV1x(i,2)*Dgtemp(i))+Cgtemp(i)*(MatGRhs1x(i,2)-MatGV1x(i,3)*Cgtemp(i)-MatGV1x(i,4)*Dgtemp(i)))/deltadNew(i)];
        end
        
    end
    Jgf=sparse(jgfi,jgfj,jgfv,3*nbus,nState);
    
    % construct Jff & Jfg
    
    if ~isempty(ind)
        % s
        jffi=[jffi;idxs.sIdx];
        jffj=[jffj;idxs.sIdx];
        divs=((Rm2+Rm1.*sNew).*(Rm2+Rm1.*sNew)+(Xm2+Xm1).*(Xm2+Xm1).*sNew.*sNew);
        jffv=[jffv;(Tm1+2*Tm2.*sNew+Vtemp(indIdx).*conj(Vtemp(indIdx)).*Rm2.*(1./divs-2*sNew.*(sNew.*(Xm2+Xm1).*(Xm2+Xm1)+Rm1.*(Rm2+Rm1.*sNew))./divs./divs))/2./Hm];
        %             il=V(indIdx).*yind;
        %             ir=il-(V(indIdx)-il.*Zm1)./Zmm;
        %             dx(sIdx)=dt*((Tm0+Tm1.*s+Tm2.*s.*s)-ir.*conj(ir).*Rm2./s)/2./Hm;
        
        yind=1./(Zm1+Zmm.*(Rm2+1j*Xm2.*sNew)./(Rm2+sNew.*(Zmm+1j*Xm2)));
        yInds=(1+Zm1./Zmm).*yind-1./Zmm;
        yIndsRs=yInds.*conj(yInds).*Rm2./sNew;
        yIndsRs(isnan(yIndsRs))=0;
        jfgi=[jfgi;idxs.sIdx;idxs.sIdx];
        jfgj=[jfgj;indIdx;indIdx+nbus];
        jfgv=[jfgv;2*real(Vtemp(indIdx)).*yIndsRs/2./Hm;2*imag(Vtemp(indIdx)).*yIndsRs/2./Hm];
    end
    
    if ~isempty(syn)
       if ~isempty(tg)
        pActiveInSyn=ones(nSyn,1);        
        tgovmActive=ones(nTg,1);
        tgovmActive((tgovm0-Ttgmax)>0)=0;
        tgovmActive((tgovm0-Ttgmin)<0)=0;
        pActiveInSyn(tgIdx)=tgovmActive;
%         else
%             tgIdx=[];
%              tgovmActive=zeros(nTg,1);
%              pActiveInSyn=zeros(nTg,1);
       end
       if ~isempty(exc)
        texcActive=ones(nExc,1);
        texcActive((Vavrf0-VavrMax)>0)=0;
        texcActive((Vavrf0-VavrMin)<0)=0;
        texcActiveInSyn=zeros(nSyn,1);
        texcActiveInSyn(tgIdx)=texcActive;
       end
        numSynOnBus=accumarray(syn(:,1),1,[nbus,1]);
        %delta
        jffi=[jffi;idxs.deltaIdx];
        jffj=[jffj;idxs.omegaIdx];
        jffv=[jffv;wgb];
        %omega
        jffi=[jffi;idxs.omegaIdx;idxs.omegaIdx(tgIdx);idxs.omegaIdx(tgIdx);idxs.omegaIdx(tgIdx);idxs.omegaIdx];
        jffj=[jffj;idxs.omegaIdx;idxs.tgovgIdx;idxs.omegaIdx(tgIdx);idxs.tmechIdx;idxs.dpgIdx(synIdx)];
        jffv=[jffv;-Dg/2./Mg;tgovmActive;-tgovmActive.*Ttg1./Ttg2./Rtg;tgovmActive;pActiveInSyn./numSynOnBus(synIdx)];
        %psiq
        jffi=[jffi;idxs.psidIdx(model>=8)];
        jffj=[jffj;idxs.psiqIdx(model>=8)];
        jffv=[jffv;wgb(model>=8)];
        %psid
        jffi=[jffi;idxs.psiqIdx(model>=8)];
        jffj=[jffj;idxs.psidIdx(model>=8)];
        jffv=[jffv;-wgb(model>=8)];
        %ed''
        jffi=[jffi;idxs.ed2Idx(model>=5);idxs.ed2Idx(model>=5)];
        jffj=[jffj;idxs.ed2Idx(model>=5);idxs.ed1Idx(model>=5)];
        jffv=[jffv;-1./Tgq2(model>=5);1./Tgq2(model>=5)];
        %eq''
        jffi=[jffi;idxs.eq2Idx(model>=5);idxs.eq2Idx(model>=5);idxs.eq2Idx(model>=5&texcActiveInSyn==1)];
        jffj=[jffj;idxs.eq2Idx(model>=5);idxs.eq1Idx(model>=5);idxs.vavrfIdx(texcActive==1&model(excIdx)>=5)];
        jffv=[jffv;-1./Tgd2(model>=5);1./Tgd2(model>=5);TgAA(model>=5&texcActiveInSyn==1)./Tgd1(model>=5&texcActiveInSyn==1)./Tgd2(model>=5&texcActiveInSyn==1)];
        %ed'
        jffi=[jffi;idxs.ed1Idx(model>=4)];
        jffj=[jffj;idxs.ed1Idx(model>=4)];
        jffv=[jffv;-1./Tgq1(model>=4)];
        %eq'
        jffi=[jffi;idxs.eq1Idx(model>=3)];
        jffj=[jffj;idxs.eq1Idx(model>=3)];
        jffv=[jffv;-1./Tgd1(model>=3)];
        jffi=[jffi;idxs.eq1Idx(model>=5&texcActiveInSyn==1)];
        jffj=[jffj;idxs.vavrfIdx(texcActive==1&model(excIdx)>=5)];
        jffv=[jffv;(1-TgAA(model>=5&texcActiveInSyn==1)./Tgd1(model>=5&texcActiveInSyn==1))./Tgd2(model>=5&texcActiveInSyn==1)];
        jffi=[jffi;idxs.eq1Idx((model==3|model==4)&texcActiveInSyn==1)];
        jffj=[jffj;idxs.vavrfIdx(texcActive==1&(model(excIdx)==3|model(excIdx)==4))];
        jffv=[jffv;(1-TgAA((model==3|model==4)&texcActiveInSyn==1)./Tgd1((model==3|model==4)&texcActiveInSyn==1))./Tgd2((model==3|model==4)&texcActiveInSyn==1)];
    end
    
    if ~isempty(tg)
        jffi=[jffi;idxs.tgovgIdx;idxs.tgovgIdx];
        jffj=[jffj;idxs.omegaIdx;idxs.tgovgIdx];
        jffv=[jffv;-(1-Ttg1./Ttg2)./Rtg./Ttg2;-1./Ttg2];
    end
    
    if ~isempty(exc)
        VmagAvr=abs(Vtemp(synIdx(excIdx)));
        
        jffi=[jffi;idxs.vavrmIdx];
        jffj=[jffj;idxs.vavrmIdx];
        jffv=[jffv;-1./Tavrr];
        
        jffi=[jffi;idxs.vavrrIdx;idxs.vavrrIdx;idxs.vavrrIdx];
        jffj=[jffj;idxs.vavrmIdx;idxs.vavrrefIdx;idxs.vavrrIdx];
        jffv=[jffv;-muavr0.*(1-Tavr1./Tavr2)./Tavr2;muavr0.*(1-Tavr1./Tavr2)./Tavr2;-1./Tavr2];
        
        jffi=[jffi;idxs.vavrfIdx;idxs.vavrfIdx;idxs.vavrfIdx;idxs.vavrfIdx];
        jffj=[jffj;idxs.vavrmIdx;idxs.vavrrefIdx;idxs.vavrrIdx;idxs.vavrfIdx];
        jffv=[jffv;-muavr0.*Tavr1./Tavr2./Tavre.*VmagAvr./Vavr0;muavr0.*Tavr1./Tavr2./Tavre.*VmagAvr./Vavr0;VmagAvr./Vavr0./Tavre;-1./Tavre];
    end
    
    if ~isempty(agc)
        agcExt=zeros(nbus,size(agc,2));
        agcExt(agc(:,1),:)=agc;
        
        jfgi=[jfgi;idxs.dpgIdx(freqKeptTag==1)];
        jfgj=[jfgj;idxs.fIngIdx(freqKeptTag==1)];%f
        jfgv=[jfgv;-agcExt(freqKeptTag==1,4)];
        
        numSynOnBus=accumarray(syn(:,1),1,[nbus,1]);
        jffi=[jffi;idxs.dpgIdx(syn(:,1))];
        jffj=[jffj;idxs.omegaIdx];%f
        jffv=[jffv;-agcExt(syn(:,1),4)./numSynOnBus(syn(:,1))];
    end
    
    Jff=dt*sparse(jffi,jffj,jffv,nState,nState)/2-speye(nState,nState);
    
    % construct Jfg
    Jfg=dt*sparse(jfgi,jfgj,jfgv,nState,3*nbus)/2; %P,Q,f
    
    jac=[Jff(stateFilterTag==1,stateFilterTag==1),Jfg(stateFilterTag==1,jggColTag==1);Jgf(jggRowTag==1,stateFilterTag==1),Jgg(jggRowTag==1,jggColTag==1)];
    
    correction=-jac\[diffF(stateFilterTag==1);diffGx(jggRowTag==1)];
    nStateCorr=sum(stateFilterTag);
    xTemp(stateFilterTag==1)=xTemp(stateFilterTag==1)+correction(1:nStateCorr);
    corrFull=zeros(3*nbus,1);
    corrFull(jggRowTag==1)=correction((nStateCorr+1):end);
    Vtemp=Vtemp+corrFull(1:nbus)+1j*corrFull((nbus+1):(2*nbus));
    xTemp(idxs.vIdx)=Vtemp;
    xTemp(idxs.fIdx)=corrFull(fIngIdx);
    xTemp=adjustAlgebraic(SysDataxx,xTemp);
end

dxc=(xTemp-xt)/dt;
SysParaxx=foldSysPara([],[],[],[],[],[],[],[],Ytr0+Ytr1*(t+dt),[],Ysh0+Ysh1*(t+dt),[],[VspSq2(:,1)+VspSq2(:,2)*(t+dt),zeros(size(VspSq2,1),1)],[],[],[],[]);
diff=checkEquationBalanceSynDyn(SysDataxx,SysParaxx,xTemp,dxc);

end
% % function [xTemp,diff,exitflag]=singleStepTRAP(SimData,SysData,SysPara,xt,t,dt,iden,AEmethod)
% % % Single step of trapezoidal formulation
% % % Copyright (C) Rui Yao <yaorui.thu@gmail.com>
% % %
% % % INPUT
% % %   SimData - Simulation parameters
% % %   SysData - System data for simulation
% % %   SysPara - Parameters representing the events happening in the system
% % %   xt - Initial system state
% % %	t - Relative time compared with SysPara starting time
% % %	dt - time step length
% % %	iden - Aux identifier
% % %	AEmethod - Method for solving algebraic equations
% % %		0 - HE
% % %		1 - NR
% % %
% % % OUTPUT
% % %	xTemp - State at the new step
% % %	diff - Error vector
% % %	exitflag - 
% % %		0 - Normally exit
% % %		-1 - Fail
% % %
% % 
% % [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
% % % [nState,vIdx,qIdx,sIdx,deltaIdx,omegaIdx,eq1Idx,eq2Idx,ed1Idx,ed2Idx,psidIdx,psiqIdx,pgIdx,efIdx,...
% % %     vavrmIdx,vavrrIdx,vavrfIdx,vavrrefIdx,tgovgIdx,tgovmIdx,tmechIdx,fIdx,dpgIdx,qpltIdx,vgIdx]...
% % %     =getIndexDyn(SysData);
% % [nState,idxs]...
% %     =getIndexDyn(SysData);
% % [pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,~,~,~,~,Tmech1,Varref1,Ef1,Pm1,Eq11]=unfoldSysPara(SysPara);
% % [maxTime,segTime,~,nlvl,taylorN,alphaTol,diffTol,diffTolMax,method]=unfoldSimData(SimData);
% % 
% % if isfield(SysPara,'nIslands')&&isfield(SysPara,'islands')&&isfield(SysPara,'refs')
% %     nIslands=SysPara.nIslands;islands=SysPara.islands;refs=SysPara.refs;
% % else
% %     [nIslands,islands,refs]=searchIslands(bus(:,1),line(:,1:2));
% % end
% % 
% % nbus=size(bus,1);
% % 
% % busType=zeros(nbus,1);
% % if isempty(pv)
% %     pv=zeros(0,6);
% % end
% % if isempty(pq)
% %     pq=zeros(0,6);
% % end
% % if isempty(shunt)
% %     shunt=zeros(0,7);
% % end
% % if isempty(sw)
% %     sw=zeros(0,13);
% % end
% % busType(pv(:,1))=1;
% % busType(sw(:,1))=2;
% % 
% % isw=find(busType==2);
% % ipv=find(busType~=0);
% % ipq=find(busType==0);
% % npq=size(ipq,1);
% % npv=size(ipv,1);
% % 
% % % Determine the frequency model of each island
% % freqTypeTag=zeros(nIslands,1);%0:sw,1:syn,2:steady-state f
% % freqKeptTag=zeros(nbus,1);% corresponds to freqType=2
% % frefs=refs;
% % fswTag=zeros(nbus,1);
% % fsynTag=zeros(nbus,1);
% % fswTag(isw)=1;
% % fswTagxD=fswTag;
% % fsynTag(syn(:,1))=1;
% % D0=imag(xt(vIdx));
% % for isl=1:nIslands
% %     if isempty(find(fswTag(islands==isl)==1, 1))
% %         if isempty(find(fsynTag(islands==isl)==1, 1))
% %             freqTypeTag(isl)=2;
% %             busesInIsland=find(islands==isl);
% %             [~,imin]=min(abs(D0(busesInIsland)));
% %             frefs(isl)=busesInIsland(imin(1));
% %             fswTagxD(frefs(isl))=1;
% %             freqKeptTag(busesInIsland)=1;
% %         else
% %             freqTypeTag(isl)=1;
% %         end
% %     end
% % end
% % freqKeptTagxRef=freqKeptTag;
% % freqKeptTagxRef(frefs)=0;
% % nFreqKept=sum(freqKeptTag);
% % 
% % if ~isempty(agc)
% %     agcExt=zeros(nbus,size(agc,2));
% %     agcExt(agc(:,1),:)=agc;
% %     fdk=agcExt(:,2)+agcExt(:,3); %1/R+D
% % else
% %     fdk=zeros(nbus,1);
% % end
% % fIngIdx=((2*nbus+1):(3*nbus))';
% % 
% % stateFilterTag=ones(nState,1);
% % stateFilterTag(vIdx)=0;
% % stateFilterTag(qIdx)=0;
% % stateFilterTag(pgIdx)=0;
% % stateFilterTag(efIdx)=0;
% % stateFilterTag(tgovmIdx)=0;
% % stateFilterTag(fIdx)=0;
% % % stateFilterTag(dpgIdx)=0;
% % % busFilterTag=ones(nbus,1);
% % % busFilterTag(sw(:,1))=0;
% % % jacBusFilterTag=[busFilterTag;busFilterTag];
% % 
% % busTag=ones(nbus,1);
% % busTag(sw(:,1))=0;
% % jggRowTag=[busTag;busTag;freqKeptTagxRef];
% % jggColTag=[busTag;fswTagxD==0;freqKeptTag];
% % 
% % idxNonSw=find(busType~=2);
% % idxNonSwD=find(busType~=2&fswTagxD==1);
% % 
% % nInd=size(ind,1);
% % if ~isempty(ind)
% %     indIdx=ind(:,1);
% %     Rm1=ind(:,7);
% %     Xm1=ind(:,8);
% %     Zm1=ind(:,7)+1j*ind(:,8);
% %     Rm2=ind(:,9);
% %     Xm2=ind(:,10);
% %     Zmm=1j*ind(:,13);
% %     Tm0=ind(:,15)+ind(:,16)+ind(:,17);
% %     Tm1=-ind(:,16)-2*ind(:,17);
% %     Tm2=ind(:,17);
% %     Hm=ind(:,14);
% % end
% % 
% % nSyn=size(syn,1);
% % if ~isempty(syn)
% %     synIdx=syn(:,1);
% %     wgb=syn(:,4);
% %     model=syn(:,5);
% %     Xgl=syn(:,6);
% %     Rga=syn(:,7);
% %     Xgd=syn(:,8);
% %     Xgd1=syn(:,9);
% %     Xgd2=syn(:,10);
% %     Tgd1=syn(:,11);
% %     Tgd2=syn(:,12);
% %     Xgq=syn(:,13);
% %     Xgq1=syn(:,14);
% %     Xgq2=syn(:,15);
% %     Tgq1=syn(:,16);
% %     Tgq2=syn(:,17);
% %     Mg=syn(:,18);
% %     Dg=syn(:,19);
% %     TgAA=syn(:,24);
% %     gammad=Tgd2./Tgd1.*Xgd2./Xgd1.*(Xgd-Xgd1);
% %     gammaq=Tgq2./Tgq1.*Xgq2./Xgq1.*(Xgq-Xgq1);
% % end
% % 
% % nExc=size(exc,1);
% % if ~isempty(exc)
% %     excIdx=exc(:,1);
% %     VavrMax=exc(:,3);
% %     VavrMin=exc(:,4);
% %     muavr0=exc(:,5);
% %     Tavr1=exc(:,7);
% %     Tavr2=exc(:,6);
% %     vavrf0=exc(:,8);
% %     Vavr0=exc(:,9);
% %     Tavre=exc(:,10);
% %     Tavrr=exc(:,11);
% % end
% % 
% % nTg=size(tg,1);
% % if ~isempty(tg)
% %     tgIdx=tg(:,1);
% %     wtgref=tg(:,3);
% %     Rtg=tg(:,4);
% %     Ttgmax=tg(:,5);
% %     Ttgmin=tg(:,6);
% %     Ttg2=tg(:,7);
% %     Ttg1=tg(:,8);
% % end
% % 
% % [V0,Q0,s0,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,tgovmech0,f0,dpg0,qplt0,vg0]=unfoldX(xt,SysData);
% % if ~isempty(ind)
% %     [YshInd0,Yshind1]=getLinearInterpolatorInd(nbus,ind,s0,s0);
% %     Ysh0x=Ysh0+YshInd0;
% %     Ysh1x=Ysh1+Yshind1;
% % else
% %     Ysh0x=zeros(nbus,1);
% %     Ysh1x=zeros(nbus,1);
% % end
% % Ytr=Ytr0+Ytr1*t;
% % Ysh=Ysh0x+Ysh1x*t;
% % YtrNew=Ytr0+Ytr1*(t+dt);
% % YshNew=Ysh0x+Ysh1x*(t+dt);
% % 
% % SNew=-accumarray(pq(:,1),pq(:,4)+1j*pq(:,5),[nbus,1])+accumarray(pv(:,1),pv(:,4),[nbus,1]);
% % SNew=SNew-accumarray(pq(:,1),pqIncr(:,1)+1j*pqIncr(:,2),[nbus,1])*(t+dt)+accumarray(pv(:,1),pvIncr,[nbus,1])*(t+dt);
% % if size(pqIncr,2)>=4;SNew=SNew-accumarray(pq(:,1),pqIncr(:,3)+1j*pqIncr(:,4),[nbus,1])*(t+dt)*(t+dt);end
% % SNew=SNew-accumarray(zip(:,1),(Rzip0+Rzip1*(t+dt)).*(zip(:,7)+1j*zip(:,10)),[nbus,1]);
% % 
% % INew=-accumarray(zip(:,1),(Rzip0+Rzip1*(t+dt)).*(zip(:,6)-1j*zip(:,9)),[nbus,1]);
% % JNew=real(INew);
% % KNew=imag(INew);
% % 
% % YNew=YtrNew+sparse(1:nbus,1:nbus,YshNew,nbus,nbus);
% % GNew=real(YNew);BNew=imag(YNew);
% % 
% % % Calculate diff and jac, do N-R
% % Vtemp=V0;
% % pqx=pq;pqx(:,[4,5])=pq(:,[4,5])+pqIncr(:,1:2)*t;if size(pqIncr,2)>=4;pqx(:,[4,5])=pqx(:,[4,5])+pqIncr(:,3:4)*t*t;end
% % pvx=pv;pvx(:,4)=pv(:,4)+pvIncr*t;
% % indx=ind;if ~isempty(ind);indx(:,15:17)=ind(:,15:17).*repmat(Rind0+Rind1*t,1,3);indx(:,13)=ind(:,13)./(Reind0+Reind1*t);end
% % zipx=zip;if ~isempty(zip);zipx(:,5:10)=zip(:,5:10).*repmat(Rzip0+Rzip1*t,1,6);end
% % SysDatax=foldSysData(bus,sw,pvx,pqx,shunt,line,indx,zipx,syn,exc,tg,agc,cac,cluster);
% % 
% % pqxx=pq;pqxx(:,[4,5])=pq(:,[4,5])+pqIncr(:,1:2)*(t+dt);if size(pqIncr,2)>=4;pqxx(:,[4,5])=pqxx(:,[4,5])+pqIncr(:,3:4)*(t+dt)*(t+dt);end
% % pvxx=pv;pvxx(:,4)=pv(:,4)+pvIncr*(t+dt);
% % indxx=ind;if ~isempty(ind);indxx(:,15:17)=ind(:,15:17).*repmat(Rind0+Rind1*(t+dt),1,3);indxx(:,13)=ind(:,13)./(Reind0+Reind1*(t+dt));end
% % zipxx=zip;if ~isempty(zip);zipxx(:,5:10)=zip(:,5:10).*repmat(Rzip0+Rzip1*(t+dt),1,6);end
% % SysDataxx=foldSysData(bus,sw,pvxx,pqxx,shunt,line,indxx,zipxx,syn,exc,tg,agc,cac,cluster);
% % 
% % SysData=foldSysData(bus,sw,pvx,pqx,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster);
% % dx1=integ(SysDatax,SysPara,xt,dt);
% % xTemp=xt;
% % xTemp=xTemp+dx1;
% % 
% % iterMax=10;
% % exitflag=-1;
% % for iter=1:iterMax
% %     % calculate df and dg
% %     dx2=integ(SysDataxx,SysPara,xTemp,dt);
% %     diffF=xt-xTemp+(dx1+dx2)/2;
% %     diffF(stateFilterTag==0)=0;
% %     
% %     Ctemp=real(Vtemp);
% %     Dtemp=imag(Vtemp);
% %     Atemp=abs(Vtemp);
% %     Vgtemp=Vtemp(syn(:,1));
% %     ftemp=xTemp(fIdx);
% %     
% %     eq0=Ef0;ed0=zeros(size(Ef0));
% %     [~,~,sNew,dNew,wNew,eq1New,eq2New,ed1New,ed2New,psidNew,psiqNew,EdNew,EfNew,VavrmNew,VavrrNew,VavrfNew,VavrrefNew,tgovgNew,tgovmNew,tgovmechNew,fNew,dpgNew,qpltNew,vgNew]=unfoldX(xTemp,SysData);
% %     eqNew=EfNew;edNew=zeros(size(EfNew));
% %     
% %     [MatGV0,MatGV1,MatGRhs0,MatGRhs1]=getLinearInterpolatorSyn(syn,[],d0,dNew,ed0,ed10,ed20,edNew,ed1New,ed2New,eq0,eq10,eq20,eqNew,eq1New,eq2New,psid0,psiq0,psidNew,psiqNew);
% %     
% %     MatGVNew=MatGV0+MatGV1;
% %     MatGRhsNew=MatGRhs0+MatGRhs1;
% %     J1=GNew+sparse(syn(:,1),syn(:,1),MatGVNew(:,1),nbus,nbus);
% %     J2=-BNew+sparse(syn(:,1),syn(:,1),MatGVNew(:,2),nbus,nbus);
% %     J3=BNew+sparse(syn(:,1),syn(:,1),MatGVNew(:,3),nbus,nbus);
% %     J4=GNew+sparse(syn(:,1),syn(:,1),MatGVNew(:,4),nbus,nbus);
% %     P2freq=sparse(1:nbus,1:nbus,-freqKeptTag.*fdk,nbus,nbus);
% %     Freq2freq=sparse([1:nbus,1:nbus],[1:nbus,frefs(islands)'],[ones(1,nbus),-ones(1,nbus)],nbus,nbus);
% %     
% %     Ig=accumarray(syn(:,1),MatGRhsNew(:,1)+1j*MatGRhsNew(:,2),[nbus,1])-...
% %         accumarray(syn(:,1),(MatGVNew(:,1).*real(Vgtemp)+MatGVNew(:,2).*imag(Vgtemp))+1j*(MatGVNew(:,3).*real(Vgtemp)+MatGVNew(:,4).*imag(Vgtemp)),[nbus,1]);
% % %     diffG=conj(SNew)+INew.*abs(Vtemp)+Ig.*conj(Vtemp)-(YNew*Vtemp).*conj(Vtemp)+freqKeptTag.*(dpgNew-fdk.*fNew);
% % %     diffG(busFilterTag==0)=0;
% % %     diffGf=zeros(nbus,1);
% %     
% %     diff=conj(SNew)+INew.*abs(Vtemp)+Ig.*conj(Vtemp)-(YNew*Vtemp).*conj(Vtemp);
% %     diffP=real(diff)+freqKeptTag.*(-fdk.*ftemp+dpg0);
% %     diffVQ=imag(diff);
% %     diffGf=ftemp-ftemp(frefs(islands));
% %     if ~isempty(pv)
% %         diffVQ(pv(:,1))=abs(Vtemp(pv(:,1))).*abs(Vtemp(pv(:,1)))-(VspSq2(pv(:,1),1)+VspSq2(pv(:,1),2)*(t+dt));
% %     end
% %     
% %     diffGx=[diffP;diffVQ;diffGf];
% %     diffGx(jggRowTag==0)=0;
% %     
% %     if iter>1&&max(abs([diffF;diffGx]))<diffTol
% %         exitflag=0;
% %         break;
% %     end
% %     
% %     % construct Jgg
% %     dCtemp=sparse(1:nbus,1:nbus,Ctemp,nbus,nbus);
% %     dDtemp=sparse(1:nbus,1:nbus,Dtemp,nbus,nbus);
% %     Jac1=[dCtemp*J1+dDtemp*J3,dCtemp*J2+dDtemp*J4;...
% %         dCtemp*J3-dDtemp*J1,dCtemp*J4-dDtemp*J2];
% %     Jac2=[sparse(1:nbus,1:nbus,J1*Ctemp+J2*Dtemp,nbus,nbus),sparse(1:nbus,1:nbus,J3*Ctemp+J4*Dtemp,nbus,nbus);...
% %         sparse(1:nbus,1:nbus,J3*Ctemp+J4*Dtemp,nbus,nbus),-sparse(1:nbus,1:nbus,J1*Ctemp+J2*Dtemp,nbus,nbus)];
% %     Jac3=[sparse(1:nbus,1:nbus,JNew.*Ctemp./Atemp,nbus,nbus),sparse(1:nbus,1:nbus,JNew.*Dtemp./Atemp,nbus,nbus);...
% %         sparse(1:nbus,1:nbus,KNew.*Ctemp./Atemp,nbus,nbus),sparse(1:nbus,1:nbus,KNew.*Dtemp./Atemp,nbus,nbus)];
% %     Jac4=[sparse(syn(:,1),syn(:,1),MatGRhsNew(:,1),nbus,nbus),sparse(syn(:,1),syn(:,1),MatGRhsNew(:,2),nbus,nbus);...
% %         sparse(syn(:,1),syn(:,1),MatGRhsNew(:,2),nbus,nbus),-sparse(syn(:,1),syn(:,1),MatGRhsNew(:,1),nbus,nbus)];
% %     Jac=Jac3+Jac4-Jac1-Jac2;
% %     
% %     JacVc=sparse(1:nbus,1:nbus,2*Ctemp,nbus,nbus);
% %     JacVd=sparse(1:nbus,1:nbus,2*Dtemp,nbus,nbus);
% %     
% %     Jac(nbus+ipv,:)=[JacVc(ipv,:),JacVd(ipv,:)];
% %     
% %     Jgg=[Jac,[P2freq;sparse(nbus,nbus)];...
% %         sparse(nbus,2*nbus),Freq2freq];
% %     
% %     % construct Jgf
% %     
% %     jgfi=zeros(0,1);
% %     jgfj=zeros(0,1);
% %     jgfv=zeros(0,1);
% %     jfgi=zeros(0,1);
% %     jfgj=zeros(0,1);
% %     jfgv=zeros(0,1);
% %     jffi=zeros(0,1);
% %     jffj=zeros(0,1);
% %     jffv=zeros(0,1);
% %     
% %     if ~isempty(ind)
% %         % s
% %         deltasx=1e-3;
% %         Zind0=Zm1+Zmm.*(Rm2+1j*Xm2.*sNew)./(Rm2+(Zmm+1j*Xm2).*sNew);
% %         Yind0=1./Zind0;
% %         
% %         ZindNew=Zm1+Zmm.*(Rm2+1j*Xm2.*(sNew+deltasx))./(Rm2+(Zmm+1j*Xm2).*(sNew+deltasx));
% %         YindNew=1./ZindNew;
% %         
% %         Yshind1x=YindNew-Yind0;
% %         jgfi=[jgfi;indIdx;indIdx+nbus];
% %         jgfj=[jgfj;sIdx;sIdx];
% %         jgfv=[jgfv;Vtemp(indIdx).*conj(Vtemp(indIdx)).*real(Yshind1x/deltasx);Vtemp(indIdx).*conj(Vtemp(indIdx)).*imag(Yshind1x/deltasx)];
% %     end
% %     
% %     if ~isempty(syn)
% %         Vgtemp=Vtemp(synIdx);
% %         Cgtemp=real(Vgtemp);
% %         Dgtemp=imag(Vgtemp);
% %         
% %         cosd0=cos(dNew);
% %         sind0=sin(dNew);
% %         deltadNew=1e-5*abs(dNew)+1e-6;
% %         cosdx=cos(dNew+deltadNew);
% %         sindx=sin(dNew+deltadNew);
% %         
% %         MatGV0x=zeros(nSyn,4);
% %         MatGV1x=zeros(nSyn,4);
% %         
% %         MatGRhs0x=zeros(nSyn,2);
% %         MatGRhs1x=zeros(nSyn,2);
% %         
% %         for i=1:nSyn
% %             Mg0=[sind0(i),-cosd0(i);cosd0(i),sind0(i)];
% %             Mgx=[sindx(i),-cosdx(i);cosdx(i),sindx(i)];
% %             if model(i)==6||model(i)==5
% %                 Yg=[Rga(i),Xgq2(i);-Xgd2(i),Rga(i)]/(Rga(i)*Rga(i)+Xgd2(i)*Xgq2(i));
% %                 Mxye=Mg0'*Yg;
% %                 Mdqe=Yg;
% %                 dIdq=-Yg*(Mgx-Mg0)*[Cgtemp(i);Dgtemp(i)]/deltadNew(i);
% %                 dVdq=-Yg*Mg0;
% %                 jgfi=[jgfi;synIdx(i);synIdx(i);synIdx(i)+nbus;synIdx(i)+nbus];
% %                 jgfj=[jgfj;ed2Idx(i);eq2Idx(i);ed2Idx(i);eq2Idx(i)];
% %                 jgfv=[jgfv;Cgtemp(i)*Mxye(1,1)+Dgtemp(i)*Mxye(2,1);Cgtemp(i)*Mxye(1,2)+Dgtemp(i)*Mxye(2,2);Cgtemp(i)*Mxye(2,1)-Dgtemp(i)*Mxye(1,1);Cgtemp(i)*Mxye(2,2)-Dgtemp(i)*Mxye(1,2)];
% %                 
% %                 jffi=[jffi;eq2Idx(i);eq2Idx(i);eq2Idx(i)];
% %                 jffj=[jffj;ed2Idx(i);eq2Idx(i);deltaIdx(i)];
% %                 jffv=[jffv;-(Xgd1(i)-Xgd2(i)+gammad(i))*Mdqe(1,1)/Tgd2(i);-(Xgd1(i)-Xgd2(i)+gammad(i))*Mdqe(1,2)/Tgd2(i);-(Xgd1(i)-Xgd2(i)+gammad(i))*dIdq(1)/Tgd2(i)];
% %                 jffi=[jffi;ed2Idx(i);ed2Idx(i);ed2Idx(i)];
% %                 jffj=[jffj;ed2Idx(i);eq2Idx(i);deltaIdx(i)];
% %                 jffv=[jffv;(Xgq1(i)-Xgq2(i)+gammaq(i))*Mdqe(2,1)/Tgq2(i);(Xgq1(i)-Xgq2(i)+gammaq(i))*Mdqe(2,2)/Tgq2(i);(Xgq1(i)-Xgq2(i)+gammaq(i))*dIdq(2)/Tgq2(i)];
% %                 jffi=[jffi;eq1Idx(i);eq1Idx(i);eq1Idx(i)];
% %                 jffj=[jffj;ed2Idx(i);eq2Idx(i);deltaIdx(i)];
% %                 jffv=[jffv;-(Xgd(i)-Xgd1(i)-gammad(i))*Mdqe(1,1)/Tgd1(i);-(Xgd(i)-Xgd1(i)-gammad(i))*Mdqe(1,2)/Tgd1(i);-(Xgd(i)-Xgd1(i)-gammad(i))*dIdq(1)/Tgd1(i)];
% %                 if model(i)==6
% %                     jffi=[jffi;ed1Idx(i);ed1Idx(i);ed1Idx(i)];
% %                     jffj=[jffj;ed2Idx(i);eq2Idx(i);deltaIdx(i)];
% %                     jffv=[jffv;(Xgq(i)-Xgq1(i)-gammaq(i))*Mdqe(2,1)/Tgq1(i);(Xgq(i)-Xgq1(i)-gammaq(i))*Mdqe(2,2)/Tgq1(i);(Xgq(i)-Xgq1(i)-gammaq(i))*dIdq(2)/Tgq1(i)];
% %                 end
% %                 
% %                 jfgi=[jfgi;eq2Idx(i);eq2Idx(i)];
% %                 jfgj=[jfgj;synIdx(i);synIdx(i)+nbus];
% %                 jfgv=[jfgv;-(Xgd1(i)-Xgd2(i)+gammad(i))*dVdq(1,1)/Tgd2(i);-(Xgd1(i)-Xgd2(i)+gammad(i))*dVdq(1,2)/Tgd2(i)];
% %                 jfgi=[jfgi;ed2Idx(i);ed2Idx(i)];
% %                 jfgj=[jfgj;synIdx(i);synIdx(i)+nbus];
% %                 jfgv=[jfgv;(Xgq1(i)-Xgq2(i)+gammaq(i))*dVdq(2,1)/Tgq2(i);(Xgq1(i)-Xgq2(i)+gammaq(i))*dVdq(2,2)/Tgq2(i)];
% %                 jfgi=[jfgi;eq1Idx(i);eq1Idx(i)];
% %                 jfgj=[jfgj;synIdx(i);synIdx(i)+nbus];
% %                 jfgv=[jfgv;-(Xgd(i)-Xgd1(i)-gammad(i))*dVdq(1,1)/Tgd1(i);-(Xgd(i)-Xgd1(i)-gammad(i))*dVdq(1,2)/Tgd1(i)];
% %                 if model(i)==6
% %                     jfgi=[jfgi;ed1Idx(i);ed1Idx(i)];
% %                     jfgj=[jfgj;synIdx(i);synIdx(i)+nbus];
% %                     jfgv=[jfgv;(Xgq(i)-Xgq1(i)-gammaq(i))*Mdqe(2,1)/Tgq1(i);(Xgq(i)-Xgq1(i)-gammaq(i))*Mdqe(2,2)/Tgq1(i)];
% %                 end
% %                 
% %                 MatGRhs0x(i,:)=(Mg0'*Yg*[ed2New(i);eq2New(i)]).';
% %                 MatGRhs1x(i,:)=(Mgx'*Yg*[ed2New(i);eq2New(i)]).'-MatGRhs0x(i,:);
% %             elseif model(i)==4
% %                 Yg=[Rga(i),Xgq1(i);-Xgd1(i),Rga(i)]/(Rga(i)*Rga(i)+Xgd1(i)*Xgq1(i));
% %                 Mxye=Mg0'*Yg;
% %                 Mdqe=Yg;
% %                 dVdq=-Yg*Mg0;
% %                 dIdq=-Yg*(Mgx-Mg0)*[Cgtemp(i);Dgtemp(i)]/deltadNew(i);
% %                 jgfi=[jgfi;synIdx(i);synIdx(i);synIdx(i)+nbus;synIdx(i)+nbus];
% %                 jgfj=[jgfj;ed1Idx(i);eq1Idx(i);ed1Idx(i);eq1Idx(i)];
% %                 jgfv=[jgfv;Cgtemp(i)*Mxye(1,1)+Dgtemp(i)*Mxye(2,1);Cgtemp(i)*Mxye(1,2)+Dgtemp(i)*Mxye(2,2);Cgtemp(i)*Mxye(2,1)-Dgtemp(i)*Mxye(1,1);Cgtemp(i)*Mxye(2,2)-Dgtemp(i)*Mxye(1,2)];
% %                 
% %                 jffi=[jffi;eq1Idx(i);eq1Idx(i);eq1Idx(i)];
% %                 jffj=[jffj;ed1Idx(i);eq1Idx(i);deltaIdx(i)];
% %                 jffv=[jffv;-(Xgd(i)-Xgd1(i))*Mdqe(1,1)/Tgd1(i);-(Xgd(i)-Xgd1(i))*Mdqe(1,2)/Tgd1(i);-(Xgd(i)-Xgd1(i))*dIdq(1)/Tgd1(i)];
% %                 jffi=[jffi;ed1Idx(i);ed1Idx(i);ed1Idx(i)];
% %                 jffj=[jffj;ed1Idx(i);eq1Idx(i);deltaIdx(i)];
% %                 jffv=[jffv;(Xgq(i)-Xgq1(i))*Mdqe(2,1)/Tgq1(i);(Xgq(i)-Xgq1(i))*Mdqe(2,2)/Tgq1(i);(Xgq(i)-Xgq1(i))*dIdq(2)/Tgq1(i)];
% %                 
% %                 jfgi=[jfgi;eq1Idx(i);eq1Idx(i)];
% %                 jfgj=[jfgj;synIdx(i);synIdx(i)+nbus];
% %                 jfgv=[jfgv;-(Xgd(i)-Xgd1(i))*dVdq(1,1)/Tgd1(i);-(Xgd(i)-Xgd1(i))*dVdq(1,2)/Tgd1(i)];
% %                 jfgi=[jfgi;ed1Idx(i);ed1Idx(i)];
% %                 jfgj=[jfgj;synIdx(i);synIdx(i)+nbus];
% %                 jfgv=[jfgv;(Xgq(i)-Xgq1(i))*Mdqe(2,1)/Tgq1(i);(Xgq(i)-Xgq1(i))*Mdqe(2,2)/Tgq1(i)];
% %                 
% %                 MatGRhs0x(i,:)=(Mg0'*Yg*[ed1New(i);eq1New(i)]).';
% %                 MatGRhs1x(i,:)=(Mgx'*Yg*[ed1New(i);eq1New(i)]).'-MatGRhs0x(i,:);
% %             elseif model(i)==3||model(i)==2
% %                 Yg=[Rga(i),Xgq(i);-Xgd1(i),Rga(i)]/(Rga(i)*Rga(i)+Xgd1(i)*Xgq(i));
% %                 Mxye=Mg0'*Yg;
% %                 Mdqe=Yg;
% %                 dVdq=-Yg*Mg0;
% %                 dIdq=-Yg*(Mgx-Mg0)*[Cgtemp(i);Dgtemp(i)]/deltadNew(i);
% %                 if model(i)==3
% %                     jgfi=[jgfi;synIdx(i);synIdx(i)+nbus];
% %                     jgfj=[jgfj;eq1Idx(i);eq1Idx(i)];
% %                     jgfv=[jgfv;Cgtemp(i)*Mxye(1,2)+Dgtemp(i)*Mxye(2,2);Cgtemp(i)*Mxye(2,2)-Dgtemp(i)*Mxye(1,2)];
% %                     jffi=[jffi;eq1Idx(i);eq1Idx(i)];
% %                     jffj=[jffj;eq1Idx(i);deltaIdx(i)];
% %                     jffv=[jffv;-(Xgd(i)-Xgd1(i))*Mdqe(1,1)/Tgd1(i);-(Xgd(i)-Xgd1(i))*dIdq(1)/Tgd1(i)];
% %                     jfgi=[jfgi;eq1Idx(i);eq1Idx(i)];
% %                     jfgj=[jfgj;synIdx(i);synIdx(i)+nbus];
% %                     jfgv=[jfgv;-(Xgd(i)-Xgd1(i))*dVdq(1,1)/Tgd1(i);-(Xgd(i)-Xgd1(i))*dVdq(1,2)/Tgd1(i)];
% %                 end
% %                 
% %                 MatGRhs0x(i,:)=(Mg0'*Yg*[edNew(i);eq1New(i)]).';
% %                 MatGRhs1x(i,:)=(Mgx'*Yg*[edNew(i);eq1New(i)]).'-MatGRhs0x(i,:);
% %             elseif model(i)==8
% %                 Yg=zeros(2);
% %                 Mxye=Mg0'*diag([1/Xgd2(i);1/Xgq2(i)]);
% %                 Mdqe=diag([1/Xgd2(i);1/Xgq2(i)]);
% %                 dVdq=-Yg*Mg0;
% %                 dIdq=-Yg*(Mgx-Mg0)*[Cgtemp(i);Dgtemp(i)]/deltadNew(i);
% %                 jgfi=[jgfi;synIdx(i);synIdx(i);synIdx(i)+nbus;synIdx(i)+nbus];
% %                 jgfj=[jgfj;ed2Idx(i);eq2Idx(i);ed2Idx(i);eq2Idx(i)];
% %                 jgfv=[jgfv;-(Cgtemp(i)*Mxye(1,1)+Dgtemp(i)*Mxye(2,1));Cgtemp(i)*Mxye(1,2)+Dgtemp(i)*Mxye(2,2);-(Cgtemp(i)*Mxye(2,1)-Dgtemp(i)*Mxye(1,1));Cgtemp(i)*Mxye(2,2)-Dgtemp(i)*Mxye(1,2)];
% %                 jgfi=[jgfi;synIdx(i);synIdx(i);synIdx(i)+nbus;synIdx(i)+nbus];
% %                 jgfj=[jgfj;psiqIdx(i);psidIdx(i);psiqIdx(i);psidIdx(i)];
% %                 jgfv=[jgfv;-(Cgtemp(i)*Mxye(1,1)+Dgtemp(i)*Mxye(2,1));-(Cgtemp(i)*Mxye(1,2)+Dgtemp(i)*Mxye(2,2));-(Cgtemp(i)*Mxye(2,1)-Dgtemp(i)*Mxye(1,1));-(Cgtemp(i)*Mxye(2,2)-Dgtemp(i)*Mxye(1,2))];
% %                 
% %                 jffi=[jffi;eq2Idx(i);eq2Idx(i);eq2Idx(i);eq2Idx(i);eq2Idx(i)];
% %                 jffj=[jffj;ed2Idx(i);eq2Idx(i);psiqIdx(i);psidIdx(i);deltaIdx(i)];
% %                 jffv=[jffv;-(Xgd1(i)-Xgd2(i)+gammad(i))*Mdqe(1,1)/Tgd2(i);-(Xgd1(i)-Xgd2(i)+gammad(i))*Mdqe(1,2)/Tgd2(i);...
% %                     -(Xgd1(i)-Xgd2(i)+gammad(i))*Mdqe(1,1)/Tgd2(i);(Xgd1(i)-Xgd2(i)+gammad(i))*Mdqe(1,2)/Tgd2(i);-(Xgd1(i)-Xgd2(i)+gammad(i))*dIdq(1)/Tgd2(i)];
% %                 jffi=[jffi;ed2Idx(i);ed2Idx(i);ed2Idx(i);ed2Idx(i);ed2Idx(i)];
% %                 jffj=[jffj;ed2Idx(i);eq2Idx(i);psiqIdx(i);psidIdx(i);deltaIdx(i)];
% %                 jffv=[jffv;(Xgq1(i)-Xgq2(i)+gammaq(i))*Mdqe(2,1)/Tgq2(i);(Xgq1(i)-Xgq2(i)+gammaq(i))*Mdqe(2,2)/Tgq2(i);...
% %                     (Xgq1(i)-Xgq2(i)+gammaq(i))*Mdqe(2,1)/Tgq2(i);-(Xgq1(i)-Xgq2(i)+gammaq(i))*Mdqe(2,2)/Tgq2(i);(Xgq1(i)-Xgq2(i)+gammaq(i))*dIdq(2)/Tgq2(i)];
% %                 jffi=[jffi;eq1Idx(i);eq1Idx(i);eq1Idx(i);eq1Idx(i);eq1Idx(i)];
% %                 jffj=[jffj;ed2Idx(i);eq2Idx(i);psiqIdx(i);psidIdx(i);deltaIdx(i)];
% %                 jffv=[jffv;-(Xgd(i)-Xgd1(i)-gammad(i))*Mdqe(1,1)/Tgd1(i);-(Xgd(i)-Xgd1(i)-gammad(i))*Mdqe(1,2)/Tgd1(i);...
% %                     -(Xgd(i)-Xgd1(i)-gammad(i))*Mdqe(1,1)/Tgd1(i);(Xgd(i)-Xgd1(i)-gammad(i))*Mdqe(1,2)/Tgd1(i);-(Xgd(i)-Xgd1(i)-gammad(i))*dIdq(1)/Tgd1(i)];
% %                 jffi=[jffi;ed1Idx(i);ed1Idx(i);ed1Idx(i);ed1Idx(i);ed1Idx(i)];
% %                 jffj=[jffj;ed2Idx(i);eq2Idx(i);psiqIdx(i);psidIdx(i);deltaIdx(i)];
% %                 jffv=[jffv;(Xgq(i)-Xgq1(i)-gammaq(i))*Mdqe(2,1)/Tgq1(i);(Xgq(i)-Xgq1(i)-gammaq(i))*Mdqe(2,2)/Tgq1(i);...
% %                     (Xgq(i)-Xgq1(i)-gammaq(i))*Mdqe(2,1)/Tgq1(i);-(Xgq(i)-Xgq1(i)-gammaq(i))*Mdqe(2,2)/Tgq1(i);(Xgq(i)-Xgq1(i)-gammaq(i))*dIdq(2)/Tgq1(i)];
% %                 
% %                 MatGRhs0x(i,:)=(Mg0'*[(eq2New(i)-psidNew(i))/Xgd2(i);(-ed2New(i)-psiqNew(i))/Xgq2(i)]).';
% %                 MatGRhs1x(i,:)=(Mgx'*[(eq2New(i)-psidNew(i))/Xgd2(i);(-ed2New(i)-psiqNew(i))/Xgq2(i)]).'-MatGRhs0x(i,:);
% %             end
% %             MatVg0x=Mg0'*Yg*Mg0;
% %             MatVgNewx=Mgx'*Yg*Mgx;
% %             MatGV0x(i,:)=[MatVg0x(1,1),MatVg0x(1,2),MatVg0x(2,1),MatVg0x(2,2)];
% %             MatGV1x(i,:)=[MatVgNewx(1,1),MatVgNewx(1,2),MatVgNewx(2,1),MatVgNewx(2,2)]-MatGV0x(i,:);
% %             
% %             jgfi=[jgfi;synIdx(i);synIdx(i)+nbus];
% %             jgfj=[jgfj;deltaIdx(i);deltaIdx(i)];
% %             jgfv=[jgfv;...
% %                 (Cgtemp(i)*(MatGRhs1x(i,1)-MatGV1x(i,1)*Cgtemp(i)-MatGV1x(i,2)*Dgtemp(i))+Dgtemp(i)*(MatGRhs1x(i,2)-MatGV1x(i,3)*Cgtemp(i)-MatGV1x(i,4)*Dgtemp(i)))/deltadNew(i);...
% %                 (-Dgtemp(i)*(MatGRhs1x(i,1)-MatGV1x(i,1)*Cgtemp(i)-MatGV1x(i,2)*Dgtemp(i))+Cgtemp(i)*(MatGRhs1x(i,2)-MatGV1x(i,3)*Cgtemp(i)-MatGV1x(i,4)*Dgtemp(i)))/deltadNew(i)];
% %         end
% %         
% %     end
% %     Jgf=sparse(jgfi,jgfj,jgfv,3*nbus,nState);
% %     
% %     % construct Jff & Jfg
% %     
% %     if ~isempty(ind)
% %         % s
% %         jffi=[jffi;sIdx];
% %         jffj=[jffj;sIdx];
% %         divs=((Rm2+Rm1.*sNew).*(Rm2+Rm1.*sNew)+(Xm2+Xm1).*(Xm2+Xm1).*sNew.*sNew);
% %         jffv=[jffv;(Tm1+2*Tm2.*sNew+Vtemp(indIdx).*conj(Vtemp(indIdx)).*Rm2.*(1./divs-2*sNew.*(sNew.*(Xm2+Xm1).*(Xm2+Xm1)+Rm1.*(Rm2+Rm1.*sNew))./divs./divs))/2./Hm];
% %         %             il=V(indIdx).*yind;
% %         %             ir=il-(V(indIdx)-il.*Zm1)./Zmm;
% %         %             dx(sIdx)=dt*((Tm0+Tm1.*s+Tm2.*s.*s)-ir.*conj(ir).*Rm2./s)/2./Hm;
% %         
% %         yind=1./(Zm1+Zmm.*(Rm2+1j*Xm2.*sNew)./(Rm2+sNew.*(Zmm+1j*Xm2)));
% %         yInds=(1+Zm1./Zmm).*yind-1./Zmm;
% %         yIndsRs=yInds.*conj(yInds).*Rm2./sNew;
% %         yIndsRs(isnan(yIndsRs))=0;
% %         jfgi=[jfgi;sIdx;sIdx];
% %         jfgj=[jfgj;indIdx;indIdx+nbus];
% %         jfgv=[jfgv;2*real(Vtemp(indIdx)).*yIndsRs/2./Hm;2*imag(Vtemp(indIdx)).*yIndsRs/2./Hm];
% %     end
% %     
% %     if ~isempty(syn)
% %         pActiveInSyn=ones(nSyn,1);        
% %         tgovmActive=ones(nTg,1);
% %         tgovmActive((tgovm0-Ttgmax)>0)=0;
% %         tgovmActive((tgovm0-Ttgmin)<0)=0;
% %         pActiveInSyn(tgIdx)=tgovmActive;
% %         texcActive=ones(nExc,1);
% %         texcActive((Vavrf0-VavrMax)>0)=0;
% %         texcActive((Vavrf0-VavrMin)<0)=0;
% %         texcActiveInSyn=zeros(nSyn,1);
% %         texcActiveInSyn(tgIdx)=texcActive;        
% %         numSynOnBus=accumarray(syn(:,1),1,[nbus,1]);
% %         %delta
% %         jffi=[jffi;deltaIdx];
% %         jffj=[jffj;omegaIdx];
% %         jffv=[jffv;wgb];
% %         %omega
% %         jffi=[jffi;omegaIdx;omegaIdx(tgIdx);omegaIdx(tgIdx);omegaIdx(tgIdx);omegaIdx];
% %         jffj=[jffj;omegaIdx;tgovgIdx;omegaIdx(tgIdx);tmechIdx;dpgIdx(synIdx)];
% %         jffv=[jffv;-Dg/2./Mg;tgovmActive;-tgovmActive.*Ttg1./Ttg2./Rtg;tgovmActive;pActiveInSyn./numSynOnBus(synIdx)];
% %         %psiq
% %         jffi=[jffi;psidIdx(model>=8)];
% %         jffj=[jffj;psiqIdx(model>=8)];
% %         jffv=[jffv;wgb(model>=8)];
% %         %psid
% %         jffi=[jffi;psiqIdx(model>=8)];
% %         jffj=[jffj;psidIdx(model>=8)];
% %         jffv=[jffv;-wgb(model>=8)];
% %         %ed''
% %         jffi=[jffi;ed2Idx(model>=5);ed2Idx(model>=5)];
% %         jffj=[jffj;ed2Idx(model>=5);ed1Idx(model>=5)];
% %         jffv=[jffv;-1./Tgq2(model>=5);1./Tgq2(model>=5)];
% %         %eq''
% %         jffi=[jffi;eq2Idx(model>=5);eq2Idx(model>=5);eq2Idx(model>=5&texcActiveInSyn==1)];
% %         jffj=[jffj;eq2Idx(model>=5);eq1Idx(model>=5);vavrfIdx(texcActive==1&model(excIdx)>=5)];
% %         jffv=[jffv;-1./Tgd2(model>=5);1./Tgd2(model>=5);TgAA(model>=5&texcActiveInSyn==1)./Tgd1(model>=5&texcActiveInSyn==1)./Tgd2(model>=5&texcActiveInSyn==1)];
% %         %ed'
% %         jffi=[jffi;ed1Idx(model>=4)];
% %         jffj=[jffj;ed1Idx(model>=4)];
% %         jffv=[jffv;-1./Tgq1(model>=4)];
% %         %eq'
% %         jffi=[jffi;eq1Idx(model>=3)];
% %         jffj=[jffj;eq1Idx(model>=3)];
% %         jffv=[jffv;-1./Tgd1(model>=3)];
% %         jffi=[jffi;eq1Idx(model>=5&texcActiveInSyn==1)];
% %         jffj=[jffj;vavrfIdx(texcActive==1&model(excIdx)>=5)];
% %         jffv=[jffv;(1-TgAA(model>=5&texcActiveInSyn==1)./Tgd1(model>=5&texcActiveInSyn==1))./Tgd2(model>=5&texcActiveInSyn==1)];
% %         jffi=[jffi;eq1Idx((model==3|model==4)&texcActiveInSyn==1)];
% %         jffj=[jffj;vavrfIdx(texcActive==1&(model(excIdx)==3|model(excIdx)==4))];
% %         jffv=[jffv;(1-TgAA((model==3|model==4)&texcActiveInSyn==1)./Tgd1((model==3|model==4)&texcActiveInSyn==1))./Tgd2((model==3|model==4)&texcActiveInSyn==1)];
% %     end
% %     
% %     if ~isempty(tg)
% %         jffi=[jffi;tgovgIdx;tgovgIdx];
% %         jffj=[jffj;omegaIdx;tgovgIdx];
% %         jffv=[jffv;-(1-Ttg1./Ttg2)./Rtg./Ttg2;-1./Ttg2];
% %     end
% %     
% %     if ~isempty(exc)
% %         VmagAvr=abs(Vtemp(synIdx(excIdx)));
% %         
% %         jffi=[jffi;vavrmIdx];
% %         jffj=[jffj;vavrmIdx];
% %         jffv=[jffv;-1./Tavrr];
% %         
% %         jffi=[jffi;vavrrIdx;vavrrIdx;vavrrIdx];
% %         jffj=[jffj;vavrmIdx;vavrrefIdx;vavrrIdx];
% %         jffv=[jffv;-muavr0.*(1-Tavr1./Tavr2)./Tavr2;muavr0.*(1-Tavr1./Tavr2)./Tavr2;-1./Tavr2];
% %         
% %         jffi=[jffi;vavrfIdx;vavrfIdx;vavrfIdx;vavrfIdx];
% %         jffj=[jffj;vavrmIdx;vavrrefIdx;vavrrIdx;vavrfIdx];
% %         jffv=[jffv;-muavr0.*Tavr1./Tavr2./Tavre.*VmagAvr./Vavr0;muavr0.*Tavr1./Tavr2./Tavre.*VmagAvr./Vavr0;VmagAvr./Vavr0./Tavre;-1./Tavre];
% %     end
% %     
% %     if ~isempty(agc)
% %         agcExt=zeros(nbus,size(agc,2));
% %         agcExt(agc(:,1),:)=agc;
% %         
% %         jfgi=[jfgi;dpgIdx(freqKeptTag==1)];
% %         jfgj=[jfgj;fIngIdx(freqKeptTag==1)];%f
% %         jfgv=[jfgv;-agcExt(freqKeptTag==1,4)];
% %         
% %         numSynOnBus=accumarray(syn(:,1),1,[nbus,1]);
% %         jffi=[jffi;dpgIdx(syn(:,1))];
% %         jffj=[jffj;omegaIdx];%f
% %         jffv=[jffv;-agcExt(syn(:,1),4)./numSynOnBus(syn(:,1))];
% %     end
% %     
% %     Jff=dt*sparse(jffi,jffj,jffv,nState,nState)/2-speye(nState,nState);
% %     
% %     % construct Jfg
% %     Jfg=dt*sparse(jfgi,jfgj,jfgv,nState,3*nbus)/2; %P,Q,f
% %     
% %     jac=[Jff(stateFilterTag==1,stateFilterTag==1),Jfg(stateFilterTag==1,jggColTag==1);Jgf(jggRowTag==1,stateFilterTag==1),Jgg(jggRowTag==1,jggColTag==1)];
% %     
% %     correction=-jac\[diffF(stateFilterTag==1);diffGx(jggRowTag==1)];
% %     nStateCorr=sum(stateFilterTag);
% %     xTemp(stateFilterTag==1)=xTemp(stateFilterTag==1)+correction(1:nStateCorr);
% %     corrFull=zeros(3*nbus,1);
% %     corrFull(jggRowTag==1)=correction((nStateCorr+1):end);
% %     Vtemp=Vtemp+corrFull(1:nbus)+1j*corrFull((nbus+1):(2*nbus));
% %     xTemp(vIdx)=Vtemp;
% %     xTemp(fIdx)=corrFull(fIngIdx);
% %     xTemp=adjustAlgebraic(SysDataxx,xTemp);
% % end
% % 
% % dxc=(xTemp-xt)/dt;
% % SysParaxx=foldSysPara([],[],[],[],[],[],[],[],Ytr0+Ytr1*(t+dt),[],Ysh0+Ysh1*(t+dt),[],[VspSq2(:,1)+VspSq2(:,2)*(t+dt),zeros(size(VspSq2,1),1)],[],[],[],[]);
% % diff=checkEquationBalanceSynDyn(SysDataxx,SysParaxx,xTemp,dxc);
% % 
% % end