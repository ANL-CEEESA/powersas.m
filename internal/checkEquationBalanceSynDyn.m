function [diff,psw,qsw]=checkEquationBalanceSynDyn(SysData,SysPara,x,dx)
% Checking the equation balance in dynamic simulation
%% FUNCTION checkEquationBalanceSynDyn
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
%   SysPara - Parameters representing the events happening in the system
%	x - System states
%	dx - Time derivative of system states
%	
%%
[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
nbus=size(bus,1);
nline=size(line,1);
[V,Q,s,d,w,eq1,eq2,ed1,ed2,psid,psiq,Pm,Ef,Vavrm,Vavrr,Vavrf,Vavrref,tgovg,tgovm,tgovmech,f,dpg,qplt,vg]=unfoldX(x,SysData);
[~,~,ds,dd,dw,deq1,deq2,ded1,ded2,dpsid,dpsiq,dPm,~,dVavrm,dVavrr,dVavrf,~,dtgovg,~,dtgovmech,df,ddpg,dqplt,dvg]=unfoldX(dx,SysData);
[pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,~,Ysh0,Ysh1,Vsp2,~,~,~,~,Tmech1,Varref1,Ef1,Pm1,Eq11]=unfoldSysPara(SysPara);

if isfield(SysPara,'nIslands')&&isfield(SysPara,'islands')&&isfield(SysPara,'refs')
    nIslands=SysPara.nIslands;islands=SysPara.islands;refs=SysPara.refs;
else
    [nIslands,islands,refs]=searchIslands(bus(:,1),line(:,1:2));
end

if ~isempty(zip)%zipMode=0
    Ysh0=Ysh0+accumarray(zip(:,1),(zip(:,5)+1j*zip(:,8)).*zip(:,12),[nbus,1]);
end
if ~isempty(shunt)    
    yShunt=zeros(nbus,1);
    yShunt(shunt(:,1))=shunt(:,5)+1j*shunt(:,6);
    Ysh0=Ysh0+yShunt;
end
Y=Ytr0+sparse(1:nbus,1:nbus,Ysh0,nbus,nbus);
%
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

D0=imag(V);

% Determine the frequency model of each island
freqTypeTag=zeros(nIslands,1);%0:sw,1:syn,2:steady-state f
freqKeptTag=zeros(nbus,1);
% frefs=refs;
fswTag=zeros(nbus,1);
fsynTag=zeros(nbus,1);
fswTag(sw(:,1))=1;
% fswTagxD=fswTag;
fsynTag(syn(:,1))=1;
for isl=1:nIslands
    if isempty(find(fswTag(islands==isl)==1, 1))
        if isempty(find(fsynTag(islands==isl)==1, 1))
            freqTypeTag(isl)=2;
            busesInIsland=find(islands==isl);
%             [~,imin]=min(abs(D0(busesInIsland)));
%             frefs(isl)=busesInIsland(imin(1));
%             fswTagxD(frefs(isl))=1;
            freqKeptTag(busesInIsland)=1;
        else
            freqTypeTag(isl)=1;
        end
    end
end
% freqKeptTagxRef=freqKeptTag;
% freqKeptTagxRef(frefs)=0;
% nFreqKept=sum(freqKeptTag);
pVec=zeros(nbus,1);
qVec=zeros(nbus,1);

pVec(pv(:,1))=pVec(pv(:,1))+pv(:,4);
pVec(pq(:,1))=pVec(pq(:,1))-pq(:,4);
qVec(pq(:,1))=qVec(pq(:,1))-pq(:,5);
if ~isempty(zip)%zipMode=0
    pVec=pVec-accumarray(zip(:,1),zip(:,7).*zip(:,12),[nbus,1]);
    qVec=qVec-accumarray(zip(:,1),zip(:,10).*zip(:,12),[nbus,1]);
end
% qVec(pv(:,1))=0;

IInj=Y*V;
if ~isempty(zip)%zipMode=0
    IInj=IInj+accumarray(zip(:,1),(zip(:,6)-1j*zip(:,9)).*zip(:,12).*V(zip(:,1))./abs(V(zip(:,1))),[nbus,1]);
end
SInjRHS=V.*conj(IInj);

nSyn=size(syn,1);
diffdSyn=zeros(nSyn,1);
diffwSyn=zeros(nSyn,1);
diffeq1Syn=zeros(nSyn,1);
diffeq2Syn=zeros(nSyn,1);
diffed1Syn=zeros(nSyn,1);
diffed2Syn=zeros(nSyn,1);
diffpsiqSyn=zeros(nSyn,1);
diffpsidSyn=zeros(nSyn,1);
if ~isempty(syn)    
    synIdx=syn(:,1);    
    CG0=real(V(synIdx));
    DG0=imag(V(synIdx));
    IGd=zeros(nSyn,1);
    IGq=zeros(nSyn,1);
    cosd=cos(d);
    sind=sin(d);
    VGd=sind.*CG0-cosd.*DG0;
    VGq=cosd.*CG0+sind.*DG0;
    
    wgb=syn(:,4);
    modSyn=syn(:,5);
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
    modelTag=accumarray(modSyn,ones(nSyn,1),[8,1]);
    
    if modelTag(8)>0
        IGd(modSyn==8)=(eq2(modSyn==8)-psid(modSyn==8))./Xgd2(modSyn==8);
        IGq(modSyn==8)=(-ed2(modSyn==8)-psiq(modSyn==8))./Xgq2(modSyn==8);
        diffdSyn(modSyn==8)=dd(modSyn==8)-(wgb(modSyn==8).*w(modSyn==8));
        diffwSyn(modSyn==8)=dw(modSyn==8)-(Pm(modSyn==8)-(psid(modSyn==8).*IGq(modSyn==8)-psiq(modSyn==8).*IGd(modSyn==8))-Dg(modSyn==8).*w(modSyn==8))./Mg(modSyn==8);
        diffpsidSyn(modSyn==8)=dpsid(modSyn==8)-wgb(modSyn==8).*(Rga.*IGd(modSyn==8)+psiq(modSyn==8)+VGd(modSyn==8));
        diffpsiqSyn(modSyn==8)=dpsiq(modSyn==8)-wgb(modSyn==8).*(Rga.*IGq(modSyn==8)-psid(modSyn==8)+VGq(modSyn==8));
        diffeq1Syn(modSyn==8)=deq1(modSyn==8)-(-eq1(modSyn==8)-(Xgd(modSyn==8)-Xgd1(modSyn==8)-gammad(modSyn==8)).*IGd(modSyn==8)+(1-TgAA(modSyn==8)./Tgd1(modSyn==8)).*Ef(modSyn==8))./Tgd1(modSyn==8);
        diffed1Syn(modSyn==8)=ded1(modSyn==8)-(-ed1(modSyn==8)+(Xgq(modSyn==8)-Xgq1(modSyn==8)-gammaq(modSyn==8)).*IGq(modSyn==8))./Tgq1(modSyn==8);
        diffeq2Syn(modSyn==8)=deq2(modSyn==8)-(-eq2(modSyn==8)+eq1(modSyn==8)-(Xgd1(modSyn==8)-Xgd2(modSyn==8)+gammad(modSyn==8)).*IGd(modSyn==8)+TgAA(modSyn==8)./Tgd1(modSyn==8).*Ef(modSyn==8))./Tgd2(modSyn==8);
        diffed2Syn(modSyn==8)=ded2(modSyn==8)-(-ed2(modSyn==8)+ed1(modSyn==8)+(Xgq1(modSyn==8)-Xgq2(modSyn==8)+gammaq(modSyn==8)).*IGq(modSyn==8))./Tgq2(modSyn==8);
    end    
    if modelTag(6)>0
        IGd(modSyn==6)=((ed2(modSyn==6)-VGd(modSyn==6)).*Rga(modSyn==6)+(eq2(modSyn==6)-VGq(modSyn==6)).*Xgq2(modSyn==6))./(Rga(modSyn==6).*Rga(modSyn==6)+Xgd2(modSyn==6).*Xgq2(modSyn==6));
        IGq(modSyn==6)=(-(ed2(modSyn==6)-VGd(modSyn==6)).*Xgd2(modSyn==6)+(eq2(modSyn==6)-VGq(modSyn==6)).*Rga(modSyn==6))./(Rga(modSyn==6).*Rga(modSyn==6)+Xgd2(modSyn==6).*Xgq2(modSyn==6));
        diffdSyn(modSyn==6)=dd(modSyn==6)-(wgb(modSyn==6).*w(modSyn==6));
        diffwSyn(modSyn==6)=dw(modSyn==6)-(Pm(modSyn==6)-(VGq(modSyn==6).*IGq(modSyn==6)+VGd(modSyn==6).*IGd(modSyn==6)+Rga(modSyn==6).*(IGq(modSyn==6).*IGq(modSyn==6)+IGd(modSyn==6).*IGd(modSyn==6)))-Dg(modSyn==6).*w(modSyn==6))./Mg(modSyn==6);
        diffeq1Syn(modSyn==6)=deq1(modSyn==6)-(-eq1(modSyn==6)-(Xgd(modSyn==6)-Xgd1(modSyn==6)-gammad(modSyn==6)).*IGd(modSyn==6)+(1-TgAA(modSyn==6)./Tgd1(modSyn==6)).*Ef(modSyn==6))./Tgd1(modSyn==6);
        diffed1Syn(modSyn==6)=ded1(modSyn==6)-(-ed1(modSyn==6)+(Xgq(modSyn==6)-Xgq1(modSyn==6)-gammaq(modSyn==6)).*IGq(modSyn==6))./Tgq1(modSyn==6);
        diffeq2Syn(modSyn==6)=deq2(modSyn==6)-(-eq2(modSyn==6)+eq1(modSyn==6)-(Xgd1(modSyn==6)-Xgd2(modSyn==6)+gammad(modSyn==6)).*IGd(modSyn==6)+TgAA(modSyn==6)./Tgd1(modSyn==6).*Ef(modSyn==6))./Tgd2(modSyn==6);
        diffed2Syn(modSyn==6)=ded2(modSyn==6)-(-ed2(modSyn==6)+ed1(modSyn==6)+(Xgq1(modSyn==6)-Xgq2(modSyn==6)+gammaq(modSyn==6)).*IGq(modSyn==6))./Tgq2(modSyn==6);
    end
    if modelTag(5)>0
        IGd(modSyn==5)=((ed2(modSyn==5)-VGd(modSyn==5)).*Rga(modSyn==5)+(eq2(modSyn==5)-VGq(modSyn==5)).*Xgq2(modSyn==5))./(Rga(modSyn==5).*Rga(modSyn==5)+Xgd2(modSyn==5).*Xgq2(modSyn==5));
        IGq(modSyn==5)=(-(ed2(modSyn==5)-VGd(modSyn==5)).*Xgd2(modSyn==5)+(eq2(modSyn==5)-VGq(modSyn==5)).*Rga(modSyn==5))./(Rga(modSyn==5).*Rga(modSyn==5)+Xgd2(modSyn==5).*Xgq2(modSyn==5));
        diffdSyn(modSyn==5)=dd(modSyn==5)-(wgb(modSyn==5).*w(modSyn==5));
        diffwSyn(modSyn==5)=dw(modSyn==5)-(Pm(modSyn==5)-(VGq(modSyn==5).*IGq(modSyn==5)+VGd(modSyn==5).*IGd(modSyn==5)+Rga(modSyn==5).*(IGq(modSyn==5).*IGq(modSyn==5)+IGd(modSyn==5).*IGd(modSyn==5)))-Dg(modSyn==5).*w(modSyn==5))./Mg(modSyn==5);
        diffeq1Syn(modSyn==5)=deq1(modSyn==5)-(-eq1(modSyn==5)-(Xgd(modSyn==5)-Xgd1(modSyn==5)-gammad(modSyn==5)).*IGd(modSyn==5)+(1-TgAA(modSyn==5)./Tgd1(modSyn==5)).*Ef(modSyn==5))./Tgd1(modSyn==5);
        diffeq2Syn(modSyn==5)=deq2(modSyn==5)-(-eq2(modSyn==5)+eq1(modSyn==5)-(Xgd1(modSyn==5)-Xgd2(modSyn==5)+gammad(modSyn==5)).*IGd(modSyn==5)+TgAA(modSyn==5)./Tgd1(modSyn==5).*Ef(modSyn==5))./Tgd2(modSyn==5);
        diffed2Syn(modSyn==5)=ded2(modSyn==5)-(-ed2(modSyn==5)+ed1(modSyn==5)+(Xgq(modSyn==5)-Xgq2(modSyn==5)).*IGq(modSyn==5))./Tgq2(modSyn==5);
    end
    if modelTag(4)>0
        IGd(modSyn==4)=((ed1(modSyn==4)-VGd(modSyn==4)).*Rga(modSyn==4)+(eq1(modSyn==4)-VGq(modSyn==4)).*Xgq1(modSyn==4))./(Rga(modSyn==4).*Rga(modSyn==4)+Xgd1(modSyn==4).*Xgq1(modSyn==4));
        IGq(modSyn==4)=(-(ed1(modSyn==4)-VGd(modSyn==4)).*Xgd1(modSyn==4)+(eq1(modSyn==4)-VGq(modSyn==4)).*Rga(modSyn==4))./(Rga(modSyn==4).*Rga(modSyn==4)+Xgd1(modSyn==4).*Xgq1(modSyn==4));
        diffdSyn(modSyn==4)=dd(modSyn==4)-(wgb(modSyn==4).*w(modSyn==4));
        diffwSyn(modSyn==4)=dw(modSyn==4)-(Pm(modSyn==4)-(VGq(modSyn==4).*IGq(modSyn==4)+VGd(modSyn==4).*IGd(modSyn==4)+Rga(modSyn==4).*(IGq(modSyn==4).*IGq(modSyn==4)+IGd(modSyn==4).*IGd(modSyn==4)))-Dg(modSyn==4).*w(modSyn==4))./Mg(modSyn==4);
        diffeq1Syn(modSyn==4)=deq1(modSyn==4)-(-eq1(modSyn==4)-(Xgd(modSyn==4)-Xgd1(modSyn==4)).*IGd(modSyn==4)+Ef(modSyn==4))./Tgd1(modSyn==4);
        diffed1Syn(modSyn==4)=ded1(modSyn==4)-(-ed1(modSyn==4)+(Xgq(modSyn==4)-Xgq1(modSyn==4)).*IGq(modSyn==4))./Tgq1(modSyn==4);
    end
    if modelTag(3)>0
        IGd(modSyn==3)=((-VGd(modSyn==3)).*Rga(modSyn==3)+(eq1(modSyn==3)-VGq(modSyn==3)).*Xgq(modSyn==3))./(Rga(modSyn==3).*Rga(modSyn==3)+Xgd1(modSyn==3).*Xgq(modSyn==3));
        IGq(modSyn==3)=(-(-VGd(modSyn==3)).*Xgd1(modSyn==3)+(eq1(modSyn==3)-VGq(modSyn==3)).*Rga(modSyn==3))./(Rga(modSyn==3).*Rga(modSyn==3)+Xgd1(modSyn==3).*Xgq(modSyn==3));
        diffdSyn(modSyn==3)=dd(modSyn==3)-(wgb(modSyn==3).*w(modSyn==3));
        diffwSyn(modSyn==3)=dw(modSyn==3)-(Pm(modSyn==3)-(VGq(modSyn==3).*IGq(modSyn==3)+VGd(modSyn==3).*IGd(modSyn==3)+Rga(modSyn==3).*(IGq(modSyn==3).*IGq(modSyn==3)+IGd(modSyn==3).*IGd(modSyn==3)))-Dg(modSyn==3).*w(modSyn==3))./Mg(modSyn==3);
        diffeq1Syn(modSyn==3)=deq1(modSyn==3)-(-eq1(modSyn==3)-(Xgd(modSyn==3)-Xgd1(modSyn==3)).*IGd(modSyn==3)+Ef(modSyn==3))./Tgd1(modSyn==3);
    end
    if modelTag(2)>0
        IGd(modSyn==2)=((-VGd(modSyn==2)).*Rga(modSyn==2)+(Ef(modSyn==2)-VGq(modSyn==2)).*Xgq(modSyn==2))./(Rga(modSyn==2).*Rga(modSyn==2)+Xgd(modSyn==2).*Xgq(modSyn==2));
        IGq(modSyn==2)=(-(-VGd(modSyn==2)).*Xgd(modSyn==2)+(Ef(modSyn==2)-VGq(modSyn==2)).*Rga(modSyn==2))./(Rga(modSyn==2).*Rga(modSyn==2)+Xgd(modSyn==2).*Xgq(modSyn==2));
        diffdSyn(modSyn==2)=dd(modSyn==2)-(wgb(modSyn==2).*w(modSyn==2));
        diffwSyn(modSyn==2)=dw(modSyn==2)-(Pm(modSyn==2)-(VGq(modSyn==2).*IGq(modSyn==2)+VGd(modSyn==2).*IGd(modSyn==2)+Rga(modSyn==2).*(IGq(modSyn==2).*IGq(modSyn==2)+IGd(modSyn==2).*IGd(modSyn==2)))-Dg(modSyn==2).*w(modSyn==2))./Mg(modSyn==2);
    end
     
    JG=sind.*IGd+cosd.*IGq;
    KG=-cosd.*IGd+sind.*IGq;
    
    SInjRHS_syn=-accumarray(synIdx,CG0.*JG+DG0.*KG+1j*(DG0.*JG-CG0.*KG),[nbus,1]);
    
    SInjRHS=SInjRHS+SInjRHS_syn;
end
diffSyn=[diffdSyn;diffwSyn;diffpsidSyn;diffpsiqSyn;diffeq1Syn;diffed1Syn;diffeq2Syn;diffed2Syn];

% qVec(busType~=0)=0;

if ~isempty(ind)
    nInd=size(ind,1);
    indIdx=ind(:,1);
    Rm1=ind(:,7);
    Xm1=ind(:,8);
    Zm1=ind(:,7)+1j*ind(:,8);
    Xmm=ind(:,13);
    Zme=1j*ind(:,13);
    Rm2=ind(:,9);
    Xm2=ind(:,10);
    Tm0=ind(:,15)+ind(:,16)+ind(:,17);
    Tm1=-ind(:,16)-2*ind(:,17);
    Tm2=ind(:,17);
    Hm=ind(:,14);
    
    Zind=Zm1+Zme.*(Rm2+1j*Xm2.*s)./(Rm2+(1j*Xm2+Zme).*s);
    ILind=V(indIdx)./Zind;
    SInjRHS=SInjRHS+V.*conj(accumarray(indIdx,ILind,[nbus,1]));
end

diffSInj=SInjRHS-(pVec+1j*qVec);
diffSInj=diffSInj-1j*Q;
ssw=diffSInj(busType==2);
psw=real(ssw);qsw=imag(ssw);
diffSInj(busType==2)=0;

diffsInd=zeros(0,1);
if ~isempty(ind)
    Vind=V(indIdx);
    VEind=Vind.*Zme./(Zm1+Zme);
    IRs=VEind.*s./(s.*Zm1.*Zme./(Zm1+Zme)+Rm2+1j*Xm2.*s);
    sIndLHS=ds*2.*Hm.*s;
    diffsInd=real(sIndLHS-(Tm0+s.*(Tm1+s.*Tm2)).*s+IRs.*conj(IRs).*Rm2);
end

Vmag=abs(V);

nExc=size(exc,1);
diffExcVm=zeros(nExc,1);
diffExcVr=zeros(nExc,1);
diffExcVf=zeros(nExc,1);
diffExcEf=zeros(nExc,1);
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
    
    VmagAvr=Vmag(synIdx(excIdx));
    diffExcVm=dVavrm-(VmagAvr-Vavrm)./Tavrr;
    diffExcVr=dVavrr-(muavr0.*(1-Tavr1./Tavr2).*(Vavrref-Vavrm)-Vavrr)./Tavr2;
    diffExcVf=dVavrf-((Vavrr+muavr0.*Tavr1./Tavr2.*(Vavrref-Vavrm)+vavrf0).*VmagAvr./Vavr0-Vavrf)./Tavre;
    
    Vavrf=real(Vavrf);
    Efx=Vavrf;
    tavrMaxDiff=Vavrf-VavrMax;
    tavrMinDiff=Vavrf-VavrMin;
    Efx(tavrMaxDiff>0)=VavrMax(tavrMaxDiff>0);    
    Efx(tavrMinDiff<0)=VavrMin(tavrMinDiff<0);
    
    diffExcEf=Ef(excIdx)-Efx;
end
diffExc=[diffExcVm;diffExcVr;diffExcVf;diffExcEf];

nTg=size(tg,1);
diffAgcm=zeros(nbus,1);
diffAgcf=zeros(nSyn,1);
diffAgcftg=zeros(nTg,1);
fx=zeros(nbus,1);
synTag=zeros(nbus,1);
synTag(syn(:,1))=1:nSyn;
numSynOnBus=accumarray(syn(:,1),1,[nbus,1]);
for islIdx=1:nIslands
    busIsland=find(islands==islIdx);
    synTagIsland=synTag(busIsland);
    wIsland=w(synTagIsland(synTagIsland~=0));
    if ~isempty(wIsland)
        fx(busIsland)=mean(wIsland); % note that here the freq can be different
    else
        % TODO: steady-state model
    end
    
    %         synIslIdx=syn(synTagIsland(synTagIsland~=0),1);
    %         diffAgcf=dPm+fx(synIslIdx).*agcExt(synIslIdx,4);
end
diffFreq=fx-f;
if ~isempty(agc)
    if isempty(Pm1);Pm1=zeros(nSyn,1);end 
    if isempty(Tmech1);Tmech1=zeros(nTg,1);end 
    agcExt=zeros(nbus,size(agc,2));
    agcExt(agc(:,1),:)=agc;
    diffAgcm=ddpg+f.*agcExt(:,4);
    if ~isempty(syn)
        diffAgcf=dPm-ddpg(syn(:,1))./numSynOnBus(syn(:,1))-Pm1;
        if ~isempty(tg)
            diffAgcf(tg(:,1))=0;
            diffAgcftg=dtgovmech-ddpg(syn(tg(:,1),1))./numSynOnBus(syn(tg(:,1),1))-Tmech1;
        end
    end
    fdk=agcExt(:,2)+agcExt(:,3); %1/R+D
    diffSInj=diffSInj+freqKeptTag.*(fdk.*f-dpg);
end
diffAgc=[diffFreq(freqKeptTag==0);diffAgcm;diffAgcf;diffAgcftg];

diffTgtg=zeros(nTg,1);
diffTgtm=zeros(nTg,1);
diffTgPm=zeros(nTg,1);
if ~isempty(tg)    
    tgIdx=tg(:,1);    
    wtgref=tg(:,3);
    Rtg=tg(:,4);
    Ttgmax=tg(:,5);
    Ttgmin=tg(:,6);
    Ttg2=tg(:,7);
    Ttg1=tg(:,8);
    
    diffTgtg=dtgovg-((1-Ttg1./Ttg2).*(wtgref-w(tgIdx))./Rtg-tgovg)./Ttg2;
    diffTgtm=tgovm-(tgovg+Ttg1./Ttg2.*(wtgref-w(tgIdx))./Rtg+tgovmech);
    
    tgovm=real(tgovm);
    Pmx=tgovm;
    tgovMaxDiff=tgovm(:,1)-Ttgmax;
    tgovMinDiff=tgovm(:,1)-Ttgmin;    
    Pmx(tgovMaxDiff>0)=Ttgmax(tgovMaxDiff>0);
    Pmx(tgovMinDiff<0)=Ttgmin(tgovMinDiff<0); 
    diffTgPm=Pm(tgIdx)-Pmx;
end
diffTg=[diffTgtg;diffTgtm;diffTgPm];

diffv=sum(Vsp2,2)-V.*conj(V);

diff=[real(diffSInj);imag(diffSInj);diffv(busType~=0);diffsInd;diffSyn;diffExc;diffTg;diffAgc];
end