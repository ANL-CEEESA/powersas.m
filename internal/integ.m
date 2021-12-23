function dx=integ(SysData,SysPara,xt,dt)
% [INTERNAL] Perform an integration
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
% FUNCTION integ
% CAUTION: this function is used internally. Do not call it directly unless you know what you are doing.
%
% INPUT
%   SysData - System data for simulation
%   SysPara - Parameters representing the events happening in the system
%	xt - Initial state
%	dt - Time step length
%
% OUTPUT
%	dx - Incremental of the system state
%
%
[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
[nState,idxs]...
    =getIndexDyn(SysData);
[V,Q,s,d,w,eq1,eq2,ed1,ed2,psid,psiq,Pg,Ef,Vavrm,Vavrr,Vavrf,Vavrref,tgovg,tgovm,tgovmech,f,dpg,qplt,vg]=unfoldX(xt,SysData);
nbus=size(bus,1);

[nIslands,islands,refs]=searchIslands(bus(:,1),line(:,1:2));

% Determine the frequency model of each island
freqTypeTag=zeros(nIslands,1);%0:sw,1:syn,2:steady-state f
freqKeptTag=zeros(nbus,1);
frefs=refs;
fswTag=zeros(nbus,1);
fsynTag=zeros(nbus,1);
fswTag(sw(:,1))=1;
fswTagxD=fswTag;
fsynTag(syn(:,1))=1;
D0=imag(V);
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

pVec=zeros(nbus,1);
qVec=zeros(nbus,1);

pVec(pq(:,1))=pVec(pq(:,1))+pq(:,4);
qVec(pq(:,1))=qVec(pq(:,1))+pq(:,5);
if ~isempty(zip)%zipMode=0, account for the PQ components in ZIP loads
    pVec=pVec+accumarray(zip(:,1),zip(:,7).*zip(:,12),[nbus,1]);
    qVec=qVec+accumarray(zip(:,1),zip(:,10).*zip(:,12),[nbus,1]);
end
SVec=pVec+1j*qVec;

[Ig,Im,Ii,Is]=calcCurrent(xt,syn,ind,zip,SVec,idxs.vIdx,idxs.sIdx,idxs.deltaIdx,idxs.omegaIdx,idxs.eq1Idx,idxs.eq2Idx,idxs.ed1Idx,idxs.ed2Idx,idxs.psidIdx,idxs.psiqIdx,idxs.efIdx);
[pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,~,~,~,~,~,Tmech1,Varref1,Ef1,Pm1,Eq11]=unfoldSysPara(SysPara);
%
dx=zeros(size(xt));
%
% d=xt(deltaIdx);
% w=xt(omegaIdx);
% eq1=xt(eq1Idx);
% eq2=xt(eq2Idx);
% ed1=xt(ed1Idx);
% ed2=xt(ed2Idx);
% psid=xt(psidIdx);
% psiq=xt(psiqIdx);
% V=xt(vIdx);
% s=xt(sIdx);

if ~isempty(ind)
    nInd=size(ind,1);
    indIdx=ind(:,1);
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
    Igx=real(Ig);
    Igy=imag(Ig);
    Igd=Igx.*sin(d)-Igy.*cos(d);
    Igq=Igx.*cos(d)+Igy.*sin(d);
    Vgx=real(V(synIdx));
    Vgy=imag(V(synIdx));
    Vgd=Vgx.*sin(d)-Vgy.*cos(d);
    Vgq=Vgx.*cos(d)+Vgy.*sin(d);
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
    Pe=Vgq.*Igq+Vgd.*Igd+Rga.*(Igq.*Igq+Igd.*Igd);
    
    for i=1:nSyn
        if model(i)==8
            dx(idxs.deltaIdx(i))=wgb(i)*w(i)*dt;
            dx(idxs.omegaIdx(i))=(Pg(i)-Pe(i)-Dg(i)*w(i))/Mg(i)*dt;
            dx(idxs.psidIdx(i))=wgb(i)*(Rga(i)*Igd(i)+psiq(i)+Vgd(i))*dt;
            dx(idxs.psiqIdx(i))=wgb(i)*(Rga(i)*Igq(i)-psid(i)+Vgq(i))*dt;
            dx(idxs.eq1Idx(i))=(-eq1(i)-(Xgd(i)-Xgd1(i)-gammad(i))*Igd(i)+(1-TgAA(i)/Tgd1(i))*Ef(i))/Tgd1(i)*dt;
            dx(idxs.ed1Idx(i))=(-ed1(i)+(Xgq(i)-Xgq1(i)-gammaq(i))*Igq(i))/Tgq1(i)*dt;
            dx(idxs.eq2Idx(i))=(-eq2(i)+eq1(i)-(Xgd1(i)-Xgd2(i)+gammad(i))*Igd(i)+TgAA(i)/Tgd1(i)*Ef(i))/Tgd2(i)*dt;
            dx(idxs.ed2Idx(i))=(-ed2(i)+ed1(i)+(Xgq1(i)-Xgq2(i)+gammaq(i))*Igq(i))/Tgq2(i)*dt;
        elseif model(i)==6
            dx(idxs.deltaIdx(i))=wgb(i)*w(i)*dt;
            dx(idxs.omegaIdx(i))=(Pg(i)-Pe(i)-Dg(i)*w(i))/Mg(i)*dt;
            dx(idxs.eq1Idx(i))=(-eq1(i)-(Xgd(i)-Xgd1(i)-gammad(i))*Igd(i)+(1-TgAA(i)/Tgd1(i))*Ef(i))/Tgd1(i)*dt;
            dx(idxs.ed1Idx(i))=(-ed1(i)+(Xgq(i)-Xgq1(i)-gammaq(i))*Igq(i))/Tgq1(i)*dt;
            dx(idxs.eq2Idx(i))=(-eq2(i)+eq1(i)-(Xgd1(i)-Xgd2(i)+gammad(i))*Igd(i)+TgAA(i)/Tgd1(i)*Ef(i))/Tgd2(i)*dt;
            dx(idxs.ed2Idx(i))=(-ed2(i)+ed1(i)+(Xgq1(i)-Xgq2(i)+gammaq(i))*Igq(i))/Tgq2(i)*dt;
        elseif model(i)==5
            dx(idxs.deltaIdx(i))=wgb(i)*w(i)*dt;
            dx(idxs.omegaIdx(i))=(Pg(i)-Pe(i)-Dg(i)*w(i))/Mg(i)*dt;
            dx(idxs.eq1Idx(i))=(-eq1(i)-(Xgd(i)-Xgd1(i)-gammad(i))*Igd(i)+(1-TgAA(i)/Tgd1(i))*Ef(i))/Tgd1(i)*dt;
            dx(idxs.eq2Idx(i))=(-eq2(i)+eq1(i)-(Xgd1(i)-Xgd2(i)+gammad(i))*Igd(i)+TgAA(i)/Tgd1(i)*Ef(i))/Tgd2(i)*dt;
            dx(idxs.ed2Idx(i))=(-ed2(i)+(Xgq(i)-Xgq2(i))*Igq(i))/Tgq2(i)*dt;
        elseif model(i)==4
            dx(idxs.deltaIdx(i))=wgb(i)*w(i)*dt;
            dx(idxs.omegaIdx(i))=(Pg(i)-Pe(i)-Dg(i)*w(i))/Mg(i)*dt;
            dx(idxs.eq1Idx(i))=(-eq1(i)-(Xgd(i)-Xgd1(i))*Igd(i)+Ef(i))/Tgd1(i)*dt;
            dx(idxs.ed1Idx(i))=(-ed1(i)+(Xgq(i)-Xgq1(i))*Igq(i))/Tgq1(i)*dt;
        elseif model(i)==3
            dx(idxs.deltaIdx(i))=wgb(i)*w(i)*dt;
            dx(idxs.omegaIdx(i))=(Pg(i)-Pe(i)-Dg(i)*w(i))/Mg(i)*dt;
            dx(idxs.eq1Idx(i))=(-eq1(i)-(Xgd(i)-Xgd1(i))*Igd(i)+Ef(i))/Tgd1(i)*dt;
        elseif model(i)==2
            dx(idxs.deltaIdx(i))=wgb(i)*w(i)*dt;
            dx(idxs.omegaIdx(i))=(Pg(i)-Pe(i)-Dg(i)*w(i))/Mg(i)*dt;
%             dx(eq1Idx(i))=Eq11(i)*dt;
        end
    end
    dx(idxs.efIdx)=Ef1*dt;
    dx(idxs.pgIdx)=Pm1*dt;
end

if ~isempty(ind)
    yind=1./(Zm1+Zmm.*(Rm2+1j*Xm2.*s)./(Rm2+s.*(Zmm+1j*Xm2)));
    il=V(indIdx).*yind;
    ir=il-(V(indIdx)-il.*Zm1)./Zmm;
    dx(idxs.sIdx)=dt*((Tm0+Tm1.*s+Tm2.*s.*s)-ir.*conj(ir).*Rm2./s)/2./Hm;
end

Vmag=abs(V);

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
    dx(idxs.vavrmIdx)=(VmagAvr-Vavrm)./Tavrr*dt;
    dx(idxs.vavrrIdx)=(muavr0.*(1-Tavr1./Tavr2).*(Vavrref-Vavrm)-Vavrr)./Tavr2*dt;
    dx(idxs.vavrfIdx)=((Vavrr+muavr0.*Tavr1./Tavr2.*(Vavrref-Vavrm)+vavrf0).*VmagAvr./Vavr0-Vavrf)./Tavre*dt;
    dx(idxs.vavrrefIdx)=Varref1*dt;
end

%TODO: AGC model

if ~isempty(tg)
    tgIdx=tg(:,1);
    wtgref=tg(:,3);
    Rtg=tg(:,4);
    Ttgmax=tg(:,5);
    Ttgmin=tg(:,6);
    Ttg2=tg(:,7);
    Ttg1=tg(:,8);
    
    dx(idxs.tgovgIdx)=((1-Ttg1./Ttg2).*(wtgref-w(tgIdx))./Rtg-tgovg)./Ttg2*dt;
    dx(idxs.tmechIdx)=Tmech1*dt;
end

nTg=size(tg,1);
fx=zeros(nbus,1);
synTag=zeros(nbus,1);
synTag(syn(:,1))=1:nSyn;
numSynOnBus=accumarray(syn(:,1),1,[nbus,1]);

if ~isempty(agc)
    agcExt=zeros(nbus,size(agc,2));
    agcExt(agc(:,1),:)=agc;
    dx(idxs.dpgIdx)=-f.*agcExt(:,4)*dt;  
    dx(idxs.tmechIdx)=dx(idxs.tmechIdx)+dx(idxs.dpgIdx(syn(tg(:,1),1)));
    dx(idxs.pgIdx)=dx(idxs.pgIdx)+dx(idxs.dpgIdx(syn(:,1)));
end

end
