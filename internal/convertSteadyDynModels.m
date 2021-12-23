function [SysDataNew,xNew]=convertSteadyDynModels(SysData,SysDataOrig,x,deltax,xOrig,SysPara,convMode)
% Conversion between QSS and dynamic model
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
% 	SysData - System data for conversion
%	SysDataOrig - Original system data for reference (for guiding conversion from QSS to dynamic model)
%	x - Current system states
%	deltax - Time derivative of system states
%	xOrig - Original system state
%	SysPara - 
%	convMode
%		0 - QSS -> Dynamic
%		1 - Dynamic -> QSS
% 	
% OUTPUT
%	SysDataNew - Converted system data
%	xNew - new system state
%	
% 


% [busBase,swBase,pvBase,pqBase,shuntBase,lineBase,indBase,zipBase,synBase,excBase,tgBase,agcBase,cacBase,clusterBase]=unfoldSysData(BaseSysData);
[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
[V,Q,s,d,w,eq1,eq2,ed1,ed2,psid,psiq,Pm,Ef,Vavrm,Vavrr,Vavrf,Vavrref,tgovg,tgovm,tgovmech,f,dpg,qplt,vg]=unfoldX(x,SysData);
[~,~,ds,dd,dw,deq1,deq2,ded1,ded2,dpsid,dpsiq,dPm,~,dVavrm,dVavrr,dVavrf,~,dtgovg,~,dtgovmech,df,ddpg,dqplt,dvg]=unfoldX(deltax,SysData);
[pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,...
    MatGV0,MatGV1,MatGRhs0,MatGRhs1,Tmech1,Varref1,Ef1,Pm1,Eq11,fault]=unfoldSysPara(SysPara);
[nIslands,islands]=searchIslands(bus(:,1),line(:,[1,2]));

if convMode==0 % To full dynamic model, i.e. convert all sw and pv to syn
    [buso,swo,pvo,pqo,shunto,lineo,indo,zipo,syno,exco,tgo,agco,caco,clustero]=unfoldSysData(SysDataOrig);% SysDataOrig must not be empty
    
    %
    nbus=size(bus,1);
    nline=size(line,1);
    [nIslands,islands,refs]=searchIslands(bus(:,1),line(:,1:2));
    yShunt=zeros(nbus,1);
    yShunt(shunt(:,1))=shunt(:,5)+1j*shunt(:,6);
    if ~isempty(zip)%zipMode=0
        yShunt=yShunt+accumarray(zip(:,1),(zip(:,5)+1j*zip(:,8)).*zip(:,12),[nbus,1]);
    end
    if isempty(Ytr0)
        [~,Ytr0,Ysh0]=getYMatrix(nbus,line,fault);
    end
    Ysh0=Ysh0+yShunt;
    Y=Ytr0+sparse(1:nbus,1:nbus,Ysh0,nbus,nbus);
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
    agcExt=zeros(nbus,4);
    if ~isempty(agc)
        agcExt(agc(:,1),:)=agc;
    end 
    fdk=agcExt(:,2)+agcExt(:,3); 
    
    D0=imag(V);
    % Determine the frequency model of each island
    freqTypeTag=zeros(nIslands,1);%0:sw,1:syn,2:steady-state f
    freqKeptTag=zeros(nIslands,1);
    frefs=refs;
    fswTag=zeros(nbus,1);
    fsynTag=zeros(nbus,1);
    fswTag(sw(:,1))=1;
    fswTagxD=fswTag;
    fsynTag(syn(:,1))=1;
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
    
%     pVec(pv(:,1))=pVec(pv(:,1))+pv(:,4);
    pVec(pq(:,1))=pVec(pq(:,1))-pq(:,4);
    qVec(pq(:,1))=qVec(pq(:,1))-pq(:,5);
    if ~isempty(zip)%zipMode=0
        pVec=pVec-accumarray(zip(:,1),zip(:,7).*zip(:,12),[nbus,1]);
        qVec=qVec-accumarray(zip(:,1),zip(:,10).*zip(:,12),[nbus,1]);
    end
    
    IInj=Y*V;
    if ~isempty(zip)%zipMode=0
        IInj=IInj+accumarray(zip(:,1),(zip(:,6)-1j*zip(:,9)).*zip(:,12).*V(zip(:,1))./abs(V(zip(:,1))),[nbus,1]);
    end
    SInjRHS=V.*conj(IInj);
    
    %     qVec(busType~=0)=0;
    
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
    %     diffSInj=diffSInj-1j*Q;
    %     diffSInj(busType==2)=0;
    sg=diffSInj;
    sgSyn=sg(syno(:,1));
    IgSyn=conj(sgSyn./V(syno(:,1)));
    w=f(syno(:,1));
    Qx=zeros(nbus,1);
    Eq=V(syno(:,1))+IgSyn.*(syno(:,7)+1j*syno(:,13));
    dg=atan2(imag(Eq),real(Eq));
    unitEq=Eq./abs(Eq);
    unitEd=-1j*unitEq;
    Id=unitEd.*(real(unitEd).*real(IgSyn)+imag(unitEd).*imag(IgSyn));
    Efx=Eq+1j*Id.*(syno(:,8)-syno(:,13));    
    
    synIdx=syno(:,1);
    Rga=syno(:,7);
    Xgd=syno(:,8);
    Xgd1=syno(:,9);
    Xgd2=syno(:,10);
    Xgq=syno(:,13);
    Xgq1=syno(:,14);
    Xgq2=syno(:,15);
    Vgx=real(V(synIdx));
    Vgy=imag(V(synIdx));
    Igx=real(IgSyn);
    Igy=imag(IgSyn);
    Vgd=Vgx.*sin(dg)-Vgy.*cos(dg);
    Vgq=Vgx.*cos(dg)+Vgy.*sin(dg);
    Igd=Igx.*sin(dg)-Igy.*cos(dg);
    Igq=Igx.*cos(dg)+Igy.*sin(dg);
    
    eq1=Vgq+Igd.*Xgd1+Igq.*Rga;
    ed1=Vgd+Igd.*Rga-Igq.*Xgq;
    
    SysDataNew=foldSysData(buso,swo,pvo,pq,shunt,lineo,ind,zip,syno,exco,tgo,agc,caco,clustero);
    [nState,idxs]...
        =getIndexDyn(SysDataNew);
    xx=zeros(nState,1);
    xx(idxs.vIdx)=V;
    xx(idxs.qIdx)=Qx;
    xx(idxs.sIdx)=s;
    xx(idxs.deltaIdx)=dg;
    xx(idxs.omegaIdx)=w;
    xx(idxs.efIdx)=abs(Efx);
    deltaxx=zeros(nState,1);
    deltaxx(idxs.omegaIdx)=df(syno(:,1));
%     deltaxx(tgovgIdx)=ddpg(syno(tgo(:,1),1))+tgo(:,8)./tgo(:,7)./tgo(:,4).*df(syno(tgo(:,1),1));
    deltaxx(idxs.tgovgIdx)=-(1-tgo(:,8)./tgo(:,7))./tgo(:,4).*df(syno(tgo(:,1),1));
    
    xNew=perpareInitialState(SysDataNew,xx,deltaxx);
    SysDataNew.agc=agco;
    xNew(idxs.fIdx)=f;
%     xNew(eq1Idx)=xOrig(eq1Idx);
%     xNew(eq2Idx)=xOrig(eq2Idx);
%     xNew(ed1Idx)=xOrig(ed1Idx);
%     xNew(ed2Idx)=xOrig(ed2Idx);
%     xNew(psidIdx)=xOrig(psidIdx);
%     xNew(psiqIdx)=xOrig(psiqIdx);
    xNew(idxs.vavrmIdx)=xOrig(idxs.vavrmIdx);
    xNew(idxs.vavrrIdx)=xOrig(idxs.vavrrIdx);
    xNew(idxs.vavrfIdx)=xOrig(idxs.vavrfIdx);
    xNew(idxs.vavrrefIdx)=xOrig(idxs.vavrrefIdx);
elseif convMode==1 % To QSS model, i.e. convert all syn to sw and pv (or pv if agc is on)
    [nState,idxs]...
        =getIndexDyn(SysData);
    
    nbus=size(bus,1);
    nline=size(line,1);
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
    agcExt=zeros(nbus,4);
    agcExt(1:nbus,1)=1:nbus;
    if ~isempty(agc)
        agcExt(agc(:,1),:)=agc;
    end  
    
    synIsland=islands(syn(:,1));
    synV=V(syn(:,1));
    
    dg=x(idxs.deltaIdx);
    wg=x(idxs.omegaIdx);
    eq1=x(idxs.eq1Idx);
    eq2=x(idxs.eq2Idx);
    ed1=x(idxs.ed1Idx);
    ed2=x(idxs.ed2Idx);
    psid=x(idxs.psidIdx);
    psiq=x(idxs.psiqIdx);
    vbus=x(idxs.vIdx);
    Ef=x(idxs.efIdx);
    
    if ~isempty(syn)
        nSyn=size(syn,1);
        synIdx=syn(:,1);
        Vgx=real(vbus(synIdx));
        Vgy=imag(vbus(synIdx));
        Vgd=Vgx.*sin(dg)-Vgy.*cos(dg);
        Vgq=Vgx.*cos(dg)+Vgy.*sin(dg);
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
        
        Igd=zeros(nSyn,1);
        Igq=zeros(nSyn,1);
        for i=1:nSyn
            if model(i)==8
                Igd(i)=(eq2(i)-psid(i))/Xgd2(i);
                Igq(i)=(-ed2(i)-psiq(i))/Xgq2(i);
            elseif model(i)==6
                Igd(i)=(Rga(i)*(ed2(i)-Vgd(i))+Xgq2(i)*(eq2(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd2(i)*Xgq2(i));
                Igq(i)=(-Xgd2(i)*(ed2(i)-Vgd(i))+Rga(i)*(eq2(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd2(i)*Xgq2(i));
            elseif model(i)==5
                Igd(i)=(Rga(i)*(ed2(i)-Vgd(i))+Xgq2(i)*(eq2(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd2(i)*Xgq2(i));
                Igq(i)=(-Xgd2(i)*(ed2(i)-Vgd(i))+Rga(i)*(eq2(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd2(i)*Xgq2(i));
            elseif model(i)==4
                Igd(i)=(Rga(i)*(ed1(i)-Vgd(i))+Xgq1(i)*(eq1(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd1(i)*Xgq1(i));
                Igq(i)=(-Xgd1(i)*(ed1(i)-Vgd(i))+Rga(i)*(eq1(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd1(i)*Xgq1(i));
            elseif model(i)==3
                Igd(i)=(Rga(i)*(-Vgd(i))+Xgq(i)*(eq1(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd1(i)*Xgq(i));
                Igq(i)=(-Xgd1(i)*(-Vgd(i))+Rga(i)*(eq1(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd1(i)*Xgq(i));
            elseif model(i)==2
                Igd(i)=(Rga(i)*(-Vgd(i))+Xgq(i)*(Ef(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd(i)*Xgq(i));
                Igq(i)=(-Xgd(i)*(-Vgd(i))+Rga(i)*(Ef(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd(i)*Xgq(i));
            end
        end
        
        Ig=(Igd.*sin(dg)+Igq.*cos(dg))+1j*(-Igd.*cos(dg)+Igq.*sin(dg));
    else
        Ig=zeros(0,1);
    end
    
    synS=synV.*conj(Ig);
    synP=real(synS);
    Qx=Q;
    Qx(syn(:,1))=Qx(syn(:,1))+imag(synS);
    
    pgTotal=zeros(nbus,1);
    pgTotal(pv(:,1))=pv(:,4);
    pgTotal(syn(:,1))=pgTotal(syn(:,1))+synP;
%     pgTotal(pq(:,1))=pgTotal(pq(:,1))-pq(:,4);
    busTypeNew=busType;
    busTypeNew(syn(:,1))=1;
    busTypeNew(sw(:,1))=2;
    
    f0=f;
    for island=1:nIslands
        busesInIsland=find(islands==island);
        agcCoeffInIsland=agcExt(busesInIsland,4);
        busTypeInIsland=busTypeNew(busesInIsland);
        pBusesInIsland=busesInIsland(busTypeInIsland>0);
        swBusesInIsland=busesInIsland(busTypeInIsland==2);
        
        f0(busesInIsland)=mean(f(busesInIsland));
        
        if sum(agcCoeffInIsland(busTypeInIsland>0))>0 % Effective AGC
            
        else % No Effective AGC in island
%             if isempty(swBusesInIsland)&&~isempty(pBusesInIsland)
%                 busTypeNew(pBusesInIsland(1))=2;
%             end
%             f0(busesInIsland)=0.0;
        end
        if ~isempty(swBusesInIsland)
            f0(busesInIsland)=0.0;
        end
    end
    
    agcExt(syn(:,1),3)=agcExt(syn(:,1),3)+syn(:,19);
    agcExt(syn(tg(:,1),1),2)=agcExt(syn(tg(:,1),1),2)+1./tg(:,4);
    
    TagcOri=1./agcExt(:,4);
    K=zeros(nbus,1);
    K(syn(:,1))=K(syn(:,1))+syn(:,19);
    K(syn(tg(:,1),1))=K(syn(tg(:,1),1))+1./tg(:,4);
    Tagc=zeros(nbus,1);
    Tagc(syn(:,1))=syn(:,18)./syn(:,19);
    Mtg=syn(tg(:,1),18);
    Dtg=syn(tg(:,1),18);
    Rtg=tg(:,4);
    T1tg=tg(:,8);
    T2tg=tg(:,7);
    Tagc(syn(tg(:,1),1))=2*Mtg.*T2tg.*Rtg./(Rtg.*Mtg+T1tg+Rtg.*Dtg.*T2tg-sqrt((Rtg.*Mtg-Rtg.*Dtg.*T2tg).*(Rtg.*Mtg-Rtg.*Dtg.*T2tg)+4*Rtg.*Mtg.*(T1tg-T2tg)));
    recTagc=(K.*TagcOri-sqrt(K.*K.*TagcOri.*TagcOri-4*K.*TagcOri.*Tagc))/2./TagcOri./Tagc;
    recTagc(isnan(recTagc))=0.0;
    agcExt(:,4)=recTagc;
    
    fdk=agcExt(:,2)+agcExt(:,3); 
    
    agcTag=zeros(nbus,1);
    agcTag(agc(:,1))=1;
    agcTag(fdk>0)=1;
    agcTag(agcExt(:,4)~=0)=1;
    agcNew=agcExt(agcTag==1,:);
    pvNewIdx=find(busTypeNew==1);
    swNewIdx=find(busTypeNew==2);
    pqDel=pq(busTypeNew(pq(:,1))~=0,:);
    pq(busTypeNew(pq(:,1))~=0,5)=0;
    Qx(pqDel(:,1))=Qx(pqDel(:,1))-pqDel(:,5);
    pvNew=[pvNewIdx,100.0*ones(size(pvNewIdx,1),1),bus(pvNewIdx,2),pgTotal(pvNewIdx),abs(V(pvNewIdx)),repmat([2.0,0.0,0,1],size(pvNewIdx,1),1)];
    swNew=[swNewIdx,100.0*ones(size(swNewIdx,1),1),bus(swNewIdx,2),abs(V(swNewIdx)),atan2(imag(V(swNewIdx)),real(V(swNewIdx)))*180/pi,repmat([99.00000 -99.00000 1.1 0.9  0 1 1 1],size(swNewIdx,1),1)];
    
    deltaFd=fdk.*f0;
    pvNew(:,4)=pvNew(:,4)+deltaFd(pvNew(:,1));
    deltaFd(pvNew(:,1))=0.0;
    
    pqVal=zeros(nbus,2);
    pqVal(pq(:,1),:)=pq(:,4:5);
    pqVal(:,1)=pqVal(:,1)-deltaFd;
    pqNew=[bus(:,1),100.0*ones(nbus,1),bus(:,2),pqVal,repmat([2.0,0.0,0,1],nbus,1)];
    pqNew=pqNew(pqNew(:,4)~=0|pqNew(:,5)~=0,:);
    
    SysDataNew=foldSysData(bus,swNew,pvNew,pqNew,shunt,line,ind,zip,zeros(0,26),zeros(0,14),zeros(0,8),agcNew,cac,cluster);
    [nState,idxs]...
        =getIndexDyn(SysDataNew);
    
    xx=zeros(nState,1);
    xx(idxs.vIdx)=V;
    xx(idxs.qIdx)=Qx;
    xx(idxs.sIdx)=s;
    deltaxx=zeros(nState,1);
    xNew=perpareInitialState(SysDataNew,xx,deltaxx);
    xNew(idxs.fIdx)=f0;
    dpg0=zeros(nbus,1);
    xNew(idxs.dpgIdx)=dpg0;
end
end