function x0=perpareInitialState(SysData,x,dx)
% Generate initial states of a system
%
% FUNCTION perpareInitialState
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
%   x - Partially filled states
%	dx - time derivative of x
%	
% OUTPUT
%	x0 - Fully filled state
%
[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
[V,Q,s,d,w,eq1,eq2,ed1,ed2,psid,psiq,Pm,Ef,Vavrm,Vavrr,Vavrf,Vavrref,tgovg,tgovm,tgovmech,f,dpg,qplt,vg]=unfoldX(x,SysData);
[~,~,ds,dd,dw,deq1,deq2,ded1,ded2,dpsid,dpsiq,dPm,~,dVavrm,dVavrr,dVavrf,~,dtgovg,~,dtgovmech,df,ddpg,dqplt,dvg]=unfoldX(dx,SysData);
[nState,idxs]...
        =getIndexDyn(SysData);
%
nbus=size(bus,1);
nInd=size(ind,1);
nSyn=size(syn,1);
nExc=size(exc,1);
nTg=size(tg,1);

nCac=size(idxs.qpltIdx,1);
nCluster=size(idxs.vgIdx,1);

x0=zeros(nbus+nInd+10*nSyn+4*nExc+3*nTg+2*nbus+nCac+nCluster,1);

x0(idxs.vIdx)=V;
x0(idxs.qIdx)=Q;
x0(idxs.sIdx)=s;
x0(idxs.deltaIdx)=d;
x0(idxs.omegaIdx)=w;

agcExt=zeros(nbus,4);
agcExt(1:nbus,1)=1:nbus;
if ~isempty(agc)
    agcExt(agc(:,1),:)=agc;
end
fdk=agcExt(:,2)+agcExt(:,3);      

if ~isempty(syn)
    nSyn=size(syn,1);
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
    Mxg=syn(:,18);
    Dxg=syn(:,19);
    TgAA=syn(:,24);
    gammad=Tgd2./Tgd1.*Xgd2./Xgd1.*(Xgd-Xgd1);
    gammaq=Tgq2./Tgq1.*Xgq2./Xgq1.*(Xgq-Xgq1);
    
    Cg=real(V(synIdx));
    Dg=imag(V(synIdx));    
    sind=sin(d);
    cosd=cos(d);    
    Efq=Ef;
    Efd=zeros(nSyn,1);
    
    Vgd=sind.*Cg-cosd.*Dg;
    Vgq=cosd.*Cg+sind.*Dg;
    Igd=(Rga.*(Efd-Vgd)+Xgq.*(Efq-Vgq))./(Rga.*Rga+Xgq.*Xgd);
    Igq=(-Xgd.*(Efd-Vgd)+Rga.*(Efq-Vgq))./(Rga.*Rga+Xgq.*Xgd);
    for i=1:nSyn
        if model(i)==8            
            x0(idxs.psiqIdx(i))=-(Rga(i)*Igd(i)+Vgd(i));
            x0(idxs.psidIdx(i))=(Rga(i)*Igq(i)+Vgq(i));
            x0(idxs.eq1Idx(i))=(-(Xgd(i)-Xgd1(i)-gammad(i))*Igd(i)+(1-TgAA(i)/Tgd1(i))*Ef(i));
            x0(idxs.ed1Idx(i))=((Xgq(i)-Xgq1(i)-gammaq(i))*Igq(i));
            x0(idxs.eq2Idx(i))=(x0(idxs.eq1Idx(i))-(Xgd1(i)-Xgd2(i)+gammad(i))*Igd(i)+TgAA(i)/Tgd1(i)*Ef(i));
            x0(idxs.ed2Idx(i))=(x0(idxs.ed1Idx(i))+(Xgq1(i)-Xgq2(i)+gammaq(i))*Igq(i));
        elseif model(i)==6           
            x0(idxs.psiqIdx(i))=-(Rga(i)*Igd(i)+Vgd(i));
            x0(idxs.psidIdx(i))=(Rga(i)*Igq(i)+Vgq(i));
            x0(idxs.eq1Idx(i))=(-(Xgd(i)-Xgd1(i)-gammad(i))*Igd(i)+(1-TgAA(i)/Tgd1(i))*Ef(i));
            x0(idxs.ed1Idx(i))=((Xgq(i)-Xgq1(i)-gammaq(i))*Igq(i));
            x0(idxs.eq2Idx(i))=(x0(idxs.eq1Idx(i))-(Xgd1(i)-Xgd2(i)+gammad(i))*Igd(i)+TgAA(i)/Tgd1(i)*Ef(i));
            x0(idxs.ed2Idx(i))=(x0(idxs.ed1Idx(i))+(Xgq1(i)-Xgq2(i)+gammaq(i))*Igq(i));
        elseif model(i)==5   
            x0(idxs.psiqIdx(i))=-(Rga(i)*Igd(i)+Vgd(i));
            x0(idxs.psidIdx(i))=(Rga(i)*Igq(i)+Vgq(i));
            x0(idxs.eq1Idx(i))=(-(Xgd(i)-Xgd1(i)-gammad(i))*Igd(i)+(1-TgAA(i)/Tgd1(i))*Ef(i));
            x0(idxs.eq2Idx(i))=(x0(idxs.eq1Idx(i))-(Xgd1(i)-Xgd2(i)+gammad(i))*Igd(i)+TgAA(i)/Tgd1(i)*Ef(i));
            x0(idxs.ed2Idx(i))=((Xgq(i)-Xgq2(i))*Igq(i));
        elseif model(i)==4  
            x0(idxs.psiqIdx(i))=-(Rga(i)*Igd(i)+Vgd(i));
            x0(idxs.psidIdx(i))=(Rga(i)*Igq(i)+Vgq(i));
            x0(idxs.eq1Idx(i))=(-(Xgd(i)-Xgd1(i))*Igd(i)+Ef(i));
            x0(idxs.ed1Idx(i))=((Xgq(i)-Xgq1(i))*Igq(i));          
        elseif model(i)==3
            x0(idxs.psiqIdx(i))=-(Rga(i)*Igd(i)+Vgd(i));
            x0(idxs.psidIdx(i))=(Rga(i)*Igq(i)+Vgq(i));
%             x0(eq1Idx(i))=(-(Xgd(i)-Xgd1(i))*Igd(i)+Ef(i)); 
            x0(idxs.eq1Idx(i))=Vgq(i)+Rga(i)*Igq(i)+Xgd1(i)*Igd(i);    
        elseif model(i)==2
            x0(idxs.psiqIdx(i))=-(Rga(i)*Igd(i)+Vgd(i));
            x0(idxs.psidIdx(i))=(Rga(i)*Igq(i)+Vgq(i));
        end
    end
    Pg=(Vgq+Igq.*Rga).*Igq+(Vgd+Igd.*Rga).*Igd+Dxg.*w+Mxg.*dw;
    x0(idxs.pgIdx)=Pg;
    x0(idxs.efIdx)=Ef;
end

if ~isempty(exc)
    excIdx=exc(:,1);
    muavr0=exc(:,5);
    Tavr1=exc(:,7);
    Tavr2=exc(:,6);
    vavrf0=exc(:,8);
    Vavr0=exc(:,9);
    Tavre=exc(:,10);
    Tavrr=exc(:,11);
    
    Vmag=abs(V);
    Vmagavr=Vmag(synIdx(excIdx));
    x0(idxs.vavrmIdx)=Vmagavr;
    x0(idxs.vavrrefIdx)=(Ef(excIdx).*Vavr0./Vmagavr-vavrf0)./muavr0+Vmagavr;
    x0(idxs.vavrrIdx)=muavr0.*(1-Tavr1./Tavr2).*(x0(idxs.vavrrefIdx)-x0(idxs.vavrmIdx));
    x0(idxs.vavrfIdx)=Ef(excIdx);
end

if ~isempty(tg)
    tgIdx=tg(:,1);
    
    wtgref=tg(:,3);
    Rtg=tg(:,4);
    Ttgmax=tg(:,5);
    Ttgmin=tg(:,6);
    Ttg2=tg(:,7);
    Ttg1=tg(:,8);
    
    x0(idxs.tgovgIdx)=(1-Ttg1./Ttg2).*(wtgref-x0(idxs.omegaIdx(tgIdx)))./Rtg-dtgovg.*Ttg2;
    x0(idxs.tgovmIdx)=Pg(tgIdx);
    x0(idxs.tmechIdx)=Pg(tgIdx)-Ttg1./Ttg2.*(wtgref-x0(idxs.omegaIdx(tgIdx)))./Rtg-x0(idxs.tgovgIdx);    
end

% the initial frequency equals nonminal freq.
x0(idxs.fIdx)=0;
x0(idxs.dpgIdx)=0;

%TODO: temporary
x0(idxs.qpltIdx)=0;
x0(idxs.vgIdx)=0;

end