function diff=checkEquationBalanceSyn(SysData,SysPara,x,dx,dyn)
% Calculate equation imbalance 
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
% FUNCTION checkEquationBalanceSyn
%
% INPUT (WILL be modified in a future version)
%
% OUTPUT
%	diff - Equation imbalance vector
%

[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
[V,Q,s,d,w,eq1,eq2,ed1,ed2,psid,psiq,Pm,Ef,Vavrm,Vavrr,Vavrf,Vavrref,tgovg,tgovm,tgovmech,f,dpg,qplt,vg]=unfoldX(x,SysData);
[~,~,ds,dd,dw,deq1,deq2,ded1,ded2,dpsid,dpsiq,dPm,~,dVavrm,dVavrr,dVavrf,~,dtgovg,~,dtgovmech,df,ddpg,dqplt,dvg]=unfoldX(dx,SysData);
[pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,~,Ysh0,~,Vsp2,~,~,~,~,Tmech1,Varref1,Ef1,Pm1,Eq11]=unfoldSysPara(SysPara);

nbus=size(bus,1);
if ~isempty(zip)%zipMode=0
    Ysh0=Ysh0+accumarray(zip(:,1),(zip(:,5)+1j*zip(:,8)).*zip(:,12),[nbus,1]);
end
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

pVec=zeros(nbus,1);
qVec=zeros(nbus,1);

pVec(pv(:,1))=pVec(pv(:,1))+pv(:,4);
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

if ~isempty(syn)    
    nSyn=size(syn,1);
    synIdx=syn(:,1);
    Rs=syn(:,7);
    Xd=syn(:,8);
    Xq=syn(:,13);
    
    Efd=zeros(nSyn,1);
    Efq=Ef;
    cosd=cos(d);
    sind=sin(d);
    Cg=real(V(synIdx));
    Dg=imag(V(synIdx));
    Vd=sind.*Cg-cosd.*Dg;
    Vq=cosd.*Cg+sind.*Dg;
    Id=(Rs.*(Efd-Vd)+Xq.*(Efq-Vq))./(Rs.*Rs+Xq.*Xd);
    Iq=(-Xd.*(Efd-Vd)+Rs.*(Efq-Vq))./(Rs.*Rs+Xq.*Xd);
    IG=(sind.*Id+cosd.*Iq)+1j*(-cosd.*Id+sind.*Iq);
%     IG=(Ef.*(cos(d)+1j*sin(d))-V(synIdx))./(Rs+1j*Xd);
    SInjRHS_syn=-V.*conj(accumarray(synIdx,IG,[nbus,1]));
    SInjRHS=SInjRHS+SInjRHS_syn;
end

if ~isempty(ind)
    nInd=size(ind,1);
    indIdx=ind(:,1);
    R1=ind(:,7);
    X1=ind(:,8);
    Z1=ind(:,7)+1j*ind(:,8);
    Xm=ind(:,13);
    Ze=1j*ind(:,13);
    R2=ind(:,9);
    X2=ind(:,10);
    T0=ind(:,15)+ind(:,16)+ind(:,17);
    T1=-ind(:,16)-2*ind(:,17);
    T2=ind(:,17);
    H=ind(:,14);
    
    Zind=Z1+Ze.*(R2+1j*X2.*s)./(R2+(1j*X2+Ze).*s);
    ILind=V(indIdx)./Zind;
    SInjRHS_ind=V.*conj(accumarray(indIdx,ILind,[nbus,1]));
    SInjRHS=SInjRHS+SInjRHS_ind;
end

% qVec(busType~=0)=0;
diffSInj=SInjRHS-(pVec+1j*qVec);
diffSInj=diffSInj-1j*Q;
diffSInj(busType==2)=0;

diffsInd=zeros(0,1);
if ~isempty(ind)
    Vind=V(indIdx);
    VEind=Vind.*Ze./(Z1+Ze);
    IRs=VEind.*s./(s.*Z1.*Ze./(Z1+Ze)+R2+1j*X2.*s);
    ds(dyn==0)=0;
    sIndLHS=ds*2.*H.*s;
    diffsInd=real(sIndLHS-(T0+s.*(T1+s.*T2)).*s+IRs.*conj(IRs).*R2);
end

diffv=sum(Vsp2,2)-V.*conj(V);
diffv(busType==2)=0;

diff=[real(diffSInj);imag(diffSInj);diffv(busType~=0);diffsInd];
end