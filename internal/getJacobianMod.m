function [jac,pIdx,qIdx,vIdx,sIdx,busType]=getJacobianMod(nbus,nline,bus,sw,pv,pq,shunt,line,zip,ind,s,V,dyn,Ytr,Ysh)
% Generate jacobian matrix
%
% FUNCTION getJacobianMod
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

[~,Ytrx,Yshx,ytrfr,ytrto,yshfr,yshto]=getYMatrix(nbus,line);

if nargin<14
    Ytr=Ytrx;
end
if nargin<15
    Ysh=Yshx;
end

yBr=[-ytrfr;-ytrto];

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
busType(pv(:,1))=1;
busType(sw(:,1))=2;

yShunt=zeros(nbus,1);
yShunt(shunt(:,1))=shunt(:,5)+1j*shunt(:,6);
Ysh=Ysh+yShunt;
if ~isempty(zip)%zipMode=0
    Ysh=Ysh+accumarray(zip(:,1),(zip(:,5)+1j*zip(:,8)).*zip(:,12),[nbus,1]);
end
Y=Ytr+sparse(1:nbus,1:nbus,Ysh,nbus,nbus);

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
    
    Zind=Z1+Ze.*(R2+1j*X2.*s)./(R2+(1j*X2+Ze).*s);
    
    Y=Y+sparse(indIdx,indIdx,1./Zind,nbus,nbus);
end

isw=find(busType==2);
ipv=find(busType==1);
ipq=find(busType==0);
busTypePvPq=busType(busType~=2);
ipvInPvPq=find(busTypePvPq==1);
npq=size(ipq,1);
npv=size(ipv,1);

i_fr=line(:,1);
i_to=line(:,2);
yDiag=diag(Y);

C=real(V);
D=imag(V);
vMag=abs(V);

G=real(Y);
B=imag(Y);

gcbdx=G*C-B*D;
gcbd=sparse(1:nbus,1:nbus,gcbdx,nbus,nbus);
bcgdx=B*C+G*D;
bcgd=sparse(1:nbus,1:nbus,bcgdx,nbus,nbus);
diagc=sparse(1:nbus,1:nbus,C,nbus,nbus);
diagd=sparse(1:nbus,1:nbus,D,nbus,nbus);
cg=diagc*G;
cb=diagc*B;
dg=diagd*G;
db=diagd*B;

Hmat=gcbd+cg+db;  % P/C
Nmat=bcgd-cb+dg;  % P/D
Jmat=-bcgd-cb+dg; % Q/C
Lmat=gcbd-cg-db;  % Q/D

if ~isempty(zip)
    Czip=C(zip(:,1));
    Dzip=D(zip(:,1));
    vMagzip=vMag(zip(:,1));
    Hmat=Hmat+sparse(zip(:,1),zip(:,1),zip(:,6).*Czip./vMagzip,nbus,nbus);
    Nmat=Nmat+sparse(zip(:,1),zip(:,1),zip(:,6).*Dzip./vMagzip,nbus,nbus);
    Jmat=Jmat+sparse(zip(:,1),zip(:,1),zip(:,9).*Czip./vMagzip,nbus,nbus);
    Lmat=Lmat+sparse(zip(:,1),zip(:,1),zip(:,9).*Dzip./vMagzip,nbus,nbus);
end

% Hdiag=sparse(ones(2*nline+nbus,1),[i_to;i_fr;(1:nbus)'],[-(real(yBr).*C([i_to;i_fr])-imag(yBr).*D([i_to;i_fr]));-2*real(yDiag).*C],1,nbus);
% Ndiag=sparse(ones(2*nline+nbus,1),[i_to;i_fr;(1:nbus)'],[-(real(yBr).*D([i_to;i_fr])+imag(yBr).*C([i_to;i_fr]));-2*real(yDiag).*D],1,nbus);
% Jdiag=sparse(ones(2*nline+nbus,1),[i_to;i_fr;(1:nbus)'],[(real(yBr).*D([i_to;i_fr])+imag(yBr).*C([i_to;i_fr]));2*imag(yDiag).*C],1,nbus);
% Ldiag=sparse(ones(2*nline+nbus,1),[i_to;i_fr;(1:nbus)'],[-(real(yBr).*C([i_to;i_fr])-imag(yBr).*D([i_to;i_fr]));2*imag(yDiag).*D],1,nbus);
%
% Hmat=sparse([i_fr;i_to;(1:nbus)'],[i_to;i_fr;(1:nbus)'],[-(real(yBr).*C([i_fr;i_to])+imag(yBr).*D([i_fr;i_to]));Hdiag.'],nbus,nbus);
% Nmat=sparse([i_fr;i_to;(1:nbus)'],[i_to;i_fr;(1:nbus)'],[(imag(yBr).*C([i_fr;i_to])-real(yBr).*D([i_fr;i_to]));Ndiag.'],nbus,nbus);
% Jmat=sparse([i_fr;i_to;(1:nbus)'],[i_to;i_fr;(1:nbus)'],[(imag(yBr).*C([i_fr;i_to])-real(yBr).*D([i_fr;i_to]));Jdiag.'],nbus,nbus);
% Lmat=sparse([i_fr;i_to;(1:nbus)'],[i_to;i_fr;(1:nbus)'],[(real(yBr).*C([i_fr;i_to])+imag(yBr).*D([i_fr;i_to]));Ldiag.'],nbus,nbus);

% VdC=sparse(1:nbus,1:nbus,-2*C,nbus,nbus);
% VdD=sparse(1:nbus,1:nbus,-2*D,nbus,nbus);
C0Mx=sparse(1:npv,ipvInPvPq,-2*C(ipv),npv,npv+npq);
D0Mx=sparse(1:npv,ipvInPvPq,-2*D(ipv),npv,npv+npq);

jac=[Hmat(busType~=2,busType~=2),Nmat(busType~=2,busType~=2);...
    Jmat(busType==0,busType~=2),Lmat(busType==0,busType~=2);...
    C0Mx,D0Mx];

pIdx=1:(npv+npq);
qIdx=npv+npq+(1:npq);
vIdx=npv+npq+npq+(1:npv);
sIdx=[];

ind=ind(dyn==0,:);
nInd=size(ind,1);
if nInd>1
    indIdx=ind(:,1);
    s=s(dyn==0);
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
    
    Vind=V(indIdx);
    C1=real(Vind);
    C2=imag(Vind);
    
    y2=-(R2.*Xm.*Xm)./(s.*(R2./s - X2*1i + (Xm.*(R1 - X1*1i)*1i)./...
        (- R1 + X1*1i + Xm*1i)).*(R1 + X1*1i + Xm*1i).*(X2*1i + R2./s + (Xm.*(R1 + X1*1i)*1i)./(R1 + X1*1i + Xm*1i)).*(- R1 + X1*1i + Xm*1i)); % Te=V^2.*y2
    y2ds=(R2.*Xm.*Xm)./(s.*s.*(R2./s - X2*1i + (Xm.*(R1 - X1*1i)*1i)./(- R1 + X1*1i + Xm*1i)).*(R1 + X1*1i + Xm*1i).*(X2*1i + R2./s + (Xm.*(R1 + X1*1i)*1i)./(R1 + X1*1i + Xm*1i))...
        .*(- R1 + X1*1i + Xm*1i)) - (R2.*R2.*Xm.*Xm)./(s.*s.*s.*(R2./s - X2*1i + (Xm.*(R1 - X1*1i)*1i)./(- R1 + X1*1i + Xm*1i)).*(R1 + X1*1i + Xm*1i)...
        .*(X2*1i + R2./s + (Xm.*(R1 + X1*1i)*1i)./(R1 + X1*1i + Xm*1i)).^2.*(- R1 + X1*1i + Xm*1i)) - (R2.*Xm.*R2.*Xm)./(s.*s.*s.*(R2./s - X2*1i + (Xm.*(R1 - X1*1i)*1i)...
        ./(- R1 + X1*1i + Xm*1i)).^2.*(R1 + X1*1i + Xm*1i).*(X2*1i + R2./s + (Xm.*(R1 + X1*1i)*1i)./(R1 + X1*1i + Xm*1i)).*(- R1 + X1*1i + Xm*1i));
    dTds=(T1+2*T2.*s-abs(Vind).*abs(Vind).*y2ds)/2./H;
    dTdC1=-2*C1.*real(y2)/2./H;
    dTdC2=-2*C2.*real(y2)/2./H;
    
    yIL=(R2 + s.*(X2 + Xm)*1i)./(R1.*R2 - s.*(X1.*X2 + X1.*Xm + X2.*Xm) + R2.*(X1 + Xm)*1i + R1.*s.*(X2 + Xm)*1i);
    yILds=(X2*1i + Xm*1i)./(R1.*R2 - s.*(X1.*X2 + X1.*Xm + X2.*Xm) + R2.*(X1 + Xm)*1i + R1.*s.*(X2 + Xm)*1i) ...
        + (R2 + s.*(X2 + Xm)*1i).*(X1.*X2 + X1.*Xm + X2.*Xm - R1.*(X2 + Xm)*1i)./(R1.*R2 - s.*(X1.*X2 + X1.*Xm + X2.*Xm) + R2.*(X1 + Xm)*1i + R1.*s.*(X2 + Xm)*1i).^2;
    % VILdC1=-2*C1.*conj(yIL);
    % VILdC2=-2*C2.*conj(yIL);
    % Hmat=Hmat+sparse(indIdx,indIdx,real(VILdC1),nbus,nbus);
    % Nmat=Nmat+sparse(indIdx,indIdx,real(VILdC2),nbus,nbus);
    % Jmat=Jmat+sparse(indIdx,indIdx,imag(VILdC1),nbus,nbus);
    % Lmat=Lmat+sparse(indIdx,indIdx,imag(VILdC2),nbus,nbus);
    
    VILds=abs(Vind).*abs(Vind).*conj(yILds);
    PILdsmat=sparse(indIdx,1:nInd,-real(VILds),nbus,nInd);
    QILdsmat=sparse(indIdx,1:nInd,-imag(VILds),nbus,nInd);
    dTdC1mat=sparse(1:nInd,indIdx,dTdC1,nInd,nbus);
    dTdC2mat=sparse(1:nInd,indIdx,dTdC2,nInd,nbus);
    ind2ind=sparse(1:nInd,1:nInd,dTds,nInd,nInd);
    
    jac=[jac,[PILdsmat(busType~=2,:);QILdsmat(busType==0,:);sparse(npv,nInd)]...
        dTdC1mat(:,busType~=2),dTdC2mat(:,busType~=2),ind2ind];
    sIdx=vIdx(end)+(1:nInd);
end

end