function [V,W,Q,f]=hemMachinePFSalientcontinueAlgebraic(SimData,SysData,SysPara,x0)
% HE approach for solving algebraic equations
%
% FUNCTION hemMachinePFSalientcontinueAlgebraic
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
%
% OUTPUT
%   V - HE coefficients of voltage
%   W - HE coefficients of reciprocal of voltage
%   Q - HE coefficients of reactive power
%   f - HE coefficients of frequency
%

[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
nbus=size(bus,1);
nline=size(line,1);
nSyn=size(syn,1);
[V0,Q0,s0,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,tgovmech0,f0,dpg0,qplt0,vg0]=unfoldX(x0,SysData);
[~,~,~,nlvl,~,~,~,~,~]=unfoldSimData(SimData);
[pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,MatGV0,MatGV1,MatGRhs0,MatGRhs1]=unfoldSysPara(SysPara);
%
[nIslands,islands,refs]=searchIslands(bus(:,1),line(:,1:2));

Pls=zeros(nbus,1);Pls(pq(:,1))=pqIncr(:,1);if ~isempty(pv);Pls(pv(:,1))=Pls(pv(:,1))-pvIncr;end
Qls=zeros(nbus,1);Qls(pq(:,1))=pqIncr(:,2);

%     [Y,Ytr,Ysh,ytrfr,ytrto,yshfr,yshto]=getYMatrix(nbus,line);
Ytr=Ytr0;

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
ipv=find(busType==1);
ipq=find(busType==0);
npq=size(ipq,1);
npv=size(ipv,1);

%     yShunt=zeros(nbus,1);
%     yShunt(shunt(:,1))=shunt(:,5)+1j*shunt(:,6);
%     if ~isempty(zip)%zipMode=0
%         yShunt=yShunt+accumarray(zip(:,1),(zip(:,5)+1j*zip(:,8)).*zip(:,12),[nbus,1]);
%     end
%     Ysh=Ysh+yShunt;
%     Y=Y+sparse(1:nbus,1:nbus,yShunt,nbus,nbus);
yShunt=zeros(nbus,1);
yShunt(shunt(:,1))=shunt(:,5)+1j*shunt(:,6);
if ~isempty(zip)%zipMode=0
    Ysh0=Ysh0+accumarray(zip(:,1),Rzip0.*(zip(:,5)+1j*zip(:,8)).*zip(:,12),[nbus,1]);
    Ysh1=Ysh1+accumarray(zip(:,1),Rzip1.*(zip(:,5)+1j*zip(:,8)).*zip(:,12),[nbus,1]);
end
Ysh0=Ysh0+yShunt;

Y=Ytr+sparse(1:nbus,1:nbus,Ysh0,nbus,nbus);

pVec=zeros(nbus,1);
qVec=zeros(nbus,1);
%     vSp=zeros(nbus,1);

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

f=zeros(nbus,nlvl+1);
f(:,1)=f0;

if ~isempty(agc)
    agcExt=zeros(nbus,size(agc,2));
    agcExt(agc(:,1),:)=agc;
    fdk=agcExt(:,2)+agcExt(:,3); %1/R+D
else
    fdk=zeros(nbus,1);
end

pVec(pv(:,1))=pVec(pv(:,1))+pv(:,4);
pVec(pq(:,1))=pVec(pq(:,1))-pq(:,4);
qVec(pq(:,1))=qVec(pq(:,1))-pq(:,5);
if ~isempty(zip)%zipMode=0, account for the PQ components in ZIP loads
    pVec=pVec-accumarray(zip(:,1),Rzip0.*zip(:,7).*zip(:,12),[nbus,1]);
    qVec=qVec-accumarray(zip(:,1),Rzip0.*zip(:,10).*zip(:,12),[nbus,1]);
end
% qVec(ipv)=Q0(ipv);
%     vSp(ipv)=pv(:,5);

V=zeros(nbus,nlvl+1);
V(:,1)=V0;
W=zeros(nbus,nlvl+1);
W(:,1)=1./V0;
P=zeros(nbus,nlvl+1);
P(:,1)=pVec;
%     P(isw,2:end)=0;
Q=zeros(nbus,nlvl+1);
Qxtra=zeros(size(Q));
Q(:,1)=Q0;
Qxtra(:,1)=qVec;
P(:,2:(size(Pls,2)+1))=-Pls;
Qxtra(:,2:(size(Qls,2)+1))=-Qls;
if ~isempty(zip)
    P(:,2)=P(:,2)-accumarray(zip(:,1),Rzip1.*zip(:,7).*zip(:,12),[nbus,1]);
    Qxtra(:,2)=Qxtra(:,2)-accumarray(zip(:,1),Rzip1.*zip(:,10).*zip(:,12),[nbus,1]);
end

C0=real(V(:,1));
D0=imag(V(:,1));
E0=real(W(:,1));
F0=imag(W(:,1));

C0M=sparse(1:nbus,1:nbus,C0,nbus,nbus);
D0M=sparse(1:nbus,1:nbus,D0,nbus,nbus);
E0M=sparse(1:nbus,1:nbus,E0,nbus,nbus);
F0M=sparse(1:nbus,1:nbus,F0,nbus,nbus);
P0M=sparse(1:nbus,1:nbus,P(:,1),nbus,nbus);
Q0M=sparse(1:nbus,1:nbus,Q(:,1)+Qxtra(:,1),nbus,nbus);

G=real(Y);
B=imag(Y);

if ~isempty(zip)
    nZip=size(zip,1);
    zipIdx=zip(:,1);
    IiL=zeros(nZip,nlvl+1);
    BiL=zeros(nZip,nlvl+1);
    
    Bi0=abs(V0(zipIdx));
    JI=zip(:,6);
    KI=-zip(:,9);
    
    Ii0L=Rzip0.*(JI+1j*KI).*V0(zipIdx)./Bi0;
    Ji0L=real(Ii0L);
    Ki0L=imag(Ii0L);
    
    IiL(:,1)=Ii0L;
    BiL(:,1)=Bi0;
    
    Ci0=real(V0(zipIdx));
    Di0=imag(V0(zipIdx));
    
    LHS_MatZip=[Rzip0.*JI./Bi0-Ci0.*Ji0L./Bi0./Bi0,-Rzip0.*KI./Bi0-Di0.*Ji0L./Bi0./Bi0,...
        Rzip0.*KI./Bi0-Ci0.*Ki0L./Bi0./Bi0,Rzip0.*JI./Bi0-Di0.*Ki0L./Bi0./Bi0];
    Mat_BZip=[Ci0./Bi0,Di0./Bi0];
else
    IiL=[];
end

if ~isempty(syn)
    synIdx=syn(:,1);
    MatGCD=-[sparse(synIdx,synIdx,MatGV0(:,1),nbus,nbus),sparse(synIdx,synIdx,MatGV0(:,2),nbus,nbus);...
        sparse(synIdx,synIdx,MatGV0(:,3),nbus,nbus),sparse(synIdx,synIdx,MatGV0(:,4),nbus,nbus)];
end

FreqReal=sparse(1:nbus,1:nbus,-freqKeptTag.*fdk.*E0,nbus,nbus);
FreqImag=sparse(1:nbus,1:nbus,-freqKeptTag.*fdk.*F0,nbus,nbus);
Freq2freq=sparse([1:nbus,1:nbus],[1:nbus,frefs(islands)'],[ones(1,nbus),-ones(1,nbus)],nbus,nbus);

Y11=-G;Y12=B;Y21=-B;Y22=-G;
YEF11=P0M+sparse(1:nbus,1:nbus,freqKeptTag.*(-fdk.*f0+dpg0),nbus,nbus);YEF12=-Q0M;YEF21=-Q0M;YEF22=-P0M-sparse(1:nbus,1:nbus,freqKeptTag.*(-fdk.*f0+dpg0),nbus,nbus);

if ~isempty(zip)
    Y11=Y11-sparse(1:nbus,1:nbus,accumarray(zipIdx,LHS_MatZip(:,1),[nbus,1]),nbus,nbus);
    Y12=Y12-sparse(1:nbus,1:nbus,accumarray(zipIdx,LHS_MatZip(:,2),[nbus,1]),nbus,nbus);
    Y21=Y21-sparse(1:nbus,1:nbus,accumarray(zipIdx,LHS_MatZip(:,3),[nbus,1]),nbus,nbus);
    Y22=Y22-sparse(1:nbus,1:nbus,accumarray(zipIdx,LHS_MatZip(:,4),[nbus,1]),nbus,nbus);
end
YLHS=[Y11,Y12;Y21,Y22];

if ~isempty(syn)
    YLHS=YLHS+MatGCD;
end
% 
% idxNonSw=find(busType~=2);
% idxStackMat=[idxNonSw;idxNonSw+nbus];
% 
% LHS_mat=[YLHS(idxStackMat,idxStackMat),...
%     [YEF11(busType~=2,busType~=2),YEF12(busType~=2,busType~=2),-F0M(busType~=2,ipv);...
%     YEF21(busType~=2,busType~=2),YEF22(busType~=2,busType~=2),-E0M(busType~=2,ipv)];...
%     C0M(ipv,busType~=2),D0M(ipv,busType~=2),sparse(npv,2*npq+3*npv);...
%     E0M(busType~=2,busType~=2),-F0M(busType~=2,busType~=2),C0M(busType~=2,busType~=2),-D0M(busType~=2,busType~=2),sparse(npq+npv,npv);...
%     F0M(busType~=2,busType~=2),E0M(busType~=2,busType~=2),D0M(busType~=2,busType~=2),C0M(busType~=2,busType~=2),sparse(npq+npv,npv);];

idxNonSw=find(busType~=2);
idxNonSwxD=find(fswTagxD==0);
idxNonSwD=find(busType~=2&fswTagxD==1);

LHS_mat=[YLHS([idxNonSw;idxNonSw+nbus],[idxNonSw;idxNonSw+nbus]),...
    [YEF11(idxNonSw,idxNonSw),YEF12(idxNonSw,idxNonSw),-F0M(idxNonSw,ipv),FreqReal(idxNonSw,freqKeptTag==1);...
    YEF21(idxNonSw,idxNonSw),YEF22(idxNonSw,idxNonSw),-E0M(idxNonSw,ipv),-FreqImag(idxNonSw,freqKeptTag==1)];...
    C0M(ipv,idxNonSw),D0M(ipv,idxNonSw),sparse(npv,2*npq+3*npv+nFreqKept);...
    E0M(idxNonSw,idxNonSw),-F0M(idxNonSw,idxNonSw),C0M(idxNonSw,idxNonSw),-D0M(idxNonSw,idxNonSw),sparse(npq+npv,npv+nFreqKept);...
    F0M(idxNonSw,idxNonSw),E0M(idxNonSw,idxNonSw),D0M(idxNonSw,idxNonSw),C0M(idxNonSw,idxNonSw),sparse(npq+npv,npv+nFreqKept);...
    sparse(sum(freqKeptTagxRef),size(idxNonSw,1)+size(idxNonSw,1)+2*npq+3*npv),Freq2freq(freqKeptTagxRef==1,freqKeptTag==1);...
    sparse(size(idxNonSwD,1),size(idxNonSw,1)),sparse(1:size(idxNonSwD,1),idxNonSwD,ones(size(idxNonSwD,1),1),size(idxNonSwD,1),size(idxNonSw,1)),sparse(size(idxNonSwD,1),2*npq+3*npv+nFreqKept)];

%
% if nbus<=500
%     [L_LHS_mat,U_LHS_mat,p_LHS_mat]=lu(LHS_mat,'vector');
% end
% 
% p_amd = colamd (LHS_mat) ;
% MxI = speye (size(LHS_mat)) ;
% MxQ = MxI (:, p_amd) ;
% [MxL,MxU,MxP] = lu (LHS_mat*MxQ) ;

useLU=isfield(SysPara,'iden')&&isfield(SysPara,'p_amd');

if useLU
    if isempty(SysPara.p_amd);
        p_amd = colamd (LHS_mat) ;
        save([SysPara.iden,'.mat'],'p_amd');
    else
        p_amd=SysPara.p_amd;
    end
    MxI = speye (size(LHS_mat)) ;
    MxQ = MxI (:, p_amd) ;
    [MxL,MxU,MxP] = lu (LHS_mat*MxQ) ;
end

for i=1:nlvl
    seq2=getseq(i,2);
    seq2p=getseq(i+1,2);
    seq3=getseq(i,3);
    idxSeq2=sum(seq2==i,2);
    seq2R=seq2(idxSeq2==0,:);
    
    RHSILr=zeros(nbus,1);
    RHSILi=zeros(nbus,1);
    
    RHSIiLr=zeros(nbus,1);
    RHSIiLi=zeros(nbus,1);
    if ~isempty(zip)
        RHS_BZip=(real(sum(V(zipIdx,seq2R(:,1)+1).*conj(V(zipIdx,seq2R(:,2)+1)),2))-sum(BiL(:,seq2R(:,1)+1).*BiL(:,seq2R(:,2)+1),2))./Bi0/2;
        RHZ_BIConv=sum(IiL(:,seq2R(:,1)+1).*BiL(:,seq2R(:,2)+1),2);
        RHSILr_full=Rzip1.*(JI.*real(V(zipIdx,i))-KI.*imag(V(zipIdx,i)))./Bi0-real(RHZ_BIConv)./Bi0-Ji0L.*RHS_BZip./Bi0;
        RHSILi_full=Rzip1.*(KI.*real(V(zipIdx,i))+JI.*imag(V(zipIdx,i)))./Bi0-imag(RHZ_BIConv)./Bi0-Ki0L.*RHS_BZip./Bi0;
        RHSIiLr=accumarray(zipIdx,RHSILr_full,[nbus,1]);
        RHSIiLi=accumarray(zipIdx,RHSILi_full,[nbus,1]);
    end
    
    RHSIGr=zeros(nbus,1);
    RHSIGi=zeros(nbus,1);
    if ~isempty(syn)
        RHSIGr=-accumarray(synIdx,MatGV1(:,1).*real(V(synIdx,i))+MatGV1(:,2).*imag(V(synIdx,i)),[nbus,1]);
        RHSIGi=-accumarray(synIdx,MatGV1(:,3).*real(V(synIdx,i))+MatGV1(:,4).*imag(V(synIdx,i)),[nbus,1]);
        if i==1
            RHSIGr=RHSIGr+accumarray(synIdx,MatGRhs1(:,1),[nbus,1]);
            RHSIGi=RHSIGi+accumarray(synIdx,MatGRhs1(:,2),[nbus,1]);
        end
    end
    
    % HEM Body
    RHS1=sum((-P(:,seq2(:,1)+1)+1j*(Q(:,seq2(:,1)+1)+Qxtra(:,seq2(:,1)+1))).*conj(W(:,seq2(:,2)+1)),2)+Ysh1.*V(:,i)+Ytr1*V(:,i);
    RHS2=-0.5*real(sum(V(:,seq2R(:,1)+1).*conj(V(:,seq2R(:,2)+1)),2));
    RHS3=sum(-W(:,seq2R(:,1)+1).*V(:,seq2R(:,2)+1),2);
    
    if i==1
        RHS2=RHS2+0.5*VspSq2(:,2);
    end
    
    compactRHS1=RHS1(busType~=2);
    if ~isempty(isw);compactRHS1=compactRHS1+Y(busType~=2,isw)*V(isw,i+1);end
    RHS=[real(compactRHS1)+RHSILr(busType~=2)+RHSIiLr(busType~=2)-RHSIGr(busType~=2);...
        imag(compactRHS1)+RHSILi(busType~=2)+RHSIiLi(busType~=2)-RHSIGi(busType~=2);...
        RHS2(ipv);...
        real(RHS3(busType~=2));...
        imag(RHS3(busType~=2));...
        zeros(sum(freqKeptTagxRef),1);...
        zeros(size(idxNonSwD,1),1)];
    
    if useLU
        x = MxQ * (MxU \ (MxL \ (MxP * RHS))) ;
    else
        x=LHS_mat\RHS;
    end
%     
%     V(busType~=2,i+1)=x(1:(npq+npv))+1j*x(((npq+npv)+1):(2*(npq+npv)));
%     W(busType~=2,i+1)=x((2*(npq+npv)+1):(3*(npq+npv)))+1j*x((3*(npq+npv)+1):(4*(npq+npv)));
%     Q(ipv,i+1)=x((4*(npq+npv)+1):end);
%     
    xC=real(V(:,i+1));
    xD=imag(V(:,i+1));
    xC(idxNonSw)=x(1:(npq+npv));
    xD(idxNonSw)=x(((npq+npv)+1):(2*(npq+npv)));
    V(:,i+1)=xC+1j*xD;
    W(busType~=2,i+1)=x((2*(npq+npv)+1):(3*(npq+npv)))+...
        1j*x((3*(npq+npv)+1):(4*(npq+npv)));
    Q(ipv,i+1)=x((4*(npq+npv)+1):(4*(npq+npv)+npv));
    f(freqKeptTag==1,i+1)=x((4*(npq+npv)+npv+1):end);
    
    if ~isempty(zip)
        IiL(:,i+1)=(LHS_MatZip(:,1)+1j*LHS_MatZip(:,3)).*real(V(zipIdx,i+1))+(LHS_MatZip(:,2)+1j*LHS_MatZip(:,4)).*imag(V(zipIdx,i+1))+(RHSILr_full+1j*RHSILi_full);
        BiL(:,i+1)=Mat_BZip(:,1).*real(V(zipIdx,i+1))+Mat_BZip(:,2).*imag(V(zipIdx,i+1))+RHS_BZip;
    end
    
end
end