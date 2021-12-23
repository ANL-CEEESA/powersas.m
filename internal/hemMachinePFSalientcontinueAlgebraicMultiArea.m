function [Vm,Wm,Qm]=hemMachinePFSalientcontinueAlgebraicMultiArea(SimData,SysData,SysPara,SysDataList,SysParaList,links,grossMaps,sysMaps,x0)
% Compute HE coefficients for solving algebraic equations in multiple areas (parallel computation enabled)
%
% FUNCTION hemMachinePFSalientcontinueAlgebraicMultiArea
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
%	SysDataList - The list of sub system data
%	SysParaList - The list of sub system parameters
%	links - the links connecting the sub systems
%	grossMaps - Mapping from agglomerated system to subsystem
%	sysMaps - Mapping from subsystem to agglomerated system
%   x0 - Initial system state
%	
% OUTPUT
%	Vm, Wm, Qm - HE oefficients of algebraic variables
%
%
nSys=length(SysDataList);
auxList=cell(nSys,1);
[~,~,~,nlvl,~,~,~,~,~]=unfoldSimData(SimData);
[busMain,swMain,pvMain,pqMain,shuntMain,lineMain,indMain,zipMain,synMain,excMain,tgMain,agcMain,cacMain,clusterMain]=unfoldSysData(SysData);
[~,~,~,~,~,~,Rzip0Main,~,~,~,~,~,~,~,~,~,~]=unfoldSysPara(SysPara);
[V0,Q0]=unfoldX(x0,SysData);
nbusMain=size(busMain,1);
Vm=zeros(nbusMain,nlvl+1);
Vm(:,1)=V0;
Wm=zeros(nbusMain,nlvl+1);
Wm(:,1)=1./V0;
busTypeMain=zeros(nbusMain,1);
if isempty(pvMain)
    pvMain=zeros(0,6);
end
if isempty(pqMain)
    pqMain=zeros(0,6);
end
if isempty(shuntMain)
    shuntMain=zeros(0,7);
end
if isempty(swMain)
    swMain=zeros(0,13);
end
busTypeMain(pvMain(:,1))=1;
busTypeMain(swMain(:,1))=2;

iswMain=find(busTypeMain==2);
ipvMain=find(busTypeMain==1);
ipqMain=find(busTypeMain==0);
npqMain=size(ipqMain,1);
npvMain=size(ipvMain,1);
qVecMain=zeros(nbusMain,1);
%     vSp=zeros(nbus,1);
qVecMain(pqMain(:,1))=qVecMain(pqMain(:,1))-pqMain(:,5);
if ~isempty(zipMain)%zipMode=0, account for the PQ components in ZIP loads
    qVecMain=qVecMain-accumarray(zipMain(:,1),Rzip0Main.*zipMain(:,10).*zipMain(:,12),[nbusMain,1]);
end
qVecMain(ipvMain)=Q0(ipvMain);
Qm=zeros(nbusMain,nlvl+1);
Qm(:,1)=qVecMain;

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end

for iSys=2:nSys
    [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg]=unfoldSysData(SysDataList{iSys});
    busInX=[sysMaps{iSys}.busMap';sysMaps{iSys}.extBusMap(:,2)];
    V0a=V0(busInX);
    Q0a=Q0(busInX);
    
    nbus=size(bus,1);
    nline=size(line,1);
    [pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,MatGV0,MatGV1,MatGRhs0,MatGRhs1]=unfoldSysPara(SysParaList{iSys});
    %
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
    if ~isempty(zip)%zipMode=0
        Ysh0=Ysh0+accumarray(zip(:,1),Rzip0.*(zip(:,5)+1j*zip(:,8)).*zip(:,12),[nbus,1]);
        Ysh1=Ysh1+accumarray(zip(:,1),Rzip1.*(zip(:,5)+1j*zip(:,8)).*zip(:,12),[nbus,1]);
    end
    
    Y=Ytr+sparse(1:nbus,1:nbus,Ysh0,nbus,nbus);
    
    pVec=zeros(nbus,1);
    qVec=zeros(nbus,1);
    %     vSp=zeros(nbus,1);
    
    pVec(pv(:,1))=pVec(pv(:,1))+pv(:,4);
    pVec(pq(:,1))=pVec(pq(:,1))-pq(:,4);
    qVec(pq(:,1))=qVec(pq(:,1))-pq(:,5);
    if ~isempty(zip)%zipMode=0, account for the PQ components in ZIP loads
        pVec=pVec-accumarray(zip(:,1),Rzip0.*zip(:,7).*zip(:,12),[nbus,1]);
        qVec=qVec-accumarray(zip(:,1),Rzip0.*zip(:,10).*zip(:,12),[nbus,1]);
    end
    qVec(ipv)=Q0a(ipv);
    %     vSp(ipv)=pv(:,5);
    
    V=zeros(nbus,nlvl+1);
    V(:,1)=V0a;
    W=zeros(nbus,nlvl+1);
    W(:,1)=1./V0a;
    P=zeros(nbus,nlvl+1);
    P(:,1)=pVec;
    %     P(isw,2:end)=0;
    Q=zeros(nbus,nlvl+1);
    Q(:,1)=qVec;
    P(:,2:(size(Pls,2)+1))=-Pls;
    Q(:,2:(size(Qls,2)+1))=-Qls;
    if ~isempty(zip)
        P(:,2)=P(:,2)-accumarray(zip(:,1),Rzip1.*zip(:,7).*zip(:,12),[nbus,1]);
        Q(:,2)=Q(:,2)-accumarray(zip(:,1),Rzip1.*zip(:,10).*zip(:,12),[nbus,1]);
    end
    Q(busType~=0,2:end)=0;
    
    C0=real(V(:,1));
    D0=imag(V(:,1));
    E0=real(W(:,1));
    F0=imag(W(:,1));
    
    C0M=sparse(1:nbus,1:nbus,C0,nbus,nbus);
    D0M=sparse(1:nbus,1:nbus,D0,nbus,nbus);
    E0M=sparse(1:nbus,1:nbus,E0,nbus,nbus);
    F0M=sparse(1:nbus,1:nbus,F0,nbus,nbus);
    
    P0x=P(:,1);Q0x=Q(:,1);
    P0M=sparse(1:nbus,1:nbus,P0x,nbus,nbus);
    Q0M=sparse(1:nbus,1:nbus,Q0x,nbus,nbus);
    
    G=real(Y);
    B=imag(Y);
    
    nZip=size(zip,1);
    zipIdx=zip(:,1);
    if ~isempty(zip)
        IiL=zeros(nZip,nlvl+1);
        BiL=zeros(nZip,nlvl+1);
        
        Bi0=abs(V0a(zipIdx));
        JI=zip(:,6);
        KI=-zip(:,9);
        
        Ii0L=Rzip0.*(JI+1j*KI).*V0a(zipIdx)./Bi0;
        Ji0L=real(Ii0L);
        Ki0L=imag(Ii0L);
        
        IiL(:,1)=Ii0L;
        BiL(:,1)=Bi0;
        
        Ci0=real(V0a(zipIdx));
        Di0=imag(V0a(zipIdx));
        
        LHS_MatZip_p=[Rzip0.*JI./Bi0-Ci0.*Ji0L./Bi0./Bi0,-Rzip0.*KI./Bi0-Di0.*Ji0L./Bi0./Bi0,...
            Rzip0.*KI./Bi0-Ci0.*Ki0L./Bi0./Bi0,Rzip0.*JI./Bi0-Di0.*Ki0L./Bi0./Bi0];
        Mat_BZip=[Ci0./Bi0,Di0./Bi0];
        
        auxList{iSys}.LHS_MatZip=LHS_MatZip_p;
        auxList{iSys}.Mat_BZip=Mat_BZip;
    else
        IiL=zeros(nZip,nlvl+1);
        BiL=zeros(nZip,nlvl+1);
    end
    
    synIdx=syn(:,1);
    if ~isempty(syn)
        MatGCD_p=-[sparse(synIdx,synIdx,MatGV0(:,1),nbus,nbus),sparse(synIdx,synIdx,MatGV0(:,2),nbus,nbus);...
            sparse(synIdx,synIdx,MatGV0(:,3),nbus,nbus),sparse(synIdx,synIdx,MatGV0(:,4),nbus,nbus)];
    end
    
    Y11=-G;Y12=B;Y21=-B;Y22=-G;
    YEF11=P0M;YEF12=-Q0M;YEF21=-Q0M;YEF22=-P0M;
    
    if ~isempty(zip)
        Y11=Y11-sparse(1:nbus,1:nbus,accumarray(zipIdx,LHS_MatZip_p(:,1),[nbus,1]),nbus,nbus);
        Y12=Y12-sparse(1:nbus,1:nbus,accumarray(zipIdx,LHS_MatZip_p(:,2),[nbus,1]),nbus,nbus);
        Y21=Y21-sparse(1:nbus,1:nbus,accumarray(zipIdx,LHS_MatZip_p(:,3),[nbus,1]),nbus,nbus);
        Y22=Y22-sparse(1:nbus,1:nbus,accumarray(zipIdx,LHS_MatZip_p(:,4),[nbus,1]),nbus,nbus);
    end
    YLHS=[Y11,Y12;Y21,Y22];
    
    if ~isempty(syn)
        YLHS=YLHS+MatGCD_p;
    end
    
    idxNonSw=find(busType~=2);
    idxStackMat=[idxNonSw;idxNonSw+nbus];
    
    %
    % if nbus<=500
    %     [L_LHS_mat,U_LHS_mat,p_LHS_mat]=lu(LHS_mat,'vector');
    % end
    %
    % p_amd = colamd (LHS_mat) ;
    % MxI = speye (size(LHS_mat)) ;
    % MxQ = MxI (:, p_amd) ;
    % [MxL,MxU,MxP] = lu (LHS_mat*MxQ) ;
    
    auxList{iSys}.nbus=nbus;
    auxList{iSys}.nline=nline;
    auxList{iSys}.busType=busType;
    auxList{iSys}.isw=isw;
    auxList{iSys}.ipv=ipv;
    auxList{iSys}.ipq=ipq;
    auxList{iSys}.npv=npv;
    auxList{iSys}.npq=npq;
    auxList{iSys}.Y=Y;
    auxList{iSys}.V=V;
    auxList{iSys}.W=W;
    auxList{iSys}.P=P;
    auxList{iSys}.Q=Q;
    auxList{iSys}.IiL=IiL;
    auxList{iSys}.BiL=BiL;
    
    if iSys>=2
        busTag=zeros(nbus,1);
        busTag(sysMaps{iSys}.extBusMap(:,4))=1;
        busTagx=busTag(busType~=2);
        nY=npv+npq;
        
        busIdxB=find(busTagx==1);
        busIdxI=find(busTagx==0);
        
        busIdxPurgeB=find(busTag==0&busType~=2);
        busIdxRemainB=find(busTag==1&busType~=2);
        nPurgeB=size(busIdxPurgeB,1);
        
        LHS_mat=[YLHS(idxStackMat,idxStackMat),...
            [YEF11(busType~=2,busIdxPurgeB),YEF12(busType~=2,busIdxPurgeB),-F0M(busType~=2,ipv);...
            YEF21(busType~=2,busIdxPurgeB),YEF22(busType~=2,busIdxPurgeB),-E0M(busType~=2,ipv)];...
            C0M(ipv,busType~=2),D0M(ipv,busType~=2),sparse(npv,2*nPurgeB+npv);...
            E0M(busIdxPurgeB,busType~=2),-F0M(busIdxPurgeB,busType~=2),C0M(busIdxPurgeB,busIdxPurgeB),-D0M(busIdxPurgeB,busIdxPurgeB),sparse(nPurgeB,npv);...
            F0M(busIdxPurgeB,busType~=2),E0M(busIdxPurgeB,busType~=2),D0M(busIdxPurgeB,busIdxPurgeB),C0M(busIdxPurgeB,busIdxPurgeB),sparse(nPurgeB,npv);];
        
        % Note that for sub-systems, the E and F terms for boundary
        % variables will be purged.
        
        lhsTag=zeros(size(LHS_mat,1),1);
        lhsTag([busIdxB;busIdxB+nY])=1;
        
        LHS_bb=LHS_mat(lhsTag==1,lhsTag==1);
        LHS_bi=LHS_mat(lhsTag==1,lhsTag==0);
        LHS_ib=LHS_mat(lhsTag==0,lhsTag==1);
        LHS_ii=LHS_mat(lhsTag==0,lhsTag==0);
        
        p_amd = colamd (LHS_ii) ;
        MxI = speye (size(LHS_ii)) ;
        MxQ = MxI (:, p_amd) ;
        [MxL,MxU,MxP] = lu (LHS_ii*MxQ) ;
        
        LHS_biii=(((LHS_bi*MxQ)/MxU)/MxL)*MxP;
        LHS_iiib=MxQ*(MxU\(MxL\(MxP*LHS_ib)));
        LHS_reduce_bb=LHS_bb-LHS_biii*LHS_ib;
        
        auxList{iSys}.MxQ=MxQ;
        auxList{iSys}.MxL=MxL;
        auxList{iSys}.MxU=MxU;
        auxList{iSys}.MxP=MxP;
        auxList{iSys}.busIdxPurgeB=busIdxPurgeB;
        auxList{iSys}.busIdxRemainB=busIdxRemainB;
        auxList{iSys}.LHS_bb=LHS_bb;
        auxList{iSys}.LHS_bi=LHS_bi;
        auxList{iSys}.LHS_ib=LHS_ib;
        auxList{iSys}.LHS_ii=LHS_ii;
        auxList{iSys}.LHS_biii=LHS_biii;
        auxList{iSys}.LHS_iiib=LHS_iiib;
        auxList{iSys}.LHS_reduce_bb=LHS_reduce_bb;
    end
end

iSys=1;

[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg]=unfoldSysData(SysDataList{iSys});
busInX=[sysMaps{iSys}.busMap';sysMaps{iSys}.extBusMap(:,2)];
V0a=V0(busInX);
Q0a=Q0(busInX);

nbus=size(bus,1);
nline=size(line,1);
[pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,MatGV0,MatGV1,MatGRhs0,MatGRhs1]=unfoldSysPara(SysParaList{iSys});
%
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
if ~isempty(zip)%zipMode=0
    Ysh0=Ysh0+accumarray(zip(:,1),Rzip0.*(zip(:,5)+1j*zip(:,8)).*zip(:,12),[nbus,1]);
    Ysh1=Ysh1+accumarray(zip(:,1),Rzip1.*(zip(:,5)+1j*zip(:,8)).*zip(:,12),[nbus,1]);
end

Y=Ytr+sparse(1:nbus,1:nbus,Ysh0,nbus,nbus);

pVec=zeros(nbus,1);
qVec=zeros(nbus,1);
%     vSp=zeros(nbus,1);

pVec(pv(:,1))=pVec(pv(:,1))+pv(:,4);
pVec(pq(:,1))=pVec(pq(:,1))-pq(:,4);
qVec(pq(:,1))=qVec(pq(:,1))-pq(:,5);
if ~isempty(zip)%zipMode=0, account for the PQ components in ZIP loads
    pVec=pVec-accumarray(zip(:,1),Rzip0.*zip(:,7).*zip(:,12),[nbus,1]);
    qVec=qVec-accumarray(zip(:,1),Rzip0.*zip(:,10).*zip(:,12),[nbus,1]);
end
qVec(ipv)=Q0a(ipv);
%     vSp(ipv)=pv(:,5);

V=zeros(nbus,nlvl+1);
V(:,1)=V0a;
W=zeros(nbus,nlvl+1);
W(:,1)=1./V0a;
P=zeros(nbus,nlvl+1);
P(:,1)=pVec;
%     P(isw,2:end)=0;
Q=zeros(nbus,nlvl+1);
Q(:,1)=qVec;
P(:,2:(size(Pls,2)+1))=-Pls;
Q(:,2:(size(Qls,2)+1))=-Qls;
if ~isempty(zip)
    P(:,2)=P(:,2)-accumarray(zip(:,1),Rzip1.*zip(:,7).*zip(:,12),[nbus,1]);
    Q(:,2)=Q(:,2)-accumarray(zip(:,1),Rzip1.*zip(:,10).*zip(:,12),[nbus,1]);
end
Q(busType~=0,2:end)=0;

C0=real(V(:,1));
D0=imag(V(:,1));
E0=real(W(:,1));
F0=imag(W(:,1));

C0M=sparse(1:nbus,1:nbus,C0,nbus,nbus);
D0M=sparse(1:nbus,1:nbus,D0,nbus,nbus);
E0M=sparse(1:nbus,1:nbus,E0,nbus,nbus);
F0M=sparse(1:nbus,1:nbus,F0,nbus,nbus);

P0x=P(:,1);Q0x=Q(:,1);
if iSys==1
    for iSysi=2:nSys
        P0x(sysMaps{iSysi}.extBusMap(:,2))=P0x(sysMaps{iSysi}.extBusMap(:,2))+auxList{iSysi}.P(sysMaps{iSysi}.extBusMap(:,4),1);
        Q0x(sysMaps{iSysi}.extBusMap(:,2))=Q0x(sysMaps{iSysi}.extBusMap(:,2))+auxList{iSysi}.Q(sysMaps{iSysi}.extBusMap(:,4),1);
    end
end
P0M=sparse(1:nbus,1:nbus,P0x,nbus,nbus);
Q0M=sparse(1:nbus,1:nbus,Q0x,nbus,nbus);

G=real(Y);
B=imag(Y);

nZip=size(zip,1);
zipIdx=zip(:,1);
if ~isempty(zip)
    IiL=zeros(nZip,nlvl+1);
    BiL=zeros(nZip,nlvl+1);
    
    Bi0=abs(V0a(zipIdx));
    JI=zip(:,6);
    KI=-zip(:,9);
    
    Ii0L=Rzip0.*(JI+1j*KI).*V0a(zipIdx)./Bi0;
    Ji0L=real(Ii0L);
    Ki0L=imag(Ii0L);
    
    IiL(:,1)=Ii0L;
    BiL(:,1)=Bi0;
    
    Ci0=real(V0a(zipIdx));
    Di0=imag(V0a(zipIdx));
    
    LHS_MatZip=[Rzip0.*JI./Bi0-Ci0.*Ji0L./Bi0./Bi0,-Rzip0.*KI./Bi0-Di0.*Ji0L./Bi0./Bi0,...
        Rzip0.*KI./Bi0-Ci0.*Ki0L./Bi0./Bi0,Rzip0.*JI./Bi0-Di0.*Ki0L./Bi0./Bi0];
    Mat_BZip=[Ci0./Bi0,Di0./Bi0];
    
    auxList{iSys}.LHS_MatZip=LHS_MatZip;
    auxList{iSys}.Mat_BZip=Mat_BZip;
else
    IiL=zeros(nZip,nlvl+1);
    BiL=zeros(nZip,nlvl+1);
end

synIdx=syn(:,1);
if ~isempty(syn)
    MatGCD=-[sparse(synIdx,synIdx,MatGV0(:,1),nbus,nbus),sparse(synIdx,synIdx,MatGV0(:,2),nbus,nbus);...
        sparse(synIdx,synIdx,MatGV0(:,3),nbus,nbus),sparse(synIdx,synIdx,MatGV0(:,4),nbus,nbus)];
end

Y11=-G;Y12=B;Y21=-B;Y22=-G;
YEF11=P0M;YEF12=-Q0M;YEF21=-Q0M;YEF22=-P0M;

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

idxNonSw=find(busType~=2);
idxStackMat=[idxNonSw;idxNonSw+nbus];

%
% if nbus<=500
%     [L_LHS_mat,U_LHS_mat,p_LHS_mat]=lu(LHS_mat,'vector');
% end
%
% p_amd = colamd (LHS_mat) ;
% MxI = speye (size(LHS_mat)) ;
% MxQ = MxI (:, p_amd) ;
% [MxL,MxU,MxP] = lu (LHS_mat*MxQ) ;

auxList{iSys}.nbus=nbus;
auxList{iSys}.nline=nline;
auxList{iSys}.busType=busType;
auxList{iSys}.isw=isw;
auxList{iSys}.ipv=ipv;
auxList{iSys}.ipq=ipq;
auxList{iSys}.npv=npv;
auxList{iSys}.npq=npq;
auxList{iSys}.Y=Y;
auxList{iSys}.V=V;
auxList{iSys}.W=W;
auxList{iSys}.P=P;
auxList{iSys}.Q=Q;
auxList{iSys}.IiL=IiL;
auxList{iSys}.BiL=BiL;
for iSysi=2:nSys
    idxOverlay=sysMaps{iSysi}.extBusMap(:,2);
    idxOverlayMat=[idxOverlay;idxOverlay+nbus];
    YLHS(idxOverlayMat,idxOverlayMat)=YLHS(idxOverlayMat,idxOverlayMat)+auxList{iSysi}.LHS_reduce_bb;
end

LHS_mat=[YLHS(idxStackMat,idxStackMat),...
    [YEF11(busType~=2,busType~=2),YEF12(busType~=2,busType~=2),-F0M(busType~=2,ipv);...
    YEF21(busType~=2,busType~=2),YEF22(busType~=2,busType~=2),-E0M(busType~=2,ipv)];...
    C0M(ipv,busType~=2),D0M(ipv,busType~=2),sparse(npv,2*npq+3*npv);...
    E0M(busType~=2,busType~=2),-F0M(busType~=2,busType~=2),C0M(busType~=2,busType~=2),-D0M(busType~=2,busType~=2),sparse(npq+npv,npv);...
    F0M(busType~=2,busType~=2),E0M(busType~=2,busType~=2),D0M(busType~=2,busType~=2),C0M(busType~=2,busType~=2),sparse(npq+npv,npv);];

p_amd = colamd (LHS_mat) ;
MxI = speye (size(LHS_mat)) ;
MxQ = MxI (:, p_amd) ;
[MxL,MxU,MxP] = lu (LHS_mat*MxQ) ;

auxList{iSys}.MxQ=MxQ;
auxList{iSys}.MxL=MxL;
auxList{iSys}.MxU=MxU;
auxList{iSys}.MxP=MxP;
% Merge the blocks of sub-systems into the main system's matrix

for i=1:nlvl
    seq2=getseq(i,2);
%     seq2p=getseq(i+1,2);
%     seq3=getseq(i,3);
    idxSeq2=sum(seq2==i,2);
    seq2R=seq2(idxSeq2==0,:);
    
    for iSys=2:nSys
        [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg]=unfoldSysData(SysDataList{iSys});
        [pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,MatGV0,MatGV1,MatGRhs0,MatGRhs1]=unfoldSysPara(SysParaList{iSys});
        nbus=auxList{iSys}.nbus;
        busType=auxList{iSys}.busType;
        IiL=auxList{iSys}.IiL;
        BiL=auxList{iSys}.BiL;
        V=auxList{iSys}.V;
        W=auxList{iSys}.W;
        P=auxList{iSys}.P;
        Q=auxList{iSys}.Q;
        Y=auxList{iSys}.Y;
        ipv=auxList{iSys}.ipv;
        isw=auxList{iSys}.isw;
        busIdxPurgeB=auxList{iSys}.busIdxPurgeB;
        
        RHSILr=zeros(nbus,1);
        RHSILi=zeros(nbus,1);
        
        RHSIiLr=zeros(nbus,1);
        RHSIiLi=zeros(nbus,1);
        if ~isempty(zip)
            zipIdx=zip(:,1);
            JI=zip(:,6);
            KI=-zip(:,9);
            Ji0L=real(IiL(:,1));
            Ki0L=imag(IiL(:,1));
            Bi0=BiL(:,1);
            
            RHS_BZip_p=(real(sum(V(zipIdx,seq2R(:,1)+1).*conj(V(zipIdx,seq2R(:,2)+1)),2))-sum(BiL(:,seq2R(:,1)+1).*BiL(:,seq2R(:,2)+1),2))./Bi0/2;
            RHZ_BIConv=sum(IiL(:,seq2R(:,1)+1).*BiL(:,seq2R(:,2)+1),2);
            RHSILr_full_p=Rzip1.*(JI.*real(V(zipIdx,i))-KI.*imag(V(zipIdx,i)))./Bi0-real(RHZ_BIConv)./Bi0-Ji0L.*RHS_BZip_p./Bi0;
            RHSILi_full_p=Rzip1.*(KI.*real(V(zipIdx,i))+JI.*imag(V(zipIdx,i)))./Bi0-imag(RHZ_BIConv)./Bi0-Ki0L.*RHS_BZip_p./Bi0;
            RHSIiLr=accumarray(zipIdx,RHSILr_full_p,[nbus,1]);
            RHSIiLi=accumarray(zipIdx,RHSILi_full_p,[nbus,1]);
            
            auxList{iSys}.RHSILr_full=RHSILr_full_p;
            auxList{iSys}.RHSILi_full=RHSILi_full_p;
            auxList{iSys}.RHS_BZip=RHS_BZip_p;
        end
        
        RHSIGr=zeros(nbus,1);
        RHSIGi=zeros(nbus,1);
        if ~isempty(syn)
            synIdx=syn(:,1);
            RHSIGr=-accumarray(synIdx,MatGV1(:,1).*real(V(synIdx,i))+MatGV1(:,2).*imag(V(synIdx,i)),[nbus,1]);
            RHSIGi=-accumarray(synIdx,MatGV1(:,3).*real(V(synIdx,i))+MatGV1(:,4).*imag(V(synIdx,i)),[nbus,1]);
            if i==1
                RHSIGr=RHSIGr+accumarray(synIdx,MatGRhs1(:,1),[nbus,1]);
                RHSIGi=RHSIGi+accumarray(synIdx,MatGRhs1(:,2),[nbus,1]);
            end
        end
        
        % HEM Body
        RHS1=sum((-P(:,seq2(:,1)+1)+1j*Q(:,seq2(:,1)+1)).*conj(W(:,seq2(:,2)+1)),2)+Ysh1.*V(:,i)+Ytr1*V(:,i);
        RHS2=-0.5*real(sum(V(:,seq2R(:,1)+1).*conj(V(:,seq2R(:,2)+1)),2));
        RHS3=sum(-W(:,seq2R(:,1)+1).*V(:,seq2R(:,2)+1),2);
        
        if i==1
            RHS2=RHS2+0.5*VspSq2(:,2);
        end
        
        compactRHS1=RHS1(busIdxPurgeB);
        if ~isempty(isw);compactRHS1=compactRHS1+Y(busIdxPurgeB,isw)*V(isw,i+1);end
        RHS=[real(compactRHS1)+RHSILr(busIdxPurgeB)+RHSIiLr(busIdxPurgeB)-RHSIGr(busIdxPurgeB);...
            imag(compactRHS1)+RHSILi(busIdxPurgeB)+RHSIiLi(busIdxPurgeB)-RHSIGi(busIdxPurgeB);...
            RHS2(ipv);...
            real(RHS3(busIdxPurgeB));...
            imag(RHS3(busIdxPurgeB));];
        
        auxList{iSys}.RHS=RHS;
    end
    
    iSys=1;
    [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg]=unfoldSysData(SysDataList{iSys});
    [pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,MatGV0,MatGV1,MatGRhs0,MatGRhs1]=unfoldSysPara(SysParaList{iSys});
    nbus=auxList{iSys}.nbus;
    busType=auxList{iSys}.busType;
    IiL=auxList{iSys}.IiL;
    BiL=auxList{iSys}.BiL;
    V=auxList{iSys}.V;
    W=auxList{iSys}.W;
    P=auxList{iSys}.P;
    Q=auxList{iSys}.Q;
    Y=auxList{iSys}.Y;
    ipv=auxList{iSys}.ipv;
    isw=auxList{iSys}.isw;
    npv=auxList{iSys}.npv;
    npq=auxList{iSys}.npq;
    
    MxQ=auxList{iSys}.MxQ;
    MxL=auxList{iSys}.MxL;
    MxU=auxList{iSys}.MxU;
    MxP=auxList{iSys}.MxP;
    
    RHSILr=zeros(nbus,1);
    RHSILi=zeros(nbus,1);
    
    RHSIiLr=zeros(nbus,1);
    RHSIiLi=zeros(nbus,1);
    if ~isempty(zip)
        zipIdx=zip(:,1);
        JI=zip(:,6);
        KI=-zip(:,9);
        Ji0L=real(IiL(:,1));
        Ki0L=imag(IiL(:,1));
        Bi0=BiL(:,1);
        
        RHS_BZip=(real(sum(V(zipIdx,seq2R(:,1)+1).*conj(V(zipIdx,seq2R(:,2)+1)),2))-sum(BiL(:,seq2R(:,1)+1).*BiL(:,seq2R(:,2)+1),2))./Bi0/2;
        RHZ_BIConv=sum(IiL(:,seq2R(:,1)+1).*BiL(:,seq2R(:,2)+1),2);
        RHSILr_full=Rzip1.*(JI.*real(V(zipIdx,i))-KI.*imag(V(zipIdx,i)))./Bi0-real(RHZ_BIConv)./Bi0-Ji0L.*RHS_BZip./Bi0;
        RHSILi_full=Rzip1.*(KI.*real(V(zipIdx,i))+JI.*imag(V(zipIdx,i)))./Bi0-imag(RHZ_BIConv)./Bi0-Ki0L.*RHS_BZip./Bi0;
        RHSIiLr=accumarray(zipIdx,RHSILr_full,[nbus,1]);
        RHSIiLi=accumarray(zipIdx,RHSILi_full,[nbus,1]);
        
        auxList{iSys}.RHSILr_full=RHSILr_full;
        auxList{iSys}.RHSILi_full=RHSILi_full;
        auxList{iSys}.RHS_BZip=RHS_BZip;
    end
    
    RHSIGr=zeros(nbus,1);
    RHSIGi=zeros(nbus,1);
    if ~isempty(syn)
        synIdx=syn(:,1);
        RHSIGr=-accumarray(synIdx,MatGV1(:,1).*real(V(synIdx,i))+MatGV1(:,2).*imag(V(synIdx,i)),[nbus,1]);
        RHSIGi=-accumarray(synIdx,MatGV1(:,3).*real(V(synIdx,i))+MatGV1(:,4).*imag(V(synIdx,i)),[nbus,1]);
        if i==1
            RHSIGr=RHSIGr+accumarray(synIdx,MatGRhs1(:,1),[nbus,1]);
            RHSIGi=RHSIGi+accumarray(synIdx,MatGRhs1(:,2),[nbus,1]);
        end
    end
    
    % HEM Body
    RHS1=sum((-P(:,seq2(:,1)+1)+1j*Q(:,seq2(:,1)+1)).*conj(W(:,seq2(:,2)+1)),2)+Ysh1.*V(:,i)+Ytr1*V(:,i);
    RHS2=-0.5*real(sum(V(:,seq2R(:,1)+1).*conj(V(:,seq2R(:,2)+1)),2));
    RHS3=sum(-W(:,seq2R(:,1)+1).*V(:,seq2R(:,2)+1),2);
    
    if i==1
        RHS2=RHS2+0.5*VspSq2(:,2);
    end
    
    compactRHS1=RHS1(busType~=2);
    if ~isempty(isw);compactRHS1=compactRHS1+Y(busType~=2,isw)*V(isw,i+1);end
    
    RHSIiLrExt=zeros(nbus,1);
    RHSIiLiExt=zeros(nbus,1);
    for iSys=2:nSys
        RHSIiLExtCompact=auxList{iSys}.LHS_biii*auxList{iSys}.RHS;
        RHSIiLrExt(sysMaps{iSys}.extBusMap(:,2))=RHSIiLrExt(sysMaps{iSys}.extBusMap(:,2))+RHSIiLExtCompact(1:(size(RHSIiLExtCompact,1)/2));
        RHSIiLiExt(sysMaps{iSys}.extBusMap(:,2))=RHSIiLiExt(sysMaps{iSys}.extBusMap(:,2))+RHSIiLExtCompact((size(RHSIiLExtCompact,1)/2+1):end);
    end
    
    RHS=[real(compactRHS1)+RHSILr(busType~=2)+RHSIiLr(busType~=2)-RHSIGr(busType~=2)-RHSIiLrExt(busType~=2);...
        imag(compactRHS1)+RHSILi(busType~=2)+RHSIiLi(busType~=2)-RHSIGi(busType~=2)-RHSIiLiExt(busType~=2);...
        RHS2(ipv);...
        real(RHS3(busType~=2));...
        imag(RHS3(busType~=2));];
    % Solve linear
    x = MxQ * (MxU \ (MxL \ (MxP * RHS))) ;
    
    V(busType~=2,i+1)=x(1:(npq+npv))+1j*x(((npq+npv)+1):(2*(npq+npv)));
    W(busType~=2,i+1)=x((2*(npq+npv)+1):(3*(npq+npv)))+1j*x((3*(npq+npv)+1):(4*(npq+npv)));
    Q(ipv,i+1)=x((4*(npq+npv)+1):end);
    
    iSys=1;
    if ~isempty(zip)
        zipIdx=zip(:,1);
        LHS_MatZip=auxList{iSys}.LHS_MatZip;
        Mat_BZip=auxList{iSys}.Mat_BZip;
        IiL(:,i+1)=(LHS_MatZip(:,1)+1j*LHS_MatZip(:,3)).*real(V(zipIdx,i+1))+(LHS_MatZip(:,2)+1j*LHS_MatZip(:,4)).*imag(V(zipIdx,i+1))+(RHSILr_full+1j*RHSILi_full);
        BiL(:,i+1)=Mat_BZip(:,1).*real(V(zipIdx,i+1))+Mat_BZip(:,2).*imag(V(zipIdx,i+1))+RHS_BZip;
    end
    auxList{iSys}.V=V;
    auxList{iSys}.W=W;
    auxList{iSys}.Q=Q;
    auxList{iSys}.IiL=IiL;
    auxList{iSys}.BiL=BiL;
    
    busMap=sysMaps{iSys}.busMap';
    Vm(busMap,i+1)=V(1:size(busMap,1),i+1);
    Wm(busMap,i+1)=W(1:size(busMap,1),i+1);
    Qm(busMap,i+1)=Q(1:size(busMap,1),i+1);
    
    % Finished main, start sub-systems
    auxList_p=auxList;
    for iSys=2:nSys
        aux_p=auxList_p{iSys};
        [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg]=unfoldSysData(SysDataList{iSys});
        [pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,MatGV0,MatGV1,MatGRhs0,MatGRhs1]=unfoldSysPara(SysParaList{iSys});
        nbus=aux_p.nbus;
        busType=aux_p.busType;
        IiL=aux_p.IiL;
        BiL=aux_p.BiL;
        V=aux_p.V;
        W=aux_p.W;
        Q=aux_p.Q;
        Y=aux_p.Y;
        ipv=aux_p.ipv;
        isw=aux_p.isw;
        npv=aux_p.npv;
        npq=aux_p.npq;
        
        busIdxPurgeB=aux_p.busIdxPurgeB;
        nPurgeB=size(busIdxPurgeB,1);
        busIdxRemainB=aux_p.busIdxRemainB;
        
        MxQ=aux_p.MxQ;
        MxL=aux_p.MxL;
        MxU=aux_p.MxU;
        MxP=aux_p.MxP;
        
        RHS=aux_p.RHS;
        
        Vb=auxList_p{1}.V(sysMaps{iSys}.extBusMap(:,2),i+1);
        Wb=auxList_p{1}.W(sysMaps{iSys}.extBusMap(:,2),i+1);
        
        % Solve linear
        x = MxQ * (MxU \ (MxL \ (MxP * RHS)))-auxList{iSys}.LHS_iiib*[real(Vb);imag(Vb)];
        
        V(busIdxPurgeB,i+1)=x(1:nPurgeB)+1j*x((nPurgeB+1):(2*nPurgeB));
        W(busIdxPurgeB,i+1)=x((2*nPurgeB+1):(3*nPurgeB))+1j*x((3*nPurgeB+1):(4*nPurgeB));
        Q(ipv,i+1)=x((4*nPurgeB+1):end);
        V(busIdxRemainB,i+1)=Vb;
        W(busIdxRemainB,i+1)=Wb;
        
        if ~isempty(zip)
            zipIdx=zip(:,1);
            RHSILr_full_p=auxList{iSys}.RHSILr_full;
            RHSILi_full_p=auxList{iSys}.RHSILi_full;
            RHS_BZip_p=auxList{iSys}.RHS_BZip;
            LHS_MatZip=auxList{iSys}.LHS_MatZip;
            Mat_BZip=auxList{iSys}.Mat_BZip;
            IiL(:,i+1)=(LHS_MatZip(:,1)+1j*LHS_MatZip(:,3)).*real(V(zipIdx,i+1))+(LHS_MatZip(:,2)+1j*LHS_MatZip(:,4)).*imag(V(zipIdx,i+1))+(RHSILr_full_p+1j*RHSILi_full_p);
            BiL(:,i+1)=Mat_BZip(:,1).*real(V(zipIdx,i+1))+Mat_BZip(:,2).*imag(V(zipIdx,i+1))+RHS_BZip_p;
        end
        
        aux_p.V=V;
        aux_p.W=W;
        aux_p.Q=Q;
        aux_p.IiL=IiL;
        aux_p.BiL=BiL;
        
        auxList{iSys}=aux_p;
    end
    for iSys=2:nSys
        busMap=sysMaps{iSys}.busMap';
        Vm(busMap,i+1)=auxList{iSys}.V(1:size(busMap,1),i+1);
        Wm(busMap,i+1)=auxList{iSys}.W(1:size(busMap,1),i+1);
        Qm(busMap,i+1)=auxList{iSys}.Q(1:size(busMap,1),i+1);
    end
end
end