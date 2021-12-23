function [V,W,Q,s,d,JG,KG]=hemMachinePFSalientcontinueSep(SimData,SysData,SysPara,x,pShare,Vsp)
% HE computation of power flow (extended including synchronous generators)
%
% FUNCTION hemMachinePFSalientcontinueSep
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
% INPUT (will be consolidated in a future version)
%
% OUTPUT (will be consolidated in a future version)
%

global IS_OCTAVE;

[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
[V0,Q0,s0,d0,w,eq1,eq2,ed1,ed2,psid,psiq,Pm0,Ef0,Vavrm,Vavrr,Vavrf,Vavrref,tgovg,tgovm,tgovmech,f,dpg,qplt,vg]=unfoldX(x,SysData);
[pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr,~,Ysh0,Ysh1,VspSq2,~,~,~,~,Tmech1,Varref1,Ef1,Pm1,Eq11]=unfoldSysPara(SysPara);
[maxAlpha,segAlpha,dAlpha,nlvl,taylorN,alphaTol,diffTol,diffTolMax,method]=unfoldSimData(SimData);

nbus=size(bus,1);
nline=size(line,1);
if isfield(SysPara,'nIslands')&&isfield(SysPara,'islands')&&isfield(SysPara,'refs')
    nIslands=SysPara.nIslands;islands=SysPara.islands;refs=SysPara.refs;
else
    [nIslands,islands,refs]=searchIslands(bus(:,1),line(:,1:2));
end

Pls=zeros(nbus,1);Pls(pq(:,1))=pqIncr(:,1);if ~isempty(pv);Pls(pv(:,1))=Pls(pv(:,1))-pvIncr;end
Qls=zeros(nbus,1);Qls(pq(:,1))=pqIncr(:,2);

[Y,Ytr,Ysh,ytrfr,ytrto,yshfr,yshto]=getYMatrix(nbus,line);

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
busTypePvPq=busType(busType~=2);
ipvInPvPq=find(busTypePvPq==1);
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
% qVec(ipv)=Q0(ipv);
%     vSp(ipv)=pv(:,5);

V=zeros(nbus,nlvl+1);
V(:,1)=V0;
V(isw,2)=Vsp(isw);
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
Pspec=P(:,1);
Qspec=Q(:,1)+Qxtra(:,1);
V0sq=C0.*C0+D0.*D0;

C0M=sparse(1:nbus,1:nbus,C0,nbus,nbus);
D0M=sparse(1:nbus,1:nbus,D0,nbus,nbus);
E0M=sparse(1:nbus,1:nbus,E0,nbus,nbus);
F0M=sparse(1:nbus,1:nbus,F0,nbus,nbus);
P0M=sparse(1:nbus,1:nbus,Pspec,nbus,nbus);
Q0M=sparse(1:nbus,1:nbus,Qspec,nbus,nbus);

G=real(Y);
B=imag(Y);

if ~isempty(ind)
    nInd=size(ind,1);
    indIdx=ind(:,1);
    s=zeros(nInd,nlvl+1);
    s(:,1)=s0;
    IL=zeros(nInd,nlvl+1);
    IR=zeros(nInd,nlvl+1);
    Vm=zeros(nInd,nlvl+1);
    
    R1=ind(:,7);
    X1=ind(:,8);
    Z1=ind(:,7)+1j*ind(:,8);
    Ze=1j*ind(:,13);
    R2=ind(:,9);
    X2=ind(:,10);
    T0=ind(:,15)+ind(:,16)+ind(:,17);
    T1=-ind(:,16)-2*ind(:,17);
    T2=ind(:,17);
    
    Am=sparse(indIdx,(1:nInd)',ones(1,nInd),nbus,nInd);
    
    IL(:,1)=V0(indIdx)./(Z1+Ze.*(R2+1j*X2.*s0)./(R2.*Reind0+(1j*X2.*Reind0+Ze).*s0));
    Vm(:,1)=V0(indIdx)-IL(:,1).*Z1;
    IR(:,1)=Vm(:,1).*s0./(R2+1j*X2.*s0);
    
    J0=real(IR(:,1));
    K0=imag(IR(:,1));
    JL0=real(IL(:,1));
    KL0=imag(IL(:,1));
    
    Yeind0=Reind0./Ze;
    Yeind1=Reind1./Ze;
    Ye1ind0=Reind0.*Z1./Ze;
    Ye1ind1=Reind1.*Z1./Ze;
    Ge=real(Yeind0);
    Be=imag(Yeind0);
    kg1e=real(Ye1ind0);
    kb1e=imag(Ye1ind0);
    Ge1=real(Yeind1);
    Be1=imag(Yeind1);
    kg1e1=real(Ye1ind1);
    kb1e1=imag(Ye1ind1);
    
    % The structure:
    % |C|D||x| |E|
    % |A|B||y|=|F|
    % Ly=(D-CA^-1B)y=E-CA^-1F
    % x=-A^-1By+A^-1F
    % x=[J;K;s];y=[JL;KL;C;D];
    % L=[R,S][y1;y2];
    % y1=[JL;KL];y2=[C;D];
    % y1=-R\S*y2+R\(E-CA^-1F)
    %
    %         idxE=[1,2,5];
    %         idxM=[3,4,6,7];
    %
    %         LHS_MatInd_Shr=zeros(nInd,2,2);
    %         RHS_C_Shr=zeros(nInd,2,5);
    %         LHS_MatInd_Shr2=zeros(nInd,3,4); % A^-1B
    %         LHS_MatInd_Shr3=zeros(nInd,3,3); % A^-1
    %
    %         LHS_MatInd_Shr_1=zeros(nInd,2,2);
    %         RHS_C_Shr_1=zeros(nInd,2,5);
    %         LHS_MatInd_Shr2_1=zeros(nInd,3,4); % A^-1B
    %         LHS_MatInd_Shr3_1=zeros(nInd,3,3); % A^-1
    
    %         for i=1:nInd
    %             if Rind0(i)~=0
    %                 LHS_MatInd=[R2(i),-X2(i)*s0(i),R1(i)*s0(i),-X1(i)*s0(i),-K0(i)*X2(i)-C0(indIdx(i))+JL0(i)*R1(i)-KL0(i)*X1(i),-s0(i),0;...
    %                     X2(i)*s0(i), R2(i),X1(i)*s0(i), R1(i)*s0(i), J0(i)*X2(i)-D0(indIdx(i))+JL0(i)*X1(i)+KL0(i)*R1(i),0,-s0(i);...
    %                     2*J0(i),2*K0(i),0,0,-(Rind0(i)*T0(i)+2*Rind0(i)*T1(i)*s0(i)+3*Rind0(i)*T2(i)*s0(i)*s0(i))/R2(i),0,0;...
    %                     -1,0,1+kg1e(i),-kb1e(i),0,-Ge(i), Be(i);...
    %                     0,-1,kb1e(i), 1+kg1e(i),0,-Be(i),-Ge(i);];
    %                 temp0=LHS_MatInd([3,4,5],idxE)\eye(3);                       % A^-1
    %                 LHS_MatInd_Shr2(i,:,:)=temp0*LHS_MatInd([3,4,5],idxM);       % A^-1B
    %                 LHS_MatInd_Shr3(i,:,:)=temp0;                                % A^-1
    %                 temp1=LHS_MatInd([1,2],idxE)/LHS_MatInd([3,4,5],idxE);       % CA^-1
    %                 temp2=LHS_MatInd([1,2],idxM)-temp1*LHS_MatInd([3,4,5],idxM); % L=D-CA^-1B
    %                 LHS_MatInd_Shr(i,:,:)=-temp2(:,[1,2])\temp2(:,[3,4]);        % -R\S
    %                 RHS_C_Shr(i,:,:)=temp2(:,[1,2])\[eye(2),-temp1];             % R\[I,-CA^-1]
    %
    %                 LHS_MatInd_Shr_1(i,:,:)=LHS_MatInd_Shr(i,:,:);
    %                 RHS_C_Shr_1(i,:,:)=RHS_C_Shr(i,:,:);
    %                 LHS_MatInd_Shr2_1(i,:,:)=LHS_MatInd_Shr2(i,:,:);
    %                 LHS_MatInd_Shr3_1(i,:,:)=LHS_MatInd_Shr3(i,:,:);
    %             else % singular case
    %                 LHS_MatInd=[R2(i),-X2(i)*s0(i),R1(i)*s0(i),-X1(i)*s0(i),-K0(i)*X2(i)-C0(indIdx(i))+JL0(i)*R1(i)-KL0(i)*X1(i),-s0(i),0;...
    %                     X2(i)*s0(i), R2(i),X1(i)*s0(i), R1(i)*s0(i), J0(i)*X2(i)-D0(indIdx(i))+JL0(i)*X1(i)+KL0(i)*R1(i),0,-s0(i);...
    %                     2*real(IR(i,2)),2*imag(IR(i,2)),0,0,-(Rind1(i)*T0(i)+2*Rind1(i)*T1(i)*s0(i)+3*Rind1(i)*T2(i)*s0(i)*s0(i))/R2(i),0,0;...
    %                     -1,0,1+kg1e(i),-kb1e(i),0,-Ge(i), Be(i);...
    %                     0,-1,kb1e(i), 1+kg1e(i),0,-Be(i),-Ge(i);];
    %                 temp0=LHS_MatInd([3,4,5],idxE)\eye(3);                       % A^-1
    %                 LHS_MatInd_Shr2(i,:,:)=temp0*LHS_MatInd([3,4,5],idxM);       % A^-1B
    %                 LHS_MatInd_Shr3(i,:,:)=temp0;                                % A^-1
    %                 temp1=LHS_MatInd([1,2],idxE)/LHS_MatInd([3,4,5],idxE);       % CA^-1
    %                 temp2=LHS_MatInd([1,2],idxM)-temp1*LHS_MatInd([3,4,5],idxM); % L=D-CA^-1B
    %                 LHS_MatInd_Shr(i,:,:)=-temp2(:,[1,2])\temp2(:,[3,4]);        % -R\S
    %                 RHS_C_Shr(i,:,:)=temp2(:,[1,2])\[eye(2),-temp1];             % R\[I,-CA^-1]
    %
    %                 temp3=[1+kg1e(i),-kb1e(i);kb1e(i), 1+kg1e(i)];
    %                 temp4=[Ge(i),-Be(i);Be(i),Ge(i)];
    %                 LHS_MatInd_Shr_1(i,:,:)=temp3\temp4;        % -R\S (i=1)
    %                 RHS_C_Shr_1(i,:,:)=temp3\[eye(2),zeros(2,3)];             % R\[I,-CA^-1]  (i=1)
    %             end
    %         end
    %         LHS_MatInd_Bus=zeros(nbus,2,2);                                  % \sum{-R\S} by buses
    %         LHS_MatInd_Bus(:,1,1)=accumarray(indIdx,LHS_MatInd_Shr(:,1,1),[nbus,1]);
    %         LHS_MatInd_Bus(:,1,2)=accumarray(indIdx,LHS_MatInd_Shr(:,1,2),[nbus,1]);
    %         LHS_MatInd_Bus(:,2,1)=accumarray(indIdx,LHS_MatInd_Shr(:,2,1),[nbus,1]);
    %         LHS_MatInd_Bus(:,2,2)=accumarray(indIdx,LHS_MatInd_Shr(:,2,2),[nbus,1]);
    %
    %         LHS_MatInd_Bus_1=zeros(nbus,2,2);                                  % \sum{-R\S} by buses (i=1)
    %         LHS_MatInd_Bus_1(:,1,1)=accumarray(indIdx,LHS_MatInd_Shr_1(:,1,1),[nbus,1]);
    %         LHS_MatInd_Bus_1(:,1,2)=accumarray(indIdx,LHS_MatInd_Shr_1(:,1,2),[nbus,1]);
    %         LHS_MatInd_Bus_1(:,2,1)=accumarray(indIdx,LHS_MatInd_Shr_1(:,2,1),[nbus,1]);
    %         LHS_MatInd_Bus_1(:,2,2)=accumarray(indIdx,LHS_MatInd_Shr_1(:,2,2),[nbus,1]);
    LHS_MatInd_Full=zeros(nInd,5,2);
    LHS_MatInd_Shr=zeros(nInd,2,2);
    %         LHS_MatInd_Shr_1=zeros(nInd,2,2);
    RHS_C_Shr=zeros(nInd,5,5);
    for i=1:nInd
        LHS_MatInd=[R2(i),-X2(i)*s0(i),R1(i)*s0(i),-X1(i)*s0(i),-K0(i)*X2(i)-C0(indIdx(i))+JL0(i)*R1(i)-KL0(i)*X1(i),-s0(i),0;...
            X2(i)*s0(i), R2(i),X1(i)*s0(i), R1(i)*s0(i), J0(i)*X2(i)-D0(indIdx(i))+JL0(i)*X1(i)+KL0(i)*R1(i),0,-s0(i);...
            C0(indIdx(i))-JL0(i)*R1(i)+KL0(i)*X1(i),D0(indIdx(i))-KL0(i)*R1(i)-JL0(i)*X1(i),-J0(i)*R1(i)-K0(i)*X1(i),-K0(i)*R1(i)+J0(i)*X1(i),-Rind0(i)*(T1(i)+2*T2(i)*s0(i)),J0(i),K0(i);...
            -1,0,1+kg1e(i),-kb1e(i),0,-Ge(i), Be(i);...
            0,-1,kb1e(i), 1+kg1e(i),0,-Be(i),-Ge(i);];
        MatIndA=LHS_MatInd(:,1:5);
        MatIndB=LHS_MatInd(:,[6,7]);
        MatInvA=MatIndA\eye(5);
        MatCDtoRest=-MatInvA*MatIndB;
        
        RHS_C_Shr(i,:,:)=MatInvA;
        LHS_MatInd_Full(i,:,:)=MatCDtoRest;
        LHS_MatInd_Shr(i,:,:)=MatCDtoRest([3,4],:);
        %             LHS_MatInd_Shr_1(i,:,:)=LHS_MatInd_Shr(i,:,:);
        %             temp0=LHS_MatInd([3,4,5],idxE)\eye(3);                       % A^-1
        %             LHS_MatInd_Shr2(i,:,:)=temp0*LHS_MatInd([3,4,5],idxM);       % A^-1B
        %             LHS_MatInd_Shr3(i,:,:)=temp0;                                % A^-1
        %             temp1=LHS_MatInd([1,2],idxE)/LHS_MatInd([3,4,5],idxE);       % CA^-1
        %             temp2=LHS_MatInd([1,2],idxM)-temp1*LHS_MatInd([3,4,5],idxM); % L=D-CA^-1B
        %             LHS_MatInd_Shr(i,:,:)=-temp2(:,[1,2])\temp2(:,[3,4]);        % -R\S
        %             RHS_C_Shr(i,:,:)=temp2(:,[1,2])\[eye(2),-temp1];             % R\[I,-CA^-1]
        %
        %             LHS_MatInd_Shr_1(i,:,:)=LHS_MatInd_Shr(i,:,:);
        %             RHS_C_Shr_1(i,:,:)=RHS_C_Shr(i,:,:);
        %             LHS_MatInd_Shr2_1(i,:,:)=LHS_MatInd_Shr2(i,:,:);
        %             LHS_MatInd_Shr3_1(i,:,:)=LHS_MatInd_Shr3(i,:,:);
    end
    LHS_MatInd_Bus=zeros(nbus,2,2);                                  % \sum{-R\S} by buses
    LHS_MatInd_Bus(:,1,1)=accumarray(indIdx,LHS_MatInd_Shr(:,1,1),[nbus,1]);
    LHS_MatInd_Bus(:,1,2)=accumarray(indIdx,LHS_MatInd_Shr(:,1,2),[nbus,1]);
    LHS_MatInd_Bus(:,2,1)=accumarray(indIdx,LHS_MatInd_Shr(:,2,1),[nbus,1]);
    LHS_MatInd_Bus(:,2,2)=accumarray(indIdx,LHS_MatInd_Shr(:,2,2),[nbus,1]);
    
    %         LHS_MatInd_Bus_1=zeros(nbus,2,2);                                  % \sum{-R\S} by buses (i=1)
    %         LHS_MatInd_Bus_1(:,1,1)=accumarray(indIdx,LHS_MatInd_Shr_1(:,1,1),[nbus,1]);
    %         LHS_MatInd_Bus_1(:,1,2)=accumarray(indIdx,LHS_MatInd_Shr_1(:,1,2),[nbus,1]);
    %         LHS_MatInd_Bus_1(:,2,1)=accumarray(indIdx,LHS_MatInd_Shr_1(:,2,1),[nbus,1]);
    %         LHS_MatInd_Bus_1(:,2,2)=accumarray(indIdx,LHS_MatInd_Shr_1(:,2,2),[nbus,1]);
else
    s=[];
end

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
    nSyn=size(syn,1);
    synIdx=syn(:,1);
    Rs=syn(:,7);
    Xd=syn(:,8);
    Xq=syn(:,13);
    Ms=syn(:,18);
    Ds=syn(:,19);
    
    d=zeros(nSyn,nlvl+1);
    JG=zeros(nSyn,nlvl+1);
    KG=zeros(nSyn,nlvl+1);
    Cd=zeros(nSyn,nlvl+1);
    Sd=zeros(nSyn,nlvl+1);
    Ef=zeros(nSyn,nlvl+1);
    d(:,1)=d0;
    
    Efd=zeros(nSyn,1);
    Efq=Ef0;
    cosd=cos(d0);
    sind=sin(d0);
    Cg=C0(synIdx);
    Dg=D0(synIdx);
    Vd=sind.*Cg-cosd.*Dg;
    Vq=cosd.*Cg+sind.*Dg;
    Id=(Rs.*(Efd-Vd)+Xq.*(Efq-Vq))./(Rs.*Rs+Xq.*Xd);
    Iq=(-Xd.*(Efd-Vd)+Rs.*(Efq-Vq))./(Rs.*Rs+Xq.*Xd);
    IG0=(sind.*Id+cosd.*Iq)+1j*(-cosd.*Id+sind.*Iq);
    %         IG0=(Ef0.*(cos(d0)+1j*sin(d0))-V0(synIdx))./(Rs+1j*Xs);
    JG(:,1)=real(IG0);
    KG(:,1)=imag(IG0);
    Cd(:,1)=cos(d0);
    Sd(:,1)=sin(d0);
    Ef(:,[1,2])=[Ef0,Ef1];
    
    CG0=C0(synIdx);
    DG0=D0(synIdx);
    JG0=JG(:,1);
    KG0=KG(:,1);
    Cd0=Cd(:,1);
    Sd0=Sd(:,1);
    
    [cosp,sinp,taylorN]=getTaylorPolynomials(d0,taylorN); % taylorN may be truncated
    
    Pm=zeros(nSyn,nlvl+1);
    Pm(:,1)=Pm0;
    if ~isempty(Pm1)
        Pm(:,2)=Pm1;
    end
    
    A1n=zeros(nSyn,1);
    B1n=zeros(nSyn,1);
    for i=1:taylorN
        A1n=A1n.*d0+(taylorN+1-i)*cosp(:,taylorN+2-i);
        B1n=B1n.*d0+(taylorN+1-i)*sinp(:,taylorN+2-i);
    end
    
    synByIsland=islands(synIdx);
    sumShare=accumarray(synByIsland,pShare);
    nIsland=length(sumShare);
    pShare=pShare./sumShare(synByIsland);
    mAuxIsland=sparse(1:nSyn,synByIsland,pShare,nSyn,nIsland);
    [~,idxBal]=max(mAuxIsland);
    idxBalSyn=idxBal(synByIsland);
    
    MatG1C=sparse(1:nSyn,synIdx,Cd0,nSyn,nbus);
    MatG1D=sparse(1:nSyn,synIdx,Sd0,nSyn,nbus);
    MatG2C=sparse(1:nSyn,synIdx,Sd0,nSyn,nbus);
    MatG2D=sparse(1:nSyn,synIdx,-Cd0,nSyn,nbus);
    MatG5AC=sparse([1:nSyn,1:nSyn],[synIdx;synIdx(idxBalSyn)],[-pShare(idxBalSyn).*JG0;JG0(idxBalSyn).*pShare],nSyn,nbus);
    MatG5AD=sparse([1:nSyn,1:nSyn],[synIdx;synIdx(idxBalSyn)],[-pShare(idxBalSyn).*KG0;KG0(idxBalSyn).*pShare],nSyn,nbus);
    MatG5BC=sparse(nIsland,nbus);
    MatG5BD=sparse(nIsland,nbus);
    
    MatG5AC(idxBal,:)=MatG5BC;
    MatG5AD(idxBal,:)=MatG5BD;
    
    MatG1J=sparse(1:nSyn,1:nSyn,Rs.*Cd0+Xd.*Sd0,nSyn,nSyn);
    MatG1K=sparse(1:nSyn,1:nSyn,Rs.*Sd0-Xd.*Cd0,nSyn,nSyn);
    MatG1d=sparse(1:nSyn,1:nSyn,(CG0+Rs.*JG0-Xd.*KG0).*A1n+(DG0+Rs.*KG0+Xd.*JG0).*B1n,nSyn,nSyn);
    MatG2J=sparse(1:nSyn,1:nSyn,Rs.*Sd0-Xq.*Cd0,nSyn,nSyn);
    MatG2K=sparse(1:nSyn,1:nSyn,-Rs.*Cd0-Xq.*Sd0,nSyn,nSyn);
    MatG2d=sparse(1:nSyn,1:nSyn,(-DG0-Rs.*KG0-Xq.*JG0).*A1n+(CG0+Rs.*JG0-Xq.*KG0).*B1n,nSyn,nSyn);
    MatG5AJ=sparse([1:nSyn,1:nSyn],[1:nSyn,idxBalSyn'],[-pShare(idxBalSyn).*(CG0+2*JG0.*Rs);(CG0(idxBalSyn)+2*JG0(idxBalSyn).*Rs(idxBalSyn)).*pShare],nSyn,nSyn);
    MatG5AK=sparse([1:nSyn,1:nSyn],[1:nSyn,idxBalSyn'],[-pShare(idxBalSyn).*(DG0+2*KG0.*Rs);(DG0(idxBalSyn)+2*KG0(idxBalSyn).*Rs(idxBalSyn)).*pShare],nSyn,nSyn);
    MatG5Ad=sparse(nSyn,nSyn);
    MatG5BJ=sparse(nIsland,nSyn);
    MatG5BK=sparse(nIsland,nSyn);
    MatG5Bd=sparse(synByIsland,1:nSyn,Ms,nIsland,nSyn);
    
    MatG5AJ(idxBal,:)=MatG5BJ;
    MatG5AK(idxBal,:)=MatG5BK;
    MatG5Ad(idxBal,:)=MatG5Bd;
    
    % [MatGA,MatGB]*[CD;JKd]=RHS
    % JKd=-MatGA\MatGB*CD+MatGA\RHS=MatGBiA*CD+MatGA\RHS
    MatGA=[MatG1C,MatG1D;MatG2C,MatG2D;MatG5AC,MatG5AD];
    MatGB=[MatG1J,MatG1K,MatG1d;MatG2J,MatG2K,MatG2d;MatG5AJ,MatG5AK,MatG5Ad];
    %         MatGB2=[MatG1J,MatG1K,sparse(1:nSyn,1:nSyn,(CG0+Rs.*JG0-Xd.*KG0),nSyn,nSyn),sparse(1:nSyn,1:nSyn,(DG0+Rs.*KG0+Xd.*JG0),nSyn,nSyn),sparse(nSyn,nSyn);...
    %             MatG2J,MatG2K,sparse(1:nSyn,1:nSyn,(-DG0-Rs.*KG0-Xq.*JG0),nSyn,nSyn),sparse(1:nSyn,1:nSyn,(CG0+Rs.*JG0-Xq.*KG0),nSyn,nSyn),sparse(nSyn,nSyn);...
    %             sparse(nSyn,nSyn),sparse(nSyn,nSyn),sparse(1:nSyn,1:nSyn,ones(1,nSyn),nSyn,nSyn),sparse(nSyn,nSyn),sparse(1:nSyn,1:nSyn,-A1n,nSyn,nSyn);...
    %             sparse(nSyn,nSyn),sparse(nSyn,nSyn),sparse(nSyn,nSyn),sparse(1:nSyn,1:nSyn,ones(1,nSyn),nSyn,nSyn),sparse(1:nSyn,1:nSyn,-B1n,nSyn,nSyn);...
    %             MatG5AJ,MatG5AK,sparse(nSyn,nSyn),sparse(nSyn,nSyn),MatG5Ad];
    
    MatGBiA=-MatGB\MatGA;
    GTrMat=sparse(synIdx,1:nSyn,ones(nSyn,1),nbus,nSyn);
    MatGTrans=[GTrMat,sparse(nbus,nSyn),sparse(nbus,nSyn);sparse(nbus,nSyn),GTrMat,sparse(nbus,nSyn)];
    MatGCD=MatGTrans*MatGBiA;
else
    d=[];
    JG=[];
    KG=[];
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
%     YLHS_1=[Y11,Y12;Y21,Y22];

if ~isempty(ind)
    YLHS=YLHS-...
        [sparse(1:nbus,1:nbus,LHS_MatInd_Bus(:,1,1),nbus,nbus),sparse(1:nbus,1:nbus,LHS_MatInd_Bus(:,1,2),nbus,nbus);...
        sparse(1:nbus,1:nbus,LHS_MatInd_Bus(:,2,1),nbus,nbus),sparse(1:nbus,1:nbus,LHS_MatInd_Bus(:,2,2),nbus,nbus)];
    %         YLHS_1=YLHS_1-...
    %             [sparse(1:nbus,1:nbus,LHS_MatInd_Bus_1(:,1,1),nbus,nbus),sparse(1:nbus,1:nbus,LHS_MatInd_Bus_1(:,1,2),nbus,nbus);...
    %             sparse(1:nbus,1:nbus,LHS_MatInd_Bus_1(:,2,1),nbus,nbus),sparse(1:nbus,1:nbus,LHS_MatInd_Bus_1(:,2,2),nbus,nbus)];
end

if ~isempty(syn)
    YLHS=YLHS+MatGCD;
    %         YLHS_1=YLHS_1+MatGCD;
end

% Branch 0: smaller matrix
C0i=C0./V0sq;
D0i=D0./V0sq;
% C0Mi=sparse(1:nbus,1:nbus,C0./V0sq,nbus,nbus);
% D0Mi=sparse(1:nbus,1:nbus,D0./V0sq,nbus,nbus);
PCQD=Pspec.*C0i+Qspec.*D0i;
PDQC=Pspec.*D0i-Qspec.*C0i;
YLHSx=YLHS-...
        [sparse(1:nbus,1:nbus,PCQD.*E0+PDQC.*F0,nbus,nbus),sparse(1:nbus,1:nbus,-PCQD.*F0+PDQC.*E0,nbus,nbus);...
        sparse(1:nbus,1:nbus, PDQC.*E0-PCQD.*F0,nbus,nbus),sparse(1:nbus,1:nbus,-PDQC.*F0-PCQD.*E0,nbus,nbus)];
idxNonSw=find(busType~=2);
idxStackMat=[idxNonSw;idxNonSw+nbus];
C0Mx=sparse(1:npv,ipvInPvPq,C0(ipv),npv,npv+npq);
D0Mx=sparse(1:npv,ipvInPvPq,D0(ipv),npv,npv+npq);

LHS_matx=[YLHSx(idxStackMat,idxStackMat),[-F0M(busType~=2,ipv);-E0M(busType~=2,ipv)];...
    C0Mx,D0Mx,sparse(npv,npv)];

useLU=isfield(SysPara,'iden')&&isfield(SysPara,'p_amd');

if useLU
    if isempty(SysPara.p_amd)
        p_amd = colamd (LHS_matx) ;
        save([SysPara.iden,'.mat'],'p_amd');
    else
        p_amd=SysPara.p_amd;
    end
    MxI = speye (size(LHS_matx)) ;
    MxQ = MxI (:, p_amd) ;
    if IS_OCTAVE
    [MxL,MxU,MxP,MxQx] = lu (LHS_matx*MxQ) ;
    else
    [MxL,MxU,MxP] = lu (LHS_matx*MxQ) ;
    end
end
% END OF Branch 0: smaller matrix

% Branch 1: larger matrix
% idxNonSw=find(busType~=2);
% idxStackMat=[idxNonSw;idxNonSw+nbus];
% 
% LHS_mat=[YLHS(idxStackMat,idxStackMat),...
%     [YEF11(busType~=2,busType~=2),YEF12(busType~=2,busType~=2),-F0M(busType~=2,ipv);...
%     YEF21(busType~=2,busType~=2),YEF22(busType~=2,busType~=2),-E0M(busType~=2,ipv)];...
%     C0M(ipv,busType~=2),D0M(ipv,busType~=2),sparse(npv,2*npq+3*npv);...
%     E0M(busType~=2,busType~=2),-F0M(busType~=2,busType~=2),C0M(busType~=2,busType~=2),-D0M(busType~=2,busType~=2),sparse(npq+npv,npv);...
%     F0M(busType~=2,busType~=2),E0M(busType~=2,busType~=2),D0M(busType~=2,busType~=2),C0M(busType~=2,busType~=2),sparse(npq+npv,npv);];
% 
% useLU=isfield(SysPara,'iden')&&isfield(SysPara,'p_amd');
% 
% if useLU
%     if isempty(SysPara.p_amd)
%         p_amd = colamd (LHS_mat) ;
%         save([SysPara.iden,'.mat'],'p_amd');
%     else
%         p_amd=SysPara.p_amd;
%     end
%     MxI = speye (size(LHS_mat)) ;
%     MxQ = MxI (:, p_amd) ;
%     [MxL,MxU,MxP] = lu (LHS_mat*MxQ) ;
% end
% END OF Branch 1: larger matrix

% if nbus<=LU_SIZE_TH
%     [L_LHS_mat,U_LHS_mat,p_LHS_mat]=lu(LHS_mat,'vector');
% end

%     LHS_mat_1=[YLHS_1,[YEF11,YEF12,-F0M(:,ipv);YEF21,YEF22,-E0M(:,ipv)];...
%         C0M(ipv,:),D0M(ipv,:),sparse(npv,2*nbus+npv);...
%         E0M,-F0M,C0M,-D0M,sparse(nbus,npv);...
%         F0M,E0M,D0M,C0M,sparse(nbus,npv);];

for i=1:nlvl
    
    seq2=getseq(i,2);
    seq2p=getseq(i+1,2);
    seq3=getseq(i,3);
    idxSeq2=sum(seq2==i,2);
    idxSeq2p=sum(seq2p>=i,2);
    idxSeq3=sum(seq3==i,2);
    idxSeq3x=sum(seq3(:,[2,3])==i,2);
    seq2R=seq2(idxSeq2==0,:);
    seq2px=seq2(idxSeq2p==0,:);
    seq3R=seq3(idxSeq3==0,:);
    seq3Rx=seq3(idxSeq3x==0,:);
    seq2m=getseq(i-1,2);
    seq3m=getseq(i-1,3);
    seq2mm=getseq(i-2,2);
    
    RHSILr=zeros(nbus,1);
    RHSILi=zeros(nbus,1);
    if ~isempty(ind)
        %             rhsBus=zeros(2,nInd);
        %             rhsM=sum(Vm(:,seq2R(:,1)+1).*s(:,seq2R(:,2)+1),2)-1j*X2.*sum(IR(:,seq2R(:,1)+1).*s(:,seq2R(:,2)+1),2);
        %             rhsI=-real(sum(IR(:,seq2R(:,1)+1).*conj(IR(:,seq2R(:,2)+1)),2))+...
        %                 (T1.*sum(s(:,seq2R(:,1)+1).*s(:,seq2R(:,2)+1),2)+T2.*sum(s(:,seq3R(:,1)+1).*s(:,seq3R(:,2)+1).*s(:,seq3R(:,3)+1),2)).*Rind0./R2+...
        %                 (T0.*s(:,i)+T1.*sum(s(:,seq2m(:,1)+1).*s(:,seq2m(:,2)+1),2)+T2.*sum(s(:,seq3m(:,1)+1).*s(:,seq3m(:,2)+1).*s(:,seq3m(:,3)+1),2)).*Rind1./R2;
        %             rhsIx=-real(sum(IR(:,seq2px(:,1)+1).*conj(IR(:,seq2px(:,2)+1)),2))+...
        %                 (T1.*sum(s(:,seq2R(:,1)+1).*s(:,seq2R(:,2)+1),2)+T2.*sum(s(:,seq3R(:,1)+1).*s(:,seq3R(:,2)+1).*s(:,seq3R(:,3)+1),2)).*Rind1./R2;
        %             rhsImod=Rind1.*(T1.*s(:,i)+T2.*sum(s(:,seq2m(:,1)+1).*s(:,seq2m(:,2)+1),2))+Rind0.*T2.*sum(s(:,seq2R(:,1)+1).*s(:,seq2R(:,2)+1),2)-...
        %                 real(sum(V(indIdx,seq2R(:,1)+1).*conj(IR(:,seq2R(:,2)+1)),2))+...
        %                 real(sum(IL(:,seq2R(:,1)+1).*conj(IR(:,seq2R(:,2)+1)),2).*Z1);
        %             if i==1
        %                 rhsImod=rhsImod+Rind1.*T0;
        %             end
        %             rhsIL=V(indIdx,i).*Yeind1-IL(:,i).*Ye1ind1;
        %             for j=1:nInd
        %                 %                 if Rind0(j)~=0
        %                 %                     rhsBus(:,j)=squeeze(RHS_C_Shr(j,:,:))*[real(rhsM(j));imag(rhsM(j));rhsI(j);real(rhsIL(j));imag(rhsIL(j))];
        %                 %                 else
        %                 %                     if i==1
        %                 %                         s(j,i+1)=R2(j)*Rind1(j)*T0(j)/...
        %                 %                             ((C0(indIdx(j))-R1(j)*JL0(j)+X1(j)*KL0(j))*(C0(indIdx(j))-R1(j)*JL0(j)+X1(j)*KL0(j))+...
        %                 %                             (D0(indIdx(j))-X1(j)*JL0(j)-R1(j)*KL0(j))*(D0(indIdx(j))-X1(j)*JL0(j)-R1(j)*KL0(j)));
        %                 %                         IR(j,i+1)=s(j,i+1)/R2(j)*((C0(indIdx(j))-R1(j)*JL0(j)+X1(j)*KL0(j))+1j*(D0(indIdx(j))-X1(j)*JL0(j)-R1(j)*KL0(j)));
        %                 %                         rhsBus(:,j)=squeeze(RHS_C_Shr_1(j,:,:))*[real(IR(j,i+1)+rhsIL(j));imag(IR(j,i+1)+rhsIL(j));0;0;0];
        %                 %                     else
        %                 %                         rhsBus(:,j)=squeeze(RHS_C_Shr(j,:,:))*[real(rhsM(j));imag(rhsM(j));rhsIx(j);real(rhsIL(j));imag(rhsIL(j))];
        %                 %                     end
        %                 %                 end
        %                 rhsBus(:,j)=squeeze(RHS_C_Shr(j,:,:))*[real(rhsM(j));imag(rhsM(j));rhsImod(j);real(rhsIL(j));imag(rhsIL(j))];
        %             end
        %             %accumulate currents
        %             RHSILr=accumarray(indIdx,rhsBus(1,:)',[nbus,1]);
        %             RHSILi=accumarray(indIdx,rhsBus(2,:)',[nbus,1]);
        
        rhsBus=zeros(5,nInd);
        rhsM=sum(Vm(:,seq2R(:,1)+1).*s(:,seq2R(:,2)+1),2)-1j*X2.*sum(IR(:,seq2R(:,1)+1).*s(:,seq2R(:,2)+1),2);
        rhsImod=Rind1.*(T1.*s(:,i)+T2.*sum(s(:,seq2m(:,1)+1).*s(:,seq2m(:,2)+1),2))+Rind0.*T2.*sum(s(:,seq2R(:,1)+1).*s(:,seq2R(:,2)+1),2)-...
            real(sum(V(indIdx,seq2R(:,1)+1).*conj(IR(:,seq2R(:,2)+1)),2))+...
            real(sum(IL(:,seq2R(:,1)+1).*conj(IR(:,seq2R(:,2)+1)),2).*Z1);
        if i==1
            rhsImod=rhsImod+Rind1.*T0;
        end
        rhsIL=V(indIdx,i).*Yeind1-IL(:,i).*Ye1ind1;
        for j=1:nInd
            rhsBus(:,j)=squeeze(RHS_C_Shr(j,:,:))*[real(rhsM(j));imag(rhsM(j));rhsImod(j);real(rhsIL(j));imag(rhsIL(j))];
        end
        RHSILr=accumarray(indIdx,rhsBus(3,:)',[nbus,1]);
        RHSILi=accumarray(indIdx,rhsBus(4,:)',[nbus,1]);
    end
    
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
        AG0=zeros(nSyn,1);
        BG0=zeros(nSyn,1);
        
        if taylorN>=2
            AG0=AG0+cosp(:,3).*sum(d(:,seq2R(:,1)+1).*d(:,seq2R(:,2)+1),2);
            BG0=BG0+sinp(:,3).*sum(d(:,seq2R(:,1)+1).*d(:,seq2R(:,2)+1),2);
        end
        if taylorN>=3
            AG0=AG0+cosp(:,4).*sum(d(:,seq3R(:,1)+1).*d(:,seq3R(:,2)+1).*d(:,seq3R(:,3)+1),2);
            BG0=BG0+sinp(:,4).*sum(d(:,seq3R(:,1)+1).*d(:,seq3R(:,2)+1).*d(:,seq3R(:,3)+1),2);
        end
        if taylorN>=4
            seq4=getseq(i,4);
            idxSeq4=sum(seq4==i,2);
            seq4R=seq4(idxSeq4==0,:);
            
            AG0=AG0+cosp(:,5).*sum(d(:,seq4R(:,1)+1).*d(:,seq4R(:,2)+1).*d(:,seq4R(:,3)+1).*d(:,seq4R(:,4)+1),2);
            BG0=BG0+sinp(:,5).*sum(d(:,seq4R(:,1)+1).*d(:,seq4R(:,2)+1).*d(:,seq4R(:,3)+1).*d(:,seq4R(:,4)+1),2);
        end
        
        CCr=sum(real(V(synIdx,seq2R(:,1)+1)).*Cd(:,seq2R(:,2)+1),2);
        DCr=sum(imag(V(synIdx,seq2R(:,1)+1)).*Cd(:,seq2R(:,2)+1),2);
        CSr=sum(real(V(synIdx,seq2R(:,1)+1)).*Sd(:,seq2R(:,2)+1),2);
        DSr=sum(imag(V(synIdx,seq2R(:,1)+1)).*Sd(:,seq2R(:,2)+1),2);
        JCr=sum(JG(:,seq2R(:,1)+1).*Cd(:,seq2R(:,2)+1),2);
        KCr=sum(KG(:,seq2R(:,1)+1).*Cd(:,seq2R(:,2)+1),2);
        JSr=sum(JG(:,seq2R(:,1)+1).*Sd(:,seq2R(:,2)+1),2);
        KSr=sum(KG(:,seq2R(:,1)+1).*Sd(:,seq2R(:,2)+1),2);
        
        RHSIG1=Ef(:,i+1)-(CCr+DSr+Rs.*(JCr+KSr)+Xd.*(JSr-KCr))-(CG0+Rs.*JG0-Xd.*KG0).*AG0-(DG0+Rs.*KG0+Xd.*JG0).*BG0;
        RHSIG2=-(CSr-DCr+Rs.*(JSr-KCr)-Xq.*(JCr+KSr))-(-DG0-Rs.*KG0-Xq.*JG0).*AG0-(CG0+Rs.*JG0-Xq.*KG0).*BG0;
        RHSIG3temp=-Pm(:,i+1)+...
            sum(real(V(synIdx,seq2R(:,1)+1).*(JG(:,seq2R(:,2)+1)-1j*KG(:,seq2R(:,2)+1))),2)+...
            (sum(JG(:,seq2R(:,1)+1).*JG(:,seq2R(:,2)+1),2)+sum(KG(:,seq2R(:,1)+1).*KG(:,seq2R(:,2)+1),2)).*Rs;
        
        RHSIG3=pShare(idxBalSyn).*RHSIG3temp-RHSIG3temp(idxBalSyn).*pShare;
        RHSIG3(idxBal)=0;
        RHSIG=MatGB\[RHSIG1;RHSIG2;RHSIG3];
        RHSIGJK=MatGTrans*RHSIG;
        RHSIGr=RHSIGJK(1:nbus);
        RHSIGi=RHSIGJK((nbus+1):(2*nbus));
    end
    
    % HEM Body
    RHS1=sum((-P(:,seq2(:,1)+1)+1j*(Q(:,seq2(:,1)+1)+Qxtra(:,seq2(:,1)+1))).*conj(W(:,seq2(:,2)+1)),2)+Ysh1.*V(:,i);
    RHS2=-0.5*real(sum(V(:,seq2R(:,1)+1).*conj(V(:,seq2R(:,2)+1)),2));
    RHS3=sum(-W(:,seq2R(:,1)+1).*V(:,seq2R(:,2)+1),2);
    
    if i==1
        RHS2=RHS2+0.5*VspSq2(:,2);
    end
    
    compactRHS1=RHS1(busType~=2);
    if ~isempty(isw);compactRHS1=compactRHS1+Y(busType~=2,isw)*V(isw,i+1);end
    
    RHS3r=real(RHS3);RHS3i=imag(RHS3);
    
%     RHS=[real(compactRHS1)+RHSILr(busType~=2)+RHSIiLr(busType~=2)-RHSIGr(busType~=2);...
%         imag(compactRHS1)+RHSILi(busType~=2)+RHSIiLi(busType~=2)-RHSIGi(busType~=2);...
%         RHS2(ipv);...
%         RHS3r(busType~=2);...
%         RHS3i(busType~=2);];  
%     
%     if useLU
%         x = MxQ * (MxU \ (MxL \ (MxP * RHS))) ;
%     else
%         x=LHS_mat\RHS;
%     end
%     
%     V(busType~=2,i+1)=x(1:(npq+npv))+1j*x(((npq+npv)+1):(2*(npq+npv)));
%     W(busType~=2,i+1)=x((2*(npq+npv)+1):(3*(npq+npv)))+1j*x((3*(npq+npv)+1):(4*(npq+npv)));
%     Q(ipv,i+1)=x((4*(npq+npv)+1):end);
    
    % Solve a smaller system
    
    RHS3xr=PCQD.*RHS3r+PDQC.*RHS3i;
    RHS3xi=PDQC.*RHS3r-PCQD.*RHS3i;
    RHSx=[real(compactRHS1)+RHSILr(busType~=2)+RHSIiLr(busType~=2)-RHSIGr(busType~=2)-RHS3xr(busType~=2);...
        imag(compactRHS1)+RHSILi(busType~=2)+RHSIiLi(busType~=2)-RHSIGi(busType~=2)-RHS3xi(busType~=2);...
        RHS2(ipv);];
    
    if useLU
      if IS_OCTAVE
        x = MxQ * MxQx* (MxU \ (MxL \ (MxP * RHSx))) ;
      else
        x = MxQ * (MxU \ (MxL \ (MxP * RHSx))) ;
      end
    else
        x=LHS_matx\RHSx;
    end
    
    V(busType~=2,i+1)=x(1:(npq+npv))+1j*x(((npq+npv)+1):(2*(npq+npv)));
    Q(ipv,i+1)=x((2*(npq+npv)+1):end);
    
    Cx=real(V(:,i+1));Dx=imag(V(:,i+1));
    RHS3xxr=RHS3r-E0.*Cx+F0.*Dx;RHS3xxi=RHS3i-F0.*Cx-E0.*Dx;
    Ex=C0i.*RHS3xxr+D0i.*RHS3xxi;Fx=-D0i.*RHS3xxr+C0i.*RHS3xxi;
    W(busType~=2,i+1)=Ex(busType~=2)+1j*Fx(busType~=2);
    
    if ~isempty(ind)
        for j=1:nInd
            %                 if Rind0(j)~=0
            %                     tempIL=squeeze(LHS_MatInd_Shr(j,:,:))*[real(V(indIdx(j),i+1));imag(V(indIdx(j),i+1))]+rhsBus(:,j);
            %                     tempIRs=-squeeze(LHS_MatInd_Shr2(j,:,:))*[tempIL;real(V(indIdx(j),i+1));imag(V(indIdx(j),i+1))]+...
            %                         squeeze(LHS_MatInd_Shr3(j,:,:))*[rhsI(j);real(rhsIL(j));imag(rhsIL(j))];
            %                     IL(j,i+1)=tempIL(1)+1j*tempIL(2);
            %                     IR(j,i+1)=tempIRs(1)+1j*tempIRs(2);
            %                     s(j,i+1)=tempIRs(3);
            %                     Vm(j,i+1)=V(indIdx(j),i+1)-IL(j,i+1)*Z1(j);
            %                 else % singular case
            %                     if i==1
            %                         % s and IR already generated.
            %                         %                         s(j,i+1)=R2(j)*Rind1(j)*T0(j)/...
            %                         %                             ((C0(indIdx(j))-R1(j)*JL0(j)+X1(j)*KL0(j))*(C0(indIdx(j))-R1(j)*JL0(j)+X1(j)*KL0(j))+...
            %                         %                             (D0(indIdx(j))-X1(j)*JL0(j)-R1(j)*KL0(j))*(D0(indIdx(j))-X1(j)*JL0(j)-R1(j)*KL0(j)));
            %                         %                         IR(j,i+1)=s(j,i+1)/R2(j)*((C0(indIdx(j))-R1(j)*JL0(j)+X1(j)*KL0(j))+1j*(D0(indIdx(j))-X1(j)*JL0(j)-R1(j)*KL0(j)));
            %                         IL(j,i+1)=(V(indIdx(j),i+1)*Yeind0(j)+IR(j,i+1)+V(indIdx(j),i)*Yeind1(j)-IL(j,i)*Ye1ind1(j))/(1+Ye1ind0(j));
            %                         Vm(j,i+1)=V(indIdx(j),i+1)-IL(j,i+1)*Z1(j);
            %                     else
            %                         tempIL=squeeze(LHS_MatInd_Shr(j,:,:))*[real(V(indIdx(j),i+1));imag(V(indIdx(j),i+1))]+rhsBus(:,j);
            %                         tempIRs=-squeeze(LHS_MatInd_Shr2(j,:,:))*[tempIL;real(V(indIdx(j),i+1));imag(V(indIdx(j),i+1))]+...
            %                             squeeze(LHS_MatInd_Shr3(j,:,:))*[rhsIx(j);real(rhsIL(j));imag(rhsIL(j))];
            %                         IL(j,i+1)=tempIL(1)+1j*tempIL(2);
            %                         IR(j,i+1)=tempIRs(1)+1j*tempIRs(2);
            %                         s(j,i+1)=tempIRs(3);
            %                         Vm(j,i+1)=V(indIdx(j),i+1)-IL(j,i+1)*Z1(j);
            %                     end
            %                 end
            tempx=squeeze(LHS_MatInd_Full(j,:,:))*[real(V(indIdx(j),i+1));imag(V(indIdx(j),i+1))]+rhsBus(:,j);
            IL(j,i+1)=tempx(3)+1j*tempx(4);
            IR(j,i+1)=tempx(1)+1j*tempx(2);
            s(j,i+1)=tempx(5);
            Vm(j,i+1)=V(indIdx(j),i+1)-IL(j,i+1)*Z1(j);
        end
    end
    
    if ~isempty(zip)
        IiL(:,i+1)=(LHS_MatZip(:,1)+1j*LHS_MatZip(:,3)).*real(V(zipIdx,i+1))+(LHS_MatZip(:,2)+1j*LHS_MatZip(:,4)).*imag(V(zipIdx,i+1))+(RHSILr_full+1j*RHSILi_full);
        BiL(:,i+1)=Mat_BZip(:,1).*real(V(zipIdx,i+1))+Mat_BZip(:,2).*imag(V(zipIdx,i+1))+RHS_BZip;
    end
    
    if ~isempty(syn)
        IGJKd=MatGBiA*[real(V(:,i+1));imag(V(:,i+1))]+RHSIG;
        JG(:,i+1)=IGJKd(1:nSyn);
        KG(:,i+1)=IGJKd((nSyn+1):(2*nSyn));
        d(:,i+1)=IGJKd((2*nSyn+1):(3*nSyn));
        Cd(:,i+1)=A1n.*d(:,i+1)+AG0;
        Sd(:,i+1)=B1n.*d(:,i+1)+BG0;
    end
end
end