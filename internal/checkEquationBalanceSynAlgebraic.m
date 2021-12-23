function [diff,psw,qsw]=checkEquationBalanceSynAlgebraic(SysData,SysPara,x)
% Checking the equation balance in solving algebraic equations
% FUNCTION checkEquationBalanceSynAlgebraic
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
% OUTPUT
%	diff - mismatch vector of the equations
%	psw - active power of SW
%	qsw - reactive power of SW
%
[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
nbus=size(bus,1);
nline=size(line,1);
[V,Q,s,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,f,dpg,qplt,vg]=unfoldX(x,SysData);
[pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,Vsp2,MatGV,MatGV1,MatGRhs,MatGRhs1]=unfoldSysPara(SysPara);

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

isw=find(busType==2);
ipv=find(busType~=0);
ipq=find(busType==0);
npq=size(ipq,1);
npv=size(ipv,1);

pVec=zeros(nbus,1);
qVec=zeros(nbus,1);

% Determine the frequency model of each island
freqTypeTag=zeros(nIslands,1);%0:sw,1:syn,2:steady-state f
freqKeptTag=zeros(nbus,1);
frefs=refs;
fswTag=zeros(nbus,1);
fsynTag=zeros(nbus,1);
fswTag(isw)=1;
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
if ~isempty(zip)%zipMode=0
    pVec=pVec-accumarray(zip(:,1),zip(:,7).*zip(:,12),[nbus,1]);
    qVec=qVec-accumarray(zip(:,1),zip(:,10).*zip(:,12),[nbus,1]);
end
% qVec(pv(:,1))=0;
pVec=pVec+freqKeptTag.*(-fdk.*f+dpg);

IInj=Y*V;
if ~isempty(zip)%zipMode=0
    IInj=IInj+accumarray(zip(:,1),(zip(:,6)-1j*zip(:,9)).*zip(:,12).*V(zip(:,1))./abs(V(zip(:,1))),[nbus,1]);
end
SInjRHS=V.*conj(IInj);

if ~isempty(syn)    
    nSyn=size(syn,1);
    synIdx=syn(:,1);    
    Cg=real(V(synIdx));
    Dg=imag(V(synIdx));
    IGx=-MatGV(:,1).*Cg-MatGV(:,2).*Dg+MatGRhs(:,1);
    IGy=-MatGV(:,3).*Cg-MatGV(:,4).*Dg+MatGRhs(:,2);
    IG=IGx+1j*IGy;
%     IG=(Ef.*(cos(d)+1j*sin(d))-V(synIdx))./(Rs+1j*Xd);
    SInjRHS_syn=-V.*conj(accumarray(synIdx,IG,[nbus,1]));
    SInjRHS=SInjRHS+SInjRHS_syn;
end

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

% qVec(busType~=0)=0;
diffSInj=SInjRHS-(pVec+1j*qVec);
diffSInj=diffSInj-1j*Q;
ssw=diffSInj(busType==2);
psw=real(ssw);qsw=imag(ssw);

diffSInj(busType==2)=0;

diffv=sum(Vsp2,2)-V.*conj(V);

diff=[real(diffSInj);imag(diffSInj);diffv(busType~=0)];
end