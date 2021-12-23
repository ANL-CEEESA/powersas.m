function [Vx,Wx,Qx,fx,xx,finalAlpha,alphas,diffRec]=hemMachinePFmultiStageAlgebraic(SimData,SysData,SysPara,x0,xNew,SysDataNew)
% Multi-stage scheme for solving algebraic equations
%
% FUNCTION hemMachinePFmultiStageAlgebraic
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
%   x0 - Initial system state
%   xNew - Updated system state (algebraic variables to be solved)
%   SysDataNew - (optional) updated system data
%	
% OUTPUT
%	Vx, Wx, Qx, fx - Last HE coefficients of algebraic variables
%	xx - Solved updated system state
%   finalAlpha - The ending length of this segment of simulation
%   alphas - Record of alphas
%   diffRec - A list of errors
%

setting_func_hemMachinePFmultiStageAlgebraic;

[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
if nargin>=6&&~isempty(SysDataNew)
    [~,~,~,~,~,~,~,~,synNew,~,~]=unfoldSysData(SysDataNew);
    deltaSyn=synNew(:,[7:10,13:15])-syn(:,[7:10,13:15]);
else
    synNew=[];
    deltaSyn=zeros(size(syn,1),7);
end
nbus=size(bus,1);
nline=size(line,1);
nSyn=size(syn,1);

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

[V0,Q0,s0,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,tgovmech0,f0,dpg0,qplt0,vg0]=unfoldX(x0,SysData);
eq0=Ef0;ed0=zeros(size(Ef0));
[~,~,sNew,dNew,~,eq1New,eq2New,ed1New,ed2New,psidNew,psiqNew,~,EfNew,~,~,~,~,~,~,~,fNew,dpgNew,qpltNew,vgNew]=unfoldX(xNew,SysData);
eqNew=EfNew;edNew=zeros(size(EfNew));
[maxAlpha,segAlpha,dAlpha,nlvl,taylorN,alphaTol,diffTol,diffTolMax,method]=unfoldSimData(SimData);
[pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,MatGV0,MatGV1,MatGRhs0,MatGRhs1]=unfoldSysPara(SysPara);
%
[MatGV0,MatGV1,MatGRhs0,MatGRhs1]=getLinearInterpolatorSyn(syn,synNew,d0,dNew,ed0,ed10,ed20,edNew,ed1New,ed2New,eq0,eq10,eq20,eqNew,eq1New,eq2New,psid0,psiq0,psidNew,psiqNew);

Ysh0Orig=Ysh0;
Ysh1Orig=Ysh1;
if ~isempty(ind)
    [YshInd0,Yshind1]=getLinearInterpolatorInd(nbus,ind,s0,sNew);
    Ysh0=Ysh0+YshInd0;
    Ysh1=Ysh1+Yshind1;
end

[nIslands,islands,refs]=searchIslands(bus(:,1),line(:,1:2));
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

if ~isempty(agc)
    agcExt=zeros(nbus,size(agc,2));
    agcExt(agc(:,1),:)=agc;
    fdk=agcExt(:,2)+agcExt(:,3); %1/R+D
else
    fdk=zeros(nbus,1);
end

if ~isempty(pv)
    extraPV=freqKeptTag.*(dpgNew-dpg0);
    pvIncr=pvIncr+extraPV(pv(:,1));
end

pqx=pq;
pvx=pv;
shuntx=shunt;
indx=ind;
zipx=zip;
synx=syn;

V0x=V0;
Q0x=Q0;
f0x=f0;
Rzip0x=Rzip0;
Ysh0x=Ysh0;
Ysh0Origx=Ysh0Orig;
Ytr0x=Ytr0;
Vsp0x=real(V0.*conj(V0));
MatGV0x=MatGV0;
MatGRhs0x=MatGRhs0;

alphaConfirm=0;

busType=zeros(nbus,1);
if ~isempty(pv);busType(pv(:,1))=1;end
if ~isempty(sw);busType(sw(:,1))=2;end

vc=zeros(nbus,1);
qc=zeros(size(Q0,1),1);

alphas=0;
skip=0;
heCounter=0;
noMove=0;

while alphaConfirm<1
    alpha=1-alphaConfirm;
    alphax=alpha/alphaPreMult;
    if skip==0
        SysDatax=foldSysData(bus,sw,pvx,pqx,shuntx,line,indx,zipx,synx,exc,tg,agc,cac,cluster);
        xx=foldX(SysDatax,V0x,Q0x,s0,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,tgovmech0,f0x,dpg0,qplt0,vg0);
        SysParax=foldSysPara(pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0x,Rzip1,Ytr0x,Ytr1,Ysh0x,Ysh1,[Vsp0x,VspSq2(:,2)],MatGV0x,MatGV1,MatGRhs0x,MatGRhs1);
        
        if isfield(SysPara,'iden')&&isfield(SysPara,'p_amd')
            SysParax.iden=SysPara.iden;
            SysParax.p_amd=SysPara.p_amd;
            if isempty(SysPara.p_amd)
                if exist([SysPara.iden,'.mat'],'file')
                    if ~exist('p_amd','var')
                        load([SysPara.iden,'.mat']);
                    end
                    SysParax.p_amd=p_amd;
                end
            end
        end
        [V,W,Q,f]=hemMachinePFSalientcontinueAlgebraic(SimData,SysDatax,SysParax,xx);
        heCounter=heCounter+1;
                
        [Va3,Vb3]=getPadeApproxHelm(V,0);
        [Qa3,Qb3]=getPadeApproxHelm(Q,find(busType==0));
        [fa3,fb3]=getPadeApproxHelm(f,0);
        idxNanf=find(isnan(fa3(:,1))|isnan(fb3(:,1)));fa3(idxNanf,:)=0;fb3(idxNanf,:)=0;fa3(idxNanf,1:2)=f(idxNanf,1:2);
    end
    skip=0;
    alphaLeft=0;
    alphaRight=alpha/alphaPreMult;
    
    while 1
        vc=polyvalVec(Va3(:,end:-1:1),alphax)./polyvalVec([Vb3(:,end:-1:1),ones(size(Vb3,1),1)],alphax);
        vc(isnan(vc))=V(isnan(vc),1)+alphax*V(isnan(vc),2);
        qc=polyvalVec(Qa3(:,end:-1:1),alphax)./polyvalVec([Qb3(:,end:-1:1),ones(size(Qb3,1),1)],alphax);
        qc(isnan(qc))=Q(isnan(qc),1)+alphax*Q(isnan(qc),2);
        fc=polyvalVec(fa3(:,end:-1:1),alphax)./polyvalVec([fb3(:,end:-1:1),ones(size(fb3,1),1)],alphax);
        
        Y=Ytr0x+alphax*Ytr1+sparse(1:nbus,1:nbus,Ysh0x+Ysh1*alphax,nbus,nbus);
        Vsp2=Vsp0x+alphax*VspSq2(:,2);
        
        pqxx=pqx;
        pqxx(:,[4,5])=pqx(:,[4,5])+alphax*pqIncr;
        pvxx=pvx;
        if ~isempty(pvx);pvxx(:,4)=pvx(:,4)+alphax*pvIncr;end
        
        zipxx=zipx;
        if ~isempty(zipxx)
            zipxx(:,5:10)=repmat(Rzip0x+alphax*Rzip1,1,6).*zipx(:,5:10);
        end
        
        synxx=synx;
        if ~isempty(synxx)
            synxx(:,[7:10,13:15])=synx(:,[7:10,13:15])+alphax*deltaSyn;
        end
        
        MatGVxx=MatGV0x+alphax*MatGV1;
        MatGRhsxx=MatGRhs0x+alphax*MatGRhs1;
        
        SysDataxx=foldSysData(bus,sw,pvxx,pqxx,shunt,line,indx,zipxx,synxx,exc,tg,agc,cac,cluster);
        xxx=foldX(SysDataxx,vc,qc,s0,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,tgovmech0,fc,dpg0,qplt0,vg0);
        SysParaxx=foldSysPara([],[],[],[],[],[],[],[],Ytr0x+alphax*Ytr1,[],Ysh0Origx+Ysh1Orig*alphax,[],[Vsp2,zeros(size(Vsp2))],MatGVxx,[],MatGRhsxx,[]);
        
        diff=checkEquationBalanceSynAlgebraic(SysDataxx,SysParaxx,xxx);
                
        if heCounter<=1
            Vx=V;Wx=W;Qx=Q;fx=f;
        end
        if max(abs(diff))<diffTol&&alphax==alpha
            V0x=vc;
            Q0x=qc;
            f0x=fc;
            
            alphaLeft=alphaRight;
        else
            if max(abs(diff))<diffTol
                V0x=vc;
                Q0x=qc;
                f0x=fc;
                
                alphaLeft=alphax;
                alphax=(alphaLeft+alphaRight)/2;
            else
                alphaRight=alphax;
                alphax=(alphaLeft+alphaRight)/2;
            end
        end
        if alphaRight-alphaLeft<alphaTol
            break;
        end
    end
        
    diffParadigm=max(abs(diff))+1e-12;
    countCheckAlpha=1;
    alphax=alphaLeft*alphaPreMult;
    if alpha-alphax<alphaTol/10
        alphax=alpha;
    end
    
    while countCheckAlpha<maxCount
        vc=polyvalVec(Va3(:,end:-1:1),alphax)./polyvalVec([Vb3(:,end:-1:1),ones(size(Vb3,1),1)],alphax);
        vc(isnan(vc))=V(isnan(vc),1)+alphax*V(isnan(vc),2);
        qc=polyvalVec(Qa3(:,end:-1:1),alphax)./polyvalVec([Qb3(:,end:-1:1),ones(size(Qb3,1),1)],alphax);
        qc(isnan(qc))=Q(isnan(qc),1)+alphax*Q(isnan(qc),2);
        fc=polyvalVec(fa3(:,end:-1:1),alphax)./polyvalVec([fb3(:,end:-1:1),ones(size(fb3,1),1)],alphax);
        
        Y=Ytr0x+alphax*Ytr1+sparse(1:nbus,1:nbus,Ysh0x+Ysh1*alphax,nbus,nbus);
        Vsp2=Vsp0x+alphax*VspSq2(:,2);
        
        pqxx=pqx;
        pqxx(:,[4,5])=pqx(:,[4,5])+alphax*pqIncr;
        pvxx=pvx;
        if ~isempty(pvx);pvxx(:,4)=pvx(:,4)+alphax*pvIncr;end
        
        zipxx=zipx;
        if ~isempty(zipxx)
            zipxx(:,5:10)=repmat(Rzip0x+alphax*Rzip1,1,6).*zipx(:,5:10);
        end
        
        synxx=synx;
        if ~isempty(synxx)
            synxx(:,[7:10,13:15])=synx(:,[7:10,13:15])+alphax*deltaSyn;
        end
        
        MatGVxx=MatGV0x+alphax*MatGV1;
        MatGRhsxx=MatGRhs0x+alphax*MatGRhs1;
        
        SysDataxx=foldSysData(bus,sw,pvxx,pqxx,shunt,line,indx,zipxx,synxx,exc,tg,agc,cac,cluster);
        xxx=foldX(SysDataxx,vc,qc,s0,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,tgovmech0,fc,dpg0,qplt0,vg0);
        SysParaxx=foldSysPara([],[],[],[],[],[],[],[],Ytr0x+alphax*Ytr1,[],Ysh0Origx+Ysh1Orig*alphax,[],[Vsp2,zeros(size(Vsp2))],MatGVxx,[],MatGRhsxx,[]);
        
        diff=checkEquationBalanceSynAlgebraic(SysDataxx,SysParaxx,xxx); 
        if max(abs(diff))<1.05*diffParadigm
            break;
        end
        alphax=alphax*alphaPreMult;
        countCheckAlpha=countCheckAlpha+1;
    end
    
    if countCheckAlpha>=maxCount;alphax=alphax/alphaPreMult;end
    alpha=alphax;
    pqx(:,[4,5])=pqx(:,[4,5])+alpha*pqIncr;
    synx(:,[7:10,13:15])=synx(:,[7:10,13:15])+alpha*deltaSyn;
    if ~isempty(pvx);pvx(:,4)=pvx(:,4)+alpha*pvIncr;end
    Rzip0x=Rzip0x+alpha*Rzip1;
    Ysh0x=Ysh0x+alpha*Ysh1;
    Ysh0Origx=Ysh0Origx+alpha*Ysh1Orig;
    Vsp0x=Vsp0x+alpha*VspSq2(:,2);
    Ytr0x=Ytr0x+alpha*Ytr1;
    MatGV0x=MatGV0x+alpha*MatGV1;
    MatGRhs0x=MatGRhs0x+alpha*MatGRhs1;
    
    alphaConfirm=alphaConfirm+alpha;
    addLog(['Alpha=',num2str(alphaConfirm),', da=',num2str(alpha),', (maxdiff<',num2str(diffTol),').'],'INFO');
    
    if alpha==0
        addLog('Step did not move!','INFO');
        noMove=noMove+1;
        if noMove>=MaxNoMove            
            addLog('Reached consecutive max no move, exit!','INFO');
            break;
        end
        if diffTol>=diffTolMax
            addLog('Max DiffTol reached and not move, exit!','INFO');
            break;
        end
        
        diffTol=max([diffTol*diffMul,max(abs(diff))+1e-9]);
        if diffTol>diffTolMax
            diffTol=diffTolMax;
        end
        addLog(['Enlarge tol! (Tol=',num2str(diffTol),')'],'INFO');
        skip=0;
    else    
        noMove=0;
        
        V0x=vc;
        Q0x=qc;
        f0x=fc;
        
        alphas=[alphas;alpha];
    end
end

addLog(['Total Step=',num2str(alphaConfirm),'.'],'INFO');
VSol=V0x;
QSol=Q0x;
fSol=f0x;
dpgSol=dpg0;
if ~isempty(pv)
    dpgSol=dpgSol+extraPV*alphaConfirm;
end

xx=foldX(SysData,VSol,QSol,s0,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,tgovmech0,fSol,dpgSol,qplt0,vg0);
diffRec=max(abs(diff));
finalAlpha=alphaConfirm;
end


