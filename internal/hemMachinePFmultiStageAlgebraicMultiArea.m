function [Vx,Wx,Qx,VSol,QSol,finalAlpha,alphas,diffRec]=...
    hemMachinePFmultiStageAlgebraicMultiArea(SimData,SysData,SysPara,SysDataList,SysParaList,links,grossMaps,sysMaps,x0,xNew,SysDataNew)
% Multi-stage scheme for solving algebraic equations in multiple areas (parallel computation enabled)
%
% FUNCTION hemMachinePFmultiStageAlgebraicMultiArea
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
%   xNew - Updated system state (algebraic variables to be solved)
%   SysDataNew - (optional) updated system data
%	
% OUTPUT
%	Vx, Wx, Qx - Last HE coefficients of algebraic variables
%	VSol,QSol - Updated V and Q
%   finalAlpha - The ending length of this segment of simulation
%   alphas - Record of alphas
%   diffRec - A list of errors
%
%
[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
if nargin>=11&&~isempty(SysDataNew)
    [~,~,~,~,~,~,~,~,synNew,~,~]=unfoldSysData(SysDataNew);
    deltaSyn=synNew(:,[7:10,13:15])-syn(:,[7:10,13:15]);
else
    synNew=[];
    deltaSyn=zeros(size(syn,1),7);
end
nbus=size(bus,1);
nline=size(line,1);
[V0,Q0,s0,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,tgovmech0,f0,dpg0,qplt0,vg0]=unfoldX(x0,SysData);
eq0=Ef0;ed0=zeros(size(Ef0));
[~,~,sNew,dNew,~,eq1New,eq2New,ed1New,ed2New,psidNew,psiqNew,~,EfNew,~,~,~,~,~,~,~]=unfoldX(xNew,SysData);
eqNew=EfNew;edNew=zeros(size(EfNew));
[maxAlpha,segAlpha,dAlpha,nlvl,taylorN,alphaTol,diffTol,diffTolMax,method]=unfoldSimData(SimData);
[pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,MatGV0,MatGV1,MatGRhs0,MatGRhs1]=unfoldSysPara(SysPara);
%
[MatGV0,MatGV1,MatGRhs0,MatGRhs1]=getLinearInterpolatorSyn(syn,synNew,d0,dNew,ed0,ed10,ed20,edNew,ed1New,ed2New,eq0,eq10,eq20,eqNew,eq1New,eq2New,psid0,psiq0,psidNew,psiqNew);

nSys=length(SysDataList);
SysParaLista=SysParaList;
SysDataLista=SysDataList;

for iSys=1:nSys
    SysParaLista{iSys}.MatGV0=MatGV0(grossMaps.synAreaMap==iSys,:);
    SysParaLista{iSys}.MatGV1=MatGV1(grossMaps.synAreaMap==iSys,:);
    SysParaLista{iSys}.MatGRhs0=MatGRhs0(grossMaps.synAreaMap==iSys,:);
    SysParaLista{iSys}.MatGRhs1=MatGRhs1(grossMaps.synAreaMap==iSys,:);
end

if ~isempty(ind)
    [YshInd0,Yshind1]=getLinearInterpolatorInd(nbus,ind,s0,sNew);
    Ysh0=Ysh0+YshInd0;
    Ysh1=Ysh1+Yshind1;
    
    for iSys=1:nSys
        nbusa=size(SysDataList{iSys}.bus,1);
        inda=SysDataList{iSys}.ind;
        s0a=s0(grossMaps.indAreaMap==iSys);
        sNewa=sNew(grossMaps.indAreaMap==iSys);
        [YshInd0a,Yshind1a]=getLinearInterpolatorInd(nbusa,inda,s0a,sNewa);
        SysParaLista{iSys}.Ysh0=SysParaLista{iSys}.Ysh0+YshInd0a;
        SysParaLista{iSys}.Ysh1=SysParaLista{iSys}.Ysh1+Yshind1a;
    end
end

pqx=pq;
pvx=pv;
shuntx=shunt;
indx=ind;
zipx=zip;
synx=syn;

V0x=V0;
Q0x=Q0;
Rzip0x=Rzip0;
Ysh0x=Ysh0;
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

while alphaConfirm<1
    alpha=1-alphaConfirm;
    alphax=alpha;
    if skip==0
        SysDatax=foldSysData(bus,sw,pvx,pqx,shuntx,line,indx,zipx,synx,exc,tg,agc,cac,cluster);
        xx=foldX(SysDatax,V0x,Q0x,s0,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,tgovmech0,f0,dpg0,qplt0,vg0);
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
        
        [V,W,Q]=hemMachinePFSalientcontinueAlgebraicMultiArea(SimData,SysDatax,SysParax,SysDataLista,SysParaLista,links,grossMaps,sysMaps,xx); % note parameters
        heCounter=heCounter+1;
%         [V,W,Q]=hemMachinePFSalientcontinueAlgebraic(...
%             nbus,nline,nlvl,bus,sw,pvx,pqx,shuntx,line,indx,zipx,V0x,Q0x,Pls,Qls,Rzip0x,Rzip1,syn,MatGV0x,MatGV1,MatGRhs0x,MatGRhs1,Ysh0x,Ysh1,Ytr0x,Ytr1,VspSq2);
                
        [Va3,Vb3]=getPadeApproxHelm(V,0);
        [Qa3,Qb3]=getPadeApproxHelm(Q,find(busType==0));
    end
    skip=0;
    alphaLeft=0;
    alphaRight=alpha;
    %     alphax=0.2;
    while 1
%         for i=1:size(V,1)
%             vc(i)=polyval(Va3(i,end:-1:1),alphax)./polyval([Vb3(i,end:-1:1),1],alphax);
%         end
        vc=polyvalVec(Va3(:,end:-1:1),alphax)./polyvalVec([Vb3(:,end:-1:1),ones(size(Vb3,1),1)],alphax);
        vc(isnan(vc))=V(isnan(vc),1)+alphax*V(isnan(vc),2);
%         for i=1:size(Q,1)
%             if busType(i)~=0
%                 qc(i)=polyval(Qa3(i,end:-1:1),alphax)./polyval([Qb3(i,end:-1:1),1],alphax);
%             end
%         end
        qc=polyvalVec(Qa3(:,end:-1:1),alphax)./polyvalVec([Qb3(:,end:-1:1),ones(size(Qb3,1),1)],alphax);
        qc(isnan(qc))=Q(isnan(qc),1)+alphax*Q(isnan(qc),2);
        
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
        
        %         diff=checkEquationBalanceSyn(nbus,nline,bus,sw,pvxx,pqxx,shunt,line,indxx,zipxx,vc,sc,dsc,dyn,syn,Y,Ef,dc,Vsp2,qc);
        
        %         diff=checkEquationBalanceSynAlgebraic(nbus,nline,bus,sw,pvxx,pqxx,shunt,line,zipxx,vc,syn,Y,MatGVxx,MatGRhsxx,Vsp2,qc);
        
        SysDataxx=foldSysData(bus,sw,pvxx,pqxx,shunt,line,indx,zipxx,synxx,exc,tg,agc,cac,cluster);
        xxx=foldX(SysDataxx,vc,qc,s0,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,tgovmech0,f0,dpg0,qplt0,vg0);
        SysParaxx=foldSysPara([],[],[],[],[],[],[],[],Ytr0x+alphax*Ytr1,[],Ysh0x+Ysh1*alphax,[],[Vsp2,zeros(size(Vsp2))],MatGVxx,[],MatGRhsxx,[]);
        
        diff=checkEquationBalanceSynAlgebraic(SysDataxx,SysParaxx,xxx);
        
        if max(abs(diff))<diffTol&&alphax==alpha
            V0x=vc;
            Q0x=qc;
            if heCounter<=1
                Vx=V;Wx=W;Qx=Q;
            end
            
            alphaLeft=alphaRight;
        else
            if max(abs(diff))<diffTol
                V0x=vc;
                Q0x=qc;
                
                if heCounter<=1
                    Vx=V;Wx=W;Qx=Q;
                end
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
    
    alpha=alphaLeft;
    pqx(:,[4,5])=pqx(:,[4,5])+alpha*pqIncr;
    synx(:,[7:10,13:15])=synx(:,[7:10,13:15])+alpha*deltaSyn;
    if ~isempty(pvx);pvx(:,4)=pvx(:,4)+alpha*pvIncr;end
    Rzip0x=Rzip0x+alpha*Rzip1;
    Ysh0x=Ysh0x+alpha*Ysh1;
    Vsp0x=Vsp0x+alpha*VspSq2(:,2);
    Ytr0x=Ytr0x+alpha*Ytr1;
    MatGV0x=MatGV0x+alpha*MatGV1;
    MatGRhs0x=MatGRhs0x+alpha*MatGRhs1;
    
    for iSys=1:nSys
        SysDataLista{iSys}.pq(:,[4,5])=SysDataLista{iSys}.pq(:,[4,5])+alpha*pqIncr(grossMaps.pqAreaMap==iSys,:);
        SysDataLista{iSys}.syn(:,[7:10,13:15])=SysDataLista{iSys}.syn(:,[7:10,13:15])+alpha*deltaSyn(grossMaps.synAreaMap==iSys,:);
        if ~isempty(SysDataLista{iSys}.pv);SysDataLista{iSys}.pv(:,4)=SysDataLista{iSys}.pv(:,4)+alpha*pvIncr(grossMaps.pvAreaMap==iSys,:);end
        SysParaLista{iSys}.Rzip0=SysParaLista{iSys}.Rzip0+alpha*SysParaLista{iSys}.Rzip1;
        SysParaLista{iSys}.Ysh0=SysParaLista{iSys}.Ysh0+alpha*SysParaLista{iSys}.Ysh1;
        
        SysParaLista{iSys}.Ytr0=SysParaLista{iSys}.Ytr0+alpha*SysParaLista{iSys}.Ytr1;
        SysParaLista{iSys}.MatGV0=SysParaLista{iSys}.MatGV0+alpha*SysParaLista{iSys}.MatGV1;
        SysParaLista{iSys}.MatGRhs0=SysParaLista{iSys}.MatGRhs0+alpha*SysParaLista{iSys}.MatGRhs1;
    end
    
    alphaConfirm=alphaConfirm+alpha;
    addLog(['Step=',num2str(alphaConfirm),', added=',num2str(alpha),', (maxdiff<',num2str(diffTol),').'],'INFO');
    
    if alpha==0
        addLog('Step did not move!','INFO');
%         disp('Step did not move!');
        if diffTol>=diffTolMax
            addLog('Max DiffTol reached and not move, exit!','INFO');
%             disp('Max DiffTol reached and not move, exit!');
            break;
        end
        
        diffTol=max([diffTol*1.5,max(abs(diff))+1e-9]);
        if diffTol>diffTolMax
            diffTol=diffTolMax;
        end
        addLog(['Enlarge tol! (Tol=',num2str(diffTol),')'],'INFO');
%         disp(['Enlarge tol! (Tol=',num2str(diffTol),')']);
        skip=1;
    else
        alphas=[alphas;alpha];
    end
end

addLog(['Total Step=',num2str(alphaConfirm),'.'],'INFO');
VSol=V0x;
QSol=Q0x;
% sSol=s0x;
% dSol=d0x;
diffRec=max(abs(diff));
finalAlpha=alphaConfirm;
end