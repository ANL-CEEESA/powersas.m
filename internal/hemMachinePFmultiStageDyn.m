function [t,stateCurve,finalAlpha,alphaList,diffList,exitFlag]=hemMachinePFmultiStageDyn(SimData,SysData,SysPara,x0)
% General interface for invoking dynamic simulation for a segment of time
%
% FUNCTION generalDynSimulation
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
%
% OUTPUT
%   t - A list of time points (starting with 0)
%   stateCurve - A list of states in the order of time
%   finalAlpha - The ending length of this segment of simulation
%   alphaList - Equivalent to t
%   diffList - A list of errors
%   exitFlag - 
%       0  - Success and normally exit
%       -1 - Fail to finish (due to computation errors or failures)
%       1  - System in absolute steady state
%       2  - Generator transients faded away
%       3  - Suggest using error reduction mode
%
global TaskID TaskNum

DynSimFlag=getDynSimFlags();
exitFlag=DynSimFlag.NORMAL;
%
[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
[nState,idxs]...
    =getIndexDyn(SysData);
[nIslands,islands,refs]=searchIslands(bus(:,1),line(:,1:2));
nbus=size(bus,1);
nline=size(line,1);
[V0,Q0,s0,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,Tmech0,f0,dpg0,qplt0,vg0]=unfoldX(x0,SysData);
[maxAlpha,segAlpha,dAlpha,nlvl,taylorN,alphaTol,diffTol,diffTolMax,~,varOpt]=unfoldSimData(SimData);
alphaTolOrig=alphaTol;
[pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,~,~,~,~,Tmech1,Varref1,Ef1,Pm1,Eq11]=unfoldSysPara(SysPara);

busType=zeros(nbus,1);
busType(pv(:,1))=1;
busType(sw(:,1))=2;
% zip(busType(zip(:,1))~=0,10)=0;

agcExt=zeros(nbus,4);
if ~isempty(agc)
    agcExt(agc(:,1),:)=agc;
end
agcRByIsland=accumarray(islands,agcExt(:,4),[nIslands,1]);
isAgcActiveByIsland=agcRByIsland~=0;
isFreqApchZero=isAgcActiveByIsland(islands);
% 
pqx=pq;
pvx=pv;
shuntx=shunt;
indx=ind;
zipx=zip;
% 
V0x=V0;
s0x=s0;
Q0x=Q0;
Rzip0x=Rzip0;
Rind0x=Rind0;
Reind0x=Reind0;
Ysh0x=Ysh0;
Ytr0x=Ytr0;
Vsp0x=real(V0.*conj(V0));
d0x=d0;w0x=w0;eq10x=eq10;eq20x=eq20;ed10x=ed10;ed20x=ed20;psid0x=psid0;psiq0x=psiq0;
Pm0x=Pm0;Ef0x=Ef0;Vavrm0x=Vavrm0;Vavrr0x=Vavrr0;Vavrf0x=Vavrf0;Vavrref0x=Vavrref0;
tgovg0x=tgovg0;tgovm0x=tgovm0;Tmech0x=Tmech0;
f0x=f0;dpg0x=dpg0;qplt0x=qplt0;vg0x=vg0;

VSol=V0;
QSol=Q0;
sSol=s0;
dSol=d0;
wSol=w0;
eq1Sol=eq10;
eq2Sol=eq20;
ed1Sol=ed10;
ed2Sol=ed20;
psiqSol=psiq0;
psidSol=psid0;
PmSol=Pm0;
EfSol=Ef0;
VavrmSol=Vavrm0;
VavrrSol=Vavrr0;
VavrfSol=Vavrf0;
VavrrefSol=Vavrref0;
tgovgSol=tgovg0;
tgovmSol=tgovm0;
TmechSol=Tmech0;
fSol=f0;
dpgSol=dpg0;
qpltSol=qplt0;
vgSol=vg0;

t=0;

alphaConfirm=0;

busType=zeros(nbus,1);
if ~isempty(pv);busType(pv(:,1))=1;end
if ~isempty(sw);busType(sw(:,1))=2;end

nSyn=size(syn,1);
nInd=size(ind,1);
nZip=size(zip,1);

% vc=zeros(nbus,1);
% qc=zeros(size(Q0,1),1);
% sc=zeros(nInd,1);dsc=zeros(nInd,1);
% dc=zeros(nSyn,1);ddc=zeros(nSyn,1);
% wc=zeros(nSyn,1);dwc=zeros(nSyn,1);
% eq1c=zeros(nSyn,1);deq1c=zeros(nSyn,1);
% eq2c=zeros(nSyn,1);deq2c=zeros(nSyn,1);
% ed1c=zeros(nSyn,1);ded1c=zeros(nSyn,1);
% ed2c=zeros(nSyn,1);ded2c=zeros(nSyn,1);
% psidc=zeros(nSyn,1);dpsidc=zeros(nSyn,1);
% psiqc=zeros(nSyn,1);dpsiqc=zeros(nSyn,1);

alphaList=0;
diffList=0;
diffMul=1.05;
alphaPreMult=0.91;
noMoveCnt=0;
absNoMoveCnt=0;
maxNoMoveCnt=3;
maxAbsNoMoveCnt=5;

% generate unique identifier (colamd)
if exist('TaskID','var')&&exist('TaskNum','var')&&~isempty(TaskID)&&~isempty(TaskNum)
taskTag=num2str(TaskID(1));
else
taskTag='';
end

iden=[datestr(clock,30),'_',taskTag,'K',num2str(randi([100000,899999])),'_init'];

jac=getJacobianMod(nbus,nline,bus,sw,pvx,pqx,shuntx,line,zipx,indx,s0,V0,ones(size(ind,1),1));
initJacCond=condest(jac);

if ~isempty(varOpt) && isfield(varOpt,'allowModelSwitch')
    allowModelSwitch=varOpt.allowModelSwitch;
else
    allowModelSwitch=0;
end

if ~isempty(varOpt) && isfield(varOpt,'absT')
    useAbsT=1;
    absT=varOpt.absT;
else
    useAbsT=0;
end

if ~isempty(varOpt) && isfield(varOpt,'allowExit1')
    allowExit1=varOpt.allowExit1;
else
    allowExit1=1;
end

if ~isempty(varOpt) && isfield(varOpt,'useDiffCtrl') && isfield(varOpt,'diffTolCtrl')
    useDiffCtrl=varOpt.useDiffCtrl;
    diffTolCtrl=varOpt.diffTolCtrl;
else
    useDiffCtrl=0;
    diffTolCtrl=0;
end

if ~isempty(varOpt) && isfield(varOpt,'useNomoveCtrl')
    useNomoveCtrl=varOpt.useNomoveCtrl;
else
    useNomoveCtrl=1;
end

while alphaConfirm<maxAlpha-alphaTol/1000
    alpha=min([segAlpha,maxAlpha-alphaConfirm]);
    alphax=alpha/alphaPreMult;
   
    SysDatax=foldSysData(bus,sw,pvx,pqx,shuntx,line,indx,zipx,syn,exc,tg,agc,cac,cluster);
    
    % Check the balance of AE
    
    pqxx=pqx;
    pvxx=pvx;    
    zipxx=zipx;
    if ~isempty(zipxx)
        zipxx(:,5:10)=repmat(Rzip0x,1,6).*zipx(:,5:10);
    end
    
    indxx=indx;
    if ~isempty(indxx)
        indxx(:,15)=(Rind0x).*indx(:,15);
        indxx(:,16)=(Rind0x).*indx(:,16);
        indxx(:,17)=(Rind0x).*indx(:,17);
    end
        
    if useDiffCtrl
        SysDataInit=foldSysData(bus,sw,pvxx,pqxx,shuntx,line,indxx,zipxx,syn,exc,tg,agc,cac,cluster);
        xxx=foldX(SysDataInit,V0x,Q0x,s0x,d0x,w0x,eq10x,eq20x,ed10x,ed20x,psid0x,psiq0x,Pm0x,Ef0x,Vavrm0x,Vavrr0x,Vavrf0x,Vavrref0x,tgovg0x,tgovm0x,Tmech0x,f0x,dpg0x,qplt0x,vg0x);
        [MatGV0,~,MatGRhs0,~]=getLinearInterpolatorSyn(syn,syn,d0x,d0x,Ef0x*0,ed10x,ed20x,Ef0x*0,ed10x,ed20x,Ef0x,eq10x,eq20x,Ef0x,eq10x,eq20x,psid0x,psiq0x,psid0x,psiq0x);
        SysParaInit=foldSysPara([],[],[],[],[],[],[],[],Ytr0x,[],Ysh0x,[],[Vsp0x,zeros(size(Vsp0x))],MatGV0,[],MatGRhs0,[]);
        SysParaInit.nIslands=nIslands;
        SysParaInit.islands=islands;
        SysParaInit.refs=refs; 
        
        diff=checkEquationBalanceSynAlgebraic(SysDataInit,SysParaInit,xxx);
        if max(abs(diff))>max([diffTol/2,diffTolCtrl]) % refine initial point
            SysParaNR=foldSysPara(zeros(size(pqx,1),2),zeros(size(pvx,1),1),ones(nInd,1),zeros(nInd,1),ones(nInd,1),zeros(nInd,1),ones(nZip,1),zeros(nZip,1),Ytr0x,0*Ytr0x,Ysh0x,0*Ysh0x,[Vsp0x,zeros(size(Vsp0x))],MatGV0,[],MatGRhs0,[]);
            SimDataNR=foldSimData(maxAlpha,segAlpha,dAlpha,nlvl,taylorN,alphaTol,diffTol/10,diffTolMax,[]);
            [stateNew,flag,diffRec,loop]=solveAlgebraicNR(SimDataNR,SysDataInit,SysParaNR,xxx,xxx);
            if flag==0
                V0x=stateNew(idxs.vIdx);
                Q0x=stateNew(idxs.qIdx);
                f0x=stateNew(idxs.fIdx);
            end
        end
    end
    
    xx=foldX(SysDatax,V0x,Q0x,s0x,d0x,w0x,eq10x,eq20x,ed10x,ed20x,psid0x,psiq0x,Pm0x,Ef0x,Vavrm0x,Vavrr0x,Vavrf0x,Vavrref0x,tgovg0x,tgovm0x,Tmech0x,f0x,dpg0x,qplt0x,vg0x); 
    pqIncrx=pqIncr;
    if size(pqIncrx,2)>=4
        pqIncrx(:,1:2)=pqIncr(:,1:2)+pqIncr(:,3:4)*2*alphaConfirm;
    end
    SysParax=foldSysPara(pqIncrx,pvIncr,Rind0x,Rind1,Reind0x,Reind1,Rzip0x,Rzip1,Ytr0x,Ytr1,Ysh0x,Ysh1,[Vsp0x,VspSq2(:,2)],[],[],[],[],Tmech1,Varref1,Ef1,Pm1,Eq11);
    SysParax.iden=iden;
    if exist([iden,'.mat'],'file')
        if ~exist('p_amd','var')
            load([iden,'.mat']);
        end
        SysParax.p_amd=p_amd;
    else
        SysParax.p_amd=[];
    end
    SysParax.nIslands=nIslands;
    SysParax.islands=islands;
    SysParax.refs=refs;  
    
    [V,Q,s,d,w,eq1,eq2,ed1,ed2,psid,psiq,Pm,Ef,Vavrm,Vavrr,Vavrf,Vavrref,tgovg,tgovm,Tmech,f,dpg,qplt,vg]=...
        hemMachinePFSalientcontinueDyn(SimData,SysDatax,SysParax,xx);    
    
    ds=s(:,2:end).*repmat(1:(size(s,2)-1),size(s,1),1);
    dd=d(:,2:end).*repmat(1:(size(d,2)-1),size(d,1),1);
    dw=w(:,2:end).*repmat(1:(size(w,2)-1),size(w,1),1);    
    deq1=eq1(:,2:end).*repmat(1:(size(eq1,2)-1),size(eq1,1),1);
    deq2=eq2(:,2:end).*repmat(1:(size(eq2,2)-1),size(eq2,1),1);
    ded1=ed1(:,2:end).*repmat(1:(size(ed1,2)-1),size(ed1,1),1);
    ded2=ed2(:,2:end).*repmat(1:(size(ed2,2)-1),size(ed2,1),1);
    dpsid=psid(:,2:end).*repmat(1:(size(psid,2)-1),size(psid,1),1);
    dpsiq=psiq(:,2:end).*repmat(1:(size(psiq,2)-1),size(psiq,1),1);
    dPm=Pm(:,2:end).*repmat(1:(size(Pm,2)-1),size(Pm,1),1);
    
    dVavrm=Vavrm(:,2:end).*repmat(1:(size(Vavrm,2)-1),size(Vavrm,1),1);
    dVavrr=Vavrr(:,2:end).*repmat(1:(size(Vavrr,2)-1),size(Vavrr,1),1);
    dVavrf=Vavrf(:,2:end).*repmat(1:(size(Vavrf,2)-1),size(Vavrf,1),1);
    dtgovg=tgovg(:,2:end).*repmat(1:(size(tgovg,2)-1),size(tgovg,1),1);
    dTmech=Tmech(:,2:end).*repmat(1:(size(Tmech,2)-1),size(Tmech,1),1);
    
    df=f(:,2:end).*repmat(1:(size(f,2)-1),size(f,1),1); %agc
    ddpg=dpg(:,2:end).*repmat(1:(size(dpg,2)-1),size(dpg,1),1); %agc
    % TODO: avc
    
    [Va,Vb]=getPadeApproxHelm(V,0);
    idxNanV=find(isnan(Va(:,1))|isnan(Vb(:,1)));Va(idxNanV,:)=0;Vb(idxNanV,:)=0;Va(idxNanV,1:3)=V(idxNanV,1:3);
    [Qa,Qb]=getPadeApproxHelm(Q,find(busType==0));
    idxNanQ=find(isnan(Qa(:,1))|isnan(Qb(:,1)));Qa(idxNanQ,:)=0;Qb(idxNanQ,:)=0;Qa(idxNanQ,1:3)=Q(idxNanQ,1:3);
    
    [sa,sb]=getPadeApproxHelm(s,0);
    [dsa,dsb]=getPadeApproxHelm(ds,0);
    [da,db]=getPadeApproxHelm(d,0);
    [dda,ddb]=getPadeApproxHelm(dd,0);
    [wa,wb]=getPadeApproxHelm(w,0);
    [dwa,dwb]=getPadeApproxHelm(dw,0);
    [eq1a,eq1b]=getPadeApproxHelm(eq1,0);
    [deq1a,deq1b]=getPadeApproxHelm(deq1,0);
    [eq2a,eq2b]=getPadeApproxHelm(eq2,0);
    [deq2a,deq2b]=getPadeApproxHelm(deq2,0);
    [ed1a,ed1b]=getPadeApproxHelm(ed1,0);
    [ded1a,ded1b]=getPadeApproxHelm(ded1,0);
    [ed2a,ed2b]=getPadeApproxHelm(ed2,0);
    [ded2a,ded2b]=getPadeApproxHelm(ded2,0);
    [psida,psidb]=getPadeApproxHelm(psid,0);
    [dpsida,dpsidb]=getPadeApproxHelm(dpsid,0);
    [psiqa,psiqb]=getPadeApproxHelm(psiq,0);
    [dpsiqa,dpsiqb]=getPadeApproxHelm(dpsiq,0);
    
    [Efa,Efb]=getPadeApproxHelm(Ef,0);
    [Pma,Pmb]=getPadeApproxHelm(Pm,0);
    [dPma,dPmb]=getPadeApproxHelm(dPm,0);
    [Vavrma,Vavrmb]=getPadeApproxHelm(Vavrm,0);
    [dVavrma,dVavrmb]=getPadeApproxHelm(dVavrm,0);
    [Vavrra,Vavrrb]=getPadeApproxHelm(Vavrr,0);
    [dVavrra,dVavrrb]=getPadeApproxHelm(dVavrr,0);
    [Vavrfa,Vavrfb]=getPadeApproxHelm(Vavrf,0);
    [dVavrfa,dVavrfb]=getPadeApproxHelm(dVavrf,0);
    [Vavrrefa,Vavrrefb]=getPadeApproxHelm(Vavrref,0);
    [tgovga,tgovgb]=getPadeApproxHelm(tgovg,0);
    [dtgovga,dtgovgb]=getPadeApproxHelm(dtgovg,0);
    [tgovma,tgovmb]=getPadeApproxHelm(tgovm,0);
    [Tmecha,Tmechb]=getPadeApproxHelm(Tmech,0);  
    [dTmecha,dTmechb]=getPadeApproxHelm(dTmech,0);    
    
    [fa,fb]=getPadeApproxHelm(f,0); % agc
    [dfa,dfb]=getPadeApproxHelm(df,0); % agc
    [dpga,dpgb]=getPadeApproxHelm(dpg,0); % agc
    [ddpga,ddpgb]=getPadeApproxHelm(ddpg,0); % agc
    
    [qplta,qpltb]=getPadeApproxHelm(qplt,0);
    [vga,vgb]=getPadeApproxHelm(vg,0);
    
    if ~isempty(ind)
        idxNanS=find(isnan(sa(:,1))|isnan(sb(:,1)));sa(idxNanS,:)=0;sb(idxNanS,:)=0;sa(idxNanS,1:2)=s(idxNanS,1:2);
        idxNandS=find(isnan(dsa(:,1))|isnan(dsb(:,1)));dsa(idxNandS,:)=0;dsb(idxNandS,:)=0;dsa(idxNandS,1:2)=ds(idxNandS,1:2);
    end    
    if ~isempty(syn)
        idxNand=find(isnan(da(:,1))|isnan(db(:,1)));da(idxNand,:)=0;db(idxNand,:)=0;da(idxNand,1:2)=d(idxNand,1:2);
        idxNandd=find(isnan(dda(:,1))|isnan(ddb(:,1)));dda(idxNandd,:)=0;ddb(idxNandd,:)=0;dda(idxNandd,1:2)=dd(idxNandd,1:2);
        idxNanw=find(isnan(wa(:,1))|isnan(wb(:,1)));wa(idxNanw,:)=0;wb(idxNanw,:)=0;wa(idxNanw,1:2)=w(idxNanw,1:2);
        idxNandw=find(isnan(dwa(:,1))|isnan(dwb(:,1)));dwa(idxNandw,:)=0;dwb(idxNandw,:)=0;dwa(idxNandw,1:2)=dw(idxNandw,1:2);
        idxNaneq1=find(isnan(eq1a(:,1))|isnan(eq1b(:,1)));eq1a(idxNaneq1,:)=0;eq1b(idxNaneq1,:)=0;eq1a(idxNaneq1,1:2)=eq1(idxNaneq1,1:2);
        idxNandeq1=find(isnan(deq1a(:,1))|isnan(deq1b(:,1)));deq1a(idxNandeq1,:)=0;deq1b(idxNandeq1,:)=0;deq1a(idxNandeq1,1:2)=deq1(idxNandeq1,1:2);
        idxNaneq2=find(isnan(eq2a(:,1))|isnan(eq2b(:,1)));eq2a(idxNaneq2,:)=0;eq2b(idxNaneq2,:)=0;eq2a(idxNaneq2,1:2)=eq2(idxNaneq2,1:2);
        idxNandeq2=find(isnan(deq2a(:,1))|isnan(deq2b(:,1)));deq2a(idxNandeq2,:)=0;deq2b(idxNandeq2,:)=0;deq2a(idxNandeq2,1:2)=deq2(idxNandeq2,1:2);
        idxNaned1=find(isnan(ed1a(:,1))|isnan(ed1b(:,1)));ed1a(idxNaned1,:)=0;ed1b(idxNaned1,:)=0;ed1a(idxNaned1,1:2)=ed1(idxNaned1,1:2);
        idxNanded1=find(isnan(ded1a(:,1))|isnan(ded1b(:,1)));ded1a(idxNanded1,:)=0;ded1b(idxNanded1,:)=0;ded1a(idxNanded1,1:2)=ded1(idxNanded1,1:2);
        idxNaned2=find(isnan(ed2a(:,1))|isnan(ed2b(:,1)));ed2a(idxNaned2,:)=0;ed2b(idxNaned2,:)=0;ed2a(idxNaned2,1:2)=ed2(idxNaned2,1:2);
        idxNanded2=find(isnan(ded2a(:,1))|isnan(ded2b(:,1)));ded2a(idxNanded2,:)=0;ded2b(idxNanded2,:)=0;ded2a(idxNanded2,1:2)=ded2(idxNanded2,1:2);
        idxNanpsid=find(isnan(psida(:,1))|isnan(psidb(:,1)));psida(idxNanpsid,:)=0;psidb(idxNanpsid,:)=0;psida(idxNanpsid,1:2)=psid(idxNanpsid,1:2);
        idxNandpsid=find(isnan(dpsida(:,1))|isnan(dpsidb(:,1)));dpsida(idxNandpsid,:)=0;dpsidb(idxNandpsid,:)=0;dpsida(idxNandpsid,1:2)=dpsid(idxNandpsid,1:2);
        idxNanpsiq=find(isnan(psiqa(:,1))|isnan(psiqb(:,1)));psiqa(idxNanpsiq,:)=0;psiqb(idxNanpsiq,:)=0;psiqa(idxNanpsiq,1:2)=psiq(idxNanpsiq,1:2);
        idxNandpsiq=find(isnan(dpsiqa(:,1))|isnan(dpsiqb(:,1)));dpsiqa(idxNandpsiq,:)=0;dpsiqb(idxNandpsiq,:)=0;dpsiqa(idxNandpsiq,1:2)=dpsiq(idxNandpsiq,1:2);
        idxNanEf=find(isnan(Efa(:,1))|isnan(Efb(:,1)));Efa(idxNanEf,:)=0;Efb(idxNanEf,:)=0;Efa(idxNanEf,1:2)=Ef(idxNanEf,1:2);
        idxNanPm=find(isnan(Pma(:,1))|isnan(Pmb(:,1)));Pma(idxNanPm,:)=0;Pmb(idxNanPm,:)=0;Pma(idxNanPm,1:2)=Pm(idxNanPm,1:2);
        idxNandPm=find(isnan(dPma(:,1))|isnan(dPmb(:,1)));dPma(idxNandPm,:)=0;dPmb(idxNandPm,:)=0;dPma(idxNandPm,1:2)=dPm(idxNandPm,1:2);
    end
    if ~isempty(exc)        
        idxNanVavrm=find(isnan(Vavrma(:,1))|isnan(Vavrmb(:,1)));Vavrma(idxNanVavrm,:)=0;Vavrmb(idxNanVavrm,:)=0;Vavrma(idxNanVavrm,1:2)=Vavrm(idxNanVavrm,1:2);
        idxNandVavrm=find(isnan(dVavrma(:,1))|isnan(dVavrmb(:,1)));dVavrma(idxNandVavrm,:)=0;dVavrmb(idxNandVavrm,:)=0;dVavrma(idxNandVavrm,1:2)=dVavrm(idxNandVavrm,1:2);
        idxNanVavrr=find(isnan(Vavrra(:,1))|isnan(Vavrrb(:,1)));Vavrra(idxNanVavrr,:)=0;Vavrrb(idxNanVavrr,:)=0;Vavrra(idxNanVavrr,1:2)=Vavrr(idxNanVavrr,1:2);
        idxNandVavrr=find(isnan(dVavrra(:,1))|isnan(dVavrrb(:,1)));dVavrra(idxNandVavrr,:)=0;dVavrrb(idxNandVavrr,:)=0;dVavrra(idxNandVavrr,1:2)=dVavrr(idxNandVavrr,1:2);
        idxNanVavrf=find(isnan(Vavrfa(:,1))|isnan(Vavrfb(:,1)));Vavrfa(idxNanVavrf,:)=0;Vavrfb(idxNanVavrf,:)=0;Vavrfa(idxNanVavrf,1:2)=Vavrf(idxNanVavrf,1:2);
        idxNandVavrf=find(isnan(dVavrfa(:,1))|isnan(dVavrfb(:,1)));dVavrfa(idxNandVavrf,:)=0;dVavrfb(idxNandVavrf,:)=0;dVavrfa(idxNandVavrf,1:2)=dVavrf(idxNandVavrf,1:2);
        idxNanVavrref=find(isnan(Vavrrefa(:,1))|isnan(Vavrrefb(:,1)));Vavrrefa(idxNanVavrref,:)=0;Vavrrefb(idxNanVavrref,:)=0;Vavrrefa(idxNanVavrref,1:2)=Vavrref(idxNanVavrref,1:2);
    end
    if ~isempty(tg)
        idxNantgovg=find(isnan(tgovga(:,1))|isnan(tgovgb(:,1)));tgovga(idxNantgovg,:)=0;tgovgb(idxNantgovg,:)=0;tgovga(idxNantgovg,1:2)=tgovg(idxNantgovg,1:2);
        idxNandtgovg=find(isnan(dtgovga(:,1))|isnan(dtgovgb(:,1)));dtgovga(idxNandtgovg,:)=0;dtgovgb(idxNandtgovg,:)=0;dtgovga(idxNandtgovg,1:2)=dtgovg(idxNandtgovg,1:2); 
        idxNantgovm=find(isnan(tgovma(:,1))|isnan(tgovmb(:,1)));tgovma(idxNantgovm,:)=0;tgovmb(idxNantgovm,:)=0;tgovma(idxNantgovm,1:2)=tgovm(idxNantgovm,1:2); 
        idxNanTmech=find(isnan(Tmecha(:,1))|isnan(Tmechb(:,1)));Tmecha(idxNanTmech,:)=0;Tmechb(idxNanTmech,:)=0;Tmecha(idxNanTmech,1:2)=Tmech(idxNanTmech,1:2);
        idxNandTmech=find(isnan(dTmecha(:,1))|isnan(dTmechb(:,1)));dTmecha(idxNandTmech,:)=0;dTmechb(idxNandTmech,:)=0;dTmecha(idxNandTmech,1:2)=dTmech(idxNandTmech,1:2);
    end
    
    idxNanf=find(isnan(fa(:,1))|isnan(fb(:,1)));fa(idxNanf,:)=0;fb(idxNanf,:)=0;fa(idxNanf,1:2)=f(idxNanf,1:2);
    idxNandf=find(isnan(dfa(:,1))|isnan(dfb(:,1)));dfa(idxNandf,:)=0;dfb(idxNandf,:)=0;dfa(idxNandf,1:2)=df(idxNandf,1:2);
    idxNandpg=find(isnan(dpga(:,1))|isnan(dpgb(:,1)));dpga(idxNandpg,:)=0;dpgb(idxNandpg,:)=0;dpga(idxNandpg,1:2)=dpg(idxNandpg,1:2);
    idxNanddpg=find(isnan(ddpga(:,1))|isnan(ddpgb(:,1)));ddpga(idxNanddpg,:)=0;ddpgb(idxNanddpg,:)=0;ddpga(idxNanddpg,1:2)=ddpg(idxNanddpg,1:2);
        
    if ~isempty(cac)&&~isempty(cluster)
        
    end
    
    alphaLeft=0;
    alphaRight=alpha/alphaPreMult;
    %     alphax=0.2;
    while 1
        
        vc=polyvalVec(Va(:,end:-1:1),alphax)./polyvalVec([Vb(:,end:-1:1),ones(size(Vb,1),1)],alphax);
        
        qc=polyvalVec(Qa(:,end:-1:1),alphax)./polyvalVec([Qb(:,end:-1:1),ones(size(Qb,1),1)],alphax);
        qc(busType==0,:)=0;
        
        sc=polyvalVec(sa(:,end:-1:1),alphax)./polyvalVec([sb(:,end:-1:1),ones(size(sb,1),1)],alphax);
        dsc=polyvalVec(dsa(:,end:-1:1),alphax)./polyvalVec([dsb(:,end:-1:1),ones(size(dsb,1),1)],alphax);
        
        dc=polyvalVec(da(:,end:-1:1),alphax)./polyvalVec([db(:,end:-1:1),ones(size(db,1),1)],alphax);
        ddc=polyvalVec(dda(:,end:-1:1),alphax)./polyvalVec([ddb(:,end:-1:1),ones(size(ddb,1),1)],alphax);
        wc=polyvalVec(wa(:,end:-1:1),alphax)./polyvalVec([wb(:,end:-1:1),ones(size(wb,1),1)],alphax);
        dwc=polyvalVec(dwa(:,end:-1:1),alphax)./polyvalVec([dwb(:,end:-1:1),ones(size(dwb,1),1)],alphax);
        eq1c=polyvalVec(eq1a(:,end:-1:1),alphax)./polyvalVec([eq1b(:,end:-1:1),ones(size(eq1b,1),1)],alphax);
        deq1c=polyvalVec(deq1a(:,end:-1:1),alphax)./polyvalVec([deq1b(:,end:-1:1),ones(size(deq1b,1),1)],alphax);
        eq2c=polyvalVec(eq2a(:,end:-1:1),alphax)./polyvalVec([eq2b(:,end:-1:1),ones(size(eq2b,1),1)],alphax);
        deq2c=polyvalVec(deq2a(:,end:-1:1),alphax)./polyvalVec([deq2b(:,end:-1:1),ones(size(deq2b,1),1)],alphax);
        ed1c=polyvalVec(ed1a(:,end:-1:1),alphax)./polyvalVec([ed1b(:,end:-1:1),ones(size(ed1b,1),1)],alphax);
        ded1c=polyvalVec(ded1a(:,end:-1:1),alphax)./polyvalVec([ded1b(:,end:-1:1),ones(size(ded1b,1),1)],alphax);
        ed2c=polyvalVec(ed2a(:,end:-1:1),alphax)./polyvalVec([ed2b(:,end:-1:1),ones(size(ed2b,1),1)],alphax);
        ded2c=polyvalVec(ded2a(:,end:-1:1),alphax)./polyvalVec([ded2b(:,end:-1:1),ones(size(ded2b,1),1)],alphax);
        psidc=polyvalVec(psida(:,end:-1:1),alphax)./polyvalVec([psidb(:,end:-1:1),ones(size(psidb,1),1)],alphax);
        dpsidc=polyvalVec(dpsida(:,end:-1:1),alphax)./polyvalVec([dpsidb(:,end:-1:1),ones(size(dpsidb,1),1)],alphax);
        psiqc=polyvalVec(psiqa(:,end:-1:1),alphax)./polyvalVec([psiqb(:,end:-1:1),ones(size(psiqb,1),1)],alphax);
        dpsiqc=polyvalVec(dpsiqa(:,end:-1:1),alphax)./polyvalVec([dpsiqb(:,end:-1:1),ones(size(dpsiqb,1),1)],alphax);
                
        Pmc=polyvalVec(Pma(:,end:-1:1),alphax)./polyvalVec([Pmb(:,end:-1:1),ones(size(Pmb,1),1)],alphax);
        dPmc=polyvalVec(dPma(:,end:-1:1),alphax)./polyvalVec([dPmb(:,end:-1:1),ones(size(dPmb,1),1)],alphax);
        Efc=polyvalVec(Efa(:,end:-1:1),alphax)./polyvalVec([Efb(:,end:-1:1),ones(size(Efb,1),1)],alphax);
        Vavrmc=polyvalVec(Vavrma(:,end:-1:1),alphax)./polyvalVec([Vavrmb(:,end:-1:1),ones(size(Vavrmb,1),1)],alphax);
        dVavrmc=polyvalVec(dVavrma(:,end:-1:1),alphax)./polyvalVec([dVavrmb(:,end:-1:1),ones(size(dVavrmb,1),1)],alphax);
        Vavrrc=polyvalVec(Vavrra(:,end:-1:1),alphax)./polyvalVec([Vavrrb(:,end:-1:1),ones(size(Vavrrb,1),1)],alphax);
        dVavrrc=polyvalVec(dVavrra(:,end:-1:1),alphax)./polyvalVec([dVavrrb(:,end:-1:1),ones(size(dVavrrb,1),1)],alphax);
        Vavrfc=polyvalVec(Vavrfa(:,end:-1:1),alphax)./polyvalVec([Vavrfb(:,end:-1:1),ones(size(Vavrfb,1),1)],alphax);
        dVavrfc=polyvalVec(dVavrfa(:,end:-1:1),alphax)./polyvalVec([dVavrfb(:,end:-1:1),ones(size(dVavrfb,1),1)],alphax);
        Vavrrefc=polyvalVec(Vavrrefa(:,end:-1:1),alphax)./polyvalVec([Vavrrefb(:,end:-1:1),ones(size(Vavrrefb,1),1)],alphax);
        tgovgc=polyvalVec(tgovga(:,end:-1:1),alphax)./polyvalVec([tgovgb(:,end:-1:1),ones(size(tgovgb,1),1)],alphax);
        dtgovgc=polyvalVec(dtgovga(:,end:-1:1),alphax)./polyvalVec([dtgovgb(:,end:-1:1),ones(size(dtgovgb,1),1)],alphax);
        tgovmc=polyvalVec(tgovma(:,end:-1:1),alphax)./polyvalVec([tgovmb(:,end:-1:1),ones(size(tgovmb,1),1)],alphax);
        Tmechc=polyvalVec(Tmecha(:,end:-1:1),alphax)./polyvalVec([Tmechb(:,end:-1:1),ones(size(Tmechb,1),1)],alphax);
        dTmechc=polyvalVec(dTmecha(:,end:-1:1),alphax)./polyvalVec([dTmechb(:,end:-1:1),ones(size(dTmechb,1),1)],alphax);
        
        fc=polyvalVec(fa(:,end:-1:1),alphax)./polyvalVec([fb(:,end:-1:1),ones(size(fb,1),1)],alphax);
        dfc=polyvalVec(dfa(:,end:-1:1),alphax)./polyvalVec([dfb(:,end:-1:1),ones(size(dfb,1),1)],alphax);
        dpgc=polyvalVec(dpga(:,end:-1:1),alphax)./polyvalVec([dpgb(:,end:-1:1),ones(size(dpgb,1),1)],alphax);
        ddpgc=polyvalVec(ddpga(:,end:-1:1),alphax)./polyvalVec([ddpgb(:,end:-1:1),ones(size(ddpgb,1),1)],alphax);
                
        qpltc=polyvalVec(qplta(:,end:-1:1),alphax)./polyvalVec([qpltb(:,end:-1:1),ones(size(qplt,1),1)],alphax);
        vgc=polyvalVec(vga(:,end:-1:1),alphax)./polyvalVec([vgb(:,end:-1:1),ones(size(vg,1),1)],alphax);
%         Y=Ytr0x+alphax*Ytr1+sparse(1:nbus,1:nbus,Ysh0x+Ysh1*alphax,nbus,nbus);
        Vsp2=Vsp0x+alphax*VspSq2(:,2);
        
        pqxx=pqx;
        pqxx(:,[4,5])=pqx(:,[4,5])+alphax*pqIncrx(:,1:2);
        if size(pqIncrx,2)>=4
            pqxx(:,[4,5])=pqxx(:,[4,5])+alphax*alphax*pqIncrx(:,3:4);
        end
        pvxx=pvx;
        if ~isempty(pvx);pvxx(:,4)=pvx(:,4)+alphax*pvIncr;end
        
        zipxx=zipx;
        if ~isempty(zipxx)
            zipxx(:,5:10)=repmat(Rzip0x+alphax*Rzip1,1,6).*zipx(:,5:10);
        end
        
        indxx=indx;
        if ~isempty(indxx)
            indxx(:,15)=(Rind0x+alphax*Rind1).*indx(:,15);
            indxx(:,16)=(Rind0x+alphax*Rind1).*indx(:,16);
            indxx(:,17)=(Rind0x+alphax*Rind1).*indx(:,17);         
        end
        
        SysDataxx=foldSysData(bus,sw,pvxx,pqxx,shunt,line,indxx,zipxx,syn,exc,tg,agc,cac,cluster);
        xc=foldX(SysDataxx,vc,qc,sc,dc,wc,eq1c,eq2c,ed1c,ed2c,psidc,psiqc,Pmc,Efc,Vavrmc,Vavrrc,Vavrfc,Vavrrefc,tgovgc,tgovmc,Tmechc,fc,dpgc,qpltc,vgc);
        dxc=foldX(SysDataxx,0*vc,0*qc,dsc,ddc,dwc,deq1c,deq2c,ded1c,ded2c,dpsidc,dpsiqc,dPmc,0*Efc,dVavrmc,dVavrrc,dVavrfc,0*Vavrrefc,dtgovgc,0*tgovmc,dTmechc,dfc,ddpgc,0*qpltc,0*vgc);
        SysParaxx=foldSysPara([],[],[],[],[],[],[],[],Ytr0x+alphax*Ytr1,[],Ysh0x+Ysh1*alphax,[],[Vsp2,zeros(size(Vsp2))],[],[],[],[],Tmech1,Varref1,Ef1,Pm1,Eq11);
        SysParaxx.nIslands=nIslands;
        SysParaxx.islands=islands;
        SysParaxx.refs=refs;
        
        diff=checkEquationBalanceSynDyn(SysDataxx,SysParaxx,xc,dxc);
        absDiff=max(abs(diff));
        
        if absDiff<diffTol&&alphax==alpha/alphaPreMult
            alphaLeft=alphaRight;
        else
            if absDiff<diffTol
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
    
    diffParadigm=absDiff+1e-12;
    countCheckAlpha=1;
    alphax=alphaLeft*alphaPreMult;
    if alpha-alphax<alphaTol/1000
        alphax=alpha;
    end
    maxCount=10;
    while countCheckAlpha<maxCount        
        vc=polyvalVec(Va(:,end:-1:1),alphax)./polyvalVec([Vb(:,end:-1:1),ones(size(Vb,1),1)],alphax);
        
        qc=polyvalVec(Qa(:,end:-1:1),alphax)./polyvalVec([Qb(:,end:-1:1),ones(size(Qb,1),1)],alphax);
        qc(busType==0,:)=0;
        
        sc=polyvalVec(sa(:,end:-1:1),alphax)./polyvalVec([sb(:,end:-1:1),ones(size(sb,1),1)],alphax);
        dsc=polyvalVec(dsa(:,end:-1:1),alphax)./polyvalVec([dsb(:,end:-1:1),ones(size(dsb,1),1)],alphax);
        
        dc=polyvalVec(da(:,end:-1:1),alphax)./polyvalVec([db(:,end:-1:1),ones(size(db,1),1)],alphax);
        ddc=polyvalVec(dda(:,end:-1:1),alphax)./polyvalVec([ddb(:,end:-1:1),ones(size(ddb,1),1)],alphax);
        wc=polyvalVec(wa(:,end:-1:1),alphax)./polyvalVec([wb(:,end:-1:1),ones(size(wb,1),1)],alphax);
        dwc=polyvalVec(dwa(:,end:-1:1),alphax)./polyvalVec([dwb(:,end:-1:1),ones(size(dwb,1),1)],alphax);
        eq1c=polyvalVec(eq1a(:,end:-1:1),alphax)./polyvalVec([eq1b(:,end:-1:1),ones(size(eq1b,1),1)],alphax);
        deq1c=polyvalVec(deq1a(:,end:-1:1),alphax)./polyvalVec([deq1b(:,end:-1:1),ones(size(deq1b,1),1)],alphax);
        eq2c=polyvalVec(eq2a(:,end:-1:1),alphax)./polyvalVec([eq2b(:,end:-1:1),ones(size(eq2b,1),1)],alphax);
        deq2c=polyvalVec(deq2a(:,end:-1:1),alphax)./polyvalVec([deq2b(:,end:-1:1),ones(size(deq2b,1),1)],alphax);
        ed1c=polyvalVec(ed1a(:,end:-1:1),alphax)./polyvalVec([ed1b(:,end:-1:1),ones(size(ed1b,1),1)],alphax);
        ded1c=polyvalVec(ded1a(:,end:-1:1),alphax)./polyvalVec([ded1b(:,end:-1:1),ones(size(ded1b,1),1)],alphax);
        ed2c=polyvalVec(ed2a(:,end:-1:1),alphax)./polyvalVec([ed2b(:,end:-1:1),ones(size(ed2b,1),1)],alphax);
        ded2c=polyvalVec(ded2a(:,end:-1:1),alphax)./polyvalVec([ded2b(:,end:-1:1),ones(size(ded2b,1),1)],alphax);
        psidc=polyvalVec(psida(:,end:-1:1),alphax)./polyvalVec([psidb(:,end:-1:1),ones(size(psidb,1),1)],alphax);
        dpsidc=polyvalVec(dpsida(:,end:-1:1),alphax)./polyvalVec([dpsidb(:,end:-1:1),ones(size(dpsidb,1),1)],alphax);
        psiqc=polyvalVec(psiqa(:,end:-1:1),alphax)./polyvalVec([psiqb(:,end:-1:1),ones(size(psiqb,1),1)],alphax);
        dpsiqc=polyvalVec(dpsiqa(:,end:-1:1),alphax)./polyvalVec([dpsiqb(:,end:-1:1),ones(size(dpsiqb,1),1)],alphax);
                
        Pmc=polyvalVec(Pma(:,end:-1:1),alphax)./polyvalVec([Pmb(:,end:-1:1),ones(size(Pmb,1),1)],alphax);
        dPmc=polyvalVec(dPma(:,end:-1:1),alphax)./polyvalVec([dPmb(:,end:-1:1),ones(size(dPmb,1),1)],alphax);
        Efc=polyvalVec(Efa(:,end:-1:1),alphax)./polyvalVec([Efb(:,end:-1:1),ones(size(Efb,1),1)],alphax);
        Vavrmc=polyvalVec(Vavrma(:,end:-1:1),alphax)./polyvalVec([Vavrmb(:,end:-1:1),ones(size(Vavrmb,1),1)],alphax);
        dVavrmc=polyvalVec(dVavrma(:,end:-1:1),alphax)./polyvalVec([dVavrmb(:,end:-1:1),ones(size(dVavrmb,1),1)],alphax);
        Vavrrc=polyvalVec(Vavrra(:,end:-1:1),alphax)./polyvalVec([Vavrrb(:,end:-1:1),ones(size(Vavrrb,1),1)],alphax);
        dVavrrc=polyvalVec(dVavrra(:,end:-1:1),alphax)./polyvalVec([dVavrrb(:,end:-1:1),ones(size(dVavrrb,1),1)],alphax);
        Vavrfc=polyvalVec(Vavrfa(:,end:-1:1),alphax)./polyvalVec([Vavrfb(:,end:-1:1),ones(size(Vavrfb,1),1)],alphax);
        dVavrfc=polyvalVec(dVavrfa(:,end:-1:1),alphax)./polyvalVec([dVavrfb(:,end:-1:1),ones(size(dVavrfb,1),1)],alphax);
        Vavrrefc=polyvalVec(Vavrrefa(:,end:-1:1),alphax)./polyvalVec([Vavrrefb(:,end:-1:1),ones(size(Vavrrefb,1),1)],alphax);
        tgovgc=polyvalVec(tgovga(:,end:-1:1),alphax)./polyvalVec([tgovgb(:,end:-1:1),ones(size(tgovgb,1),1)],alphax);
        dtgovgc=polyvalVec(dtgovga(:,end:-1:1),alphax)./polyvalVec([dtgovgb(:,end:-1:1),ones(size(dtgovgb,1),1)],alphax);
        tgovmc=polyvalVec(tgovma(:,end:-1:1),alphax)./polyvalVec([tgovmb(:,end:-1:1),ones(size(tgovmb,1),1)],alphax);
        Tmechc=polyvalVec(Tmecha(:,end:-1:1),alphax)./polyvalVec([Tmechb(:,end:-1:1),ones(size(Tmechb,1),1)],alphax);
        dTmechc=polyvalVec(dTmecha(:,end:-1:1),alphax)./polyvalVec([dTmechb(:,end:-1:1),ones(size(dTmechb,1),1)],alphax);
        
        fc=polyvalVec(fa(:,end:-1:1),alphax)./polyvalVec([fb(:,end:-1:1),ones(size(fb,1),1)],alphax);
        dfc=polyvalVec(dfa(:,end:-1:1),alphax)./polyvalVec([dfb(:,end:-1:1),ones(size(dfb,1),1)],alphax);
        dpgc=polyvalVec(dpga(:,end:-1:1),alphax)./polyvalVec([dpgb(:,end:-1:1),ones(size(dpgb,1),1)],alphax);
        ddpgc=polyvalVec(ddpga(:,end:-1:1),alphax)./polyvalVec([ddpgb(:,end:-1:1),ones(size(ddpgb,1),1)],alphax);
                
        qpltc=polyvalVec(qplta(:,end:-1:1),alphax)./polyvalVec([qpltb(:,end:-1:1),ones(size(qplt,1),1)],alphax);
        vgc=polyvalVec(vga(:,end:-1:1),alphax)./polyvalVec([vgb(:,end:-1:1),ones(size(vg,1),1)],alphax);
        
%         Y=Ytr0x+alphax*Ytr1+sparse(1:nbus,1:nbus,Ysh0x+Ysh1*alphax,nbus,nbus);
        Vsp2=Vsp0x+alphax*VspSq2(:,2);
        
        pqxx=pqx;
        pqxx(:,[4,5])=pqx(:,[4,5])+alphax*pqIncrx(:,1:2);
        if size(pqIncrx,2)>=4
            pqxx(:,[4,5])=pqxx(:,[4,5])+alphax*alphax*pqIncrx(:,3:4);
        end
        pvxx=pvx;
        if ~isempty(pvx);pvxx(:,4)=pvx(:,4)+alphax*pvIncr;end
        
        zipxx=zipx;
        if ~isempty(zipxx)
            zipxx(:,5:10)=repmat(Rzip0x+alphax*Rzip1,1,6).*zipx(:,5:10);
        end
        
        indxx=indx;
        if ~isempty(indxx)
            indxx(:,15)=(Rind0x+alphax*Rind1).*indx(:,15);
            indxx(:,16)=(Rind0x+alphax*Rind1).*indx(:,16);
            indxx(:,17)=(Rind0x+alphax*Rind1).*indx(:,17);         
        end
        
        SysDataxx=foldSysData(bus,sw,pvxx,pqxx,shunt,line,indxx,zipxx,syn,exc,tg,agc,cac,cluster);
        xc=foldX(SysDataxx,vc,qc,sc,dc,wc,eq1c,eq2c,ed1c,ed2c,psidc,psiqc,Pmc,Efc,Vavrmc,Vavrrc,Vavrfc,Vavrrefc,tgovgc,tgovmc,Tmechc,fc,dpgc,qpltc,vgc);
        dxc=foldX(SysDataxx,0*vc,0*qc,dsc,ddc,dwc,deq1c,deq2c,ded1c,ded2c,dpsidc,dpsiqc,dPmc,0*Efc,dVavrmc,dVavrrc,dVavrfc,0*Vavrrefc,dtgovgc,0*tgovmc,dTmechc,dfc,ddpgc,0*qpltc,0*vgc);
        SysParaxx=foldSysPara([],[],[],[],[],[],[],[],Ytr0x+alphax*Ytr1,[],Ysh0x+Ysh1*alphax,[],[Vsp2,zeros(size(Vsp2))],[],[],[],[],Tmech1,Varref1,Ef1,Pm1,Eq11); 
        SysParaxx.nIslands=nIslands;
        SysParaxx.islands=islands;
        SysParaxx.refs=refs; 
        
        diff=checkEquationBalanceSynDyn(SysDataxx,SysParaxx,xc,dxc);
        
        if max(abs(diff))<1.05*diffParadigm
            break;
        end
        alphax=alphax*alphaPreMult;
        countCheckAlpha=countCheckAlpha+1;
    end
    if countCheckAlpha>=maxCount;alphax=alphax/alphaPreMult;end
    alpha=alphax;
    
    pqx(:,[4,5])=pqx(:,[4,5])+alpha*pqIncrx(:,1:2);
    if size(pqIncrx,2)>=4
        pqx(:,[4,5])=pqx(:,[4,5])+alpha*alpha*pqIncrx(:,3:4);
    end
    if ~isempty(pvx);pvx(:,4)=pvx(:,4)+alpha*pvIncr;end
    Rzip0x=Rzip0x+alpha*Rzip1;
    Rind0x=Rind0x+alpha*Rind1;
    Reind0x=Reind0x+alpha*Reind1;
    Ysh0x=Ysh0x+alpha*Ysh1;
    Vsp0x=Vsp0x+alpha*VspSq2(:,2);
    Ytr0x=Ytr0x+alpha*Ytr1;
    
    alphaConfirm=alphaConfirm+alpha;
    if useAbsT
        strPre=['T=',num2str(absT+alphaConfirm,'%7.5f'),', '];
    else
        strPre=['Step=',num2str(alphaConfirm,'%7.5f'),', '];
    end
    addLog([strPre,'dt=',num2str(alpha,'%7.5f'),', (maxdiff<',num2str(diffTol),').'],'INFO');
%     addLog(['Fluctuation=',num2str(max(fluctuation)),' rate=',num2str(max(fluctuation)/maxTcoeff),'.'],'INFO');
    
    if alpha==0
        noMoveCnt=noMoveCnt+1;
        absNoMoveCnt=absNoMoveCnt+1;
        addLog(['Step did not move! (MaxDiff=',num2str(max(abs(diff))),')' ],'INFO');
        if diffTol>=diffTolMax
            addLog('Max DiffTol reached and not move, exit!','INFO');
            exitFlag=DynSimFlag.FAIL;
            break;
        end
        
        if noMoveCnt>=maxNoMoveCnt
            zipxx=zipx;
            if ~isempty(zipxx)
                zipxx(:,5:10)=repmat(Rzip0x,1,6).*zipx(:,5:10);
            end
            jac=getJacobianMod(nbus,nline,bus,sw,pvx,pqx,shuntx,line,zipxx,indx,sc,vc,ones(size(ind,1),1));
            jacCond=condest(jac);
            addLog(['Cond of Jacobian=',num2str(jacCond), ', relative cond=',num2str(jacCond/initJacCond)],'INFO');
            if jacCond>50*initJacCond
                addLog('Singular likely, exit!','INFO');
                exitFlag=DynSimFlag.FAIL;
                break;
            end
        end
        
        if allowModelSwitch && useNomoveCtrl && absNoMoveCnt>maxAbsNoMoveCnt && alphaConfirm>0
            addLog('Too many times of no move!','INFO');
            exitFlag=DynSimFlag.NOMOVE_CTRL;
            break;
        end
        
        alphaTol=alphaTol/2;
        if alphaTol<alphaTolOrig/100
            alphaTol=alphaTolOrig/100;
        end
        
        diffTol=max([diffTol*diffMul,max(abs(diff))]);
        if diffTol>diffTolMax
            diffTol=diffTolMax;
        end
        addLog(['Enlarge tol! (Tol=',num2str(diffTol),')'],'INFO');
    else % step moves 
        if alphaTol<alphaTolOrig
            alphaTol=alphaTol*2;
        end
        if alphaTol>alphaTolOrig
            alphaTol=alphaTolOrig;
        end
        
        if useDiffCtrl&&(length(diffList)>1)&&(diffList(2)<diffTolMax/10)&&(max(abs(diff))>diffTolMax/2)
            exitFlag=3;
        end
        
        if useAbsT && isfield(varOpt,'lastSwitchT') && isfield(varOpt,'switchDelay')
            if absT+alphaConfirm>varOpt.lastSwitchT+varOpt.switchDelay
                allowSwitchByDelay=1;
            else
                allowSwitchByDelay=0;
            end
        else
            allowSwitchByDelay=0;
        end
        if allowModelSwitch && allowSwitchByDelay          
            CC=real(V);
            DD=imag(V);
            V2=zeros(size(CC));
            for lvl=1:size(V2,2)
                V2(:,lvl)=sum(CC(:,1:lvl).*CC(:,lvl:-1:1),2)+sum(DD(:,1:lvl).*DD(:,lvl:-1:1),2);
            end
            
            if ~isempty(d)
                xd=d(2:end,:)-repmat(d(1,:),size(d,1)-1,1);
            else
                xd=d;
            end
            if ~isempty(w)
                islandSyn=islands(syn(:,1));
                refSynTag=zeros(size(syn,1),1);
                for isl=1:nIslands
                    synInIsland=find(islandSyn==isl);
                    if ~isempty(synInIsland)
                        refSynTag(islandSyn==isl)=synInIsland(1);
                    end
                end
                xw=w-w(refSynTag,:);
            else
                xw=w;
            end
            [V2a,V2b]=getPadeApproxHelm(V2,0);
            idxNanV2=find(isnan(V2a(:,1))|isnan(V2b(:,1)));V2a(idxNanV2,:)=0;V2b(idxNanV2,:)=0;V2a(idxNanV2,1:2)=V2(idxNanV2,1:2);
            [xda,xdb]=getPadeApproxHelm(xd,0);
            idxNanxd=find(isnan(xda(:,1))|isnan(xdb(:,1)));xda(idxNanxd,:)=0;xdb(idxNanxd,:)=0;xda(idxNanxd,1:2)=xd(idxNanxd,1:2);
            [xwa,xwb]=getPadeApproxHelm(xw,0);
            idxNanxw=find(isnan(xwa(:,1))|isnan(xwb(:,1)));xwa(idxNanxw,:)=0;xwb(idxNanxw,:)=0;xwa(idxNanxw,1:2)=xw(idxNanxw,1:2);
            
            %         heCoeffsToCheck=[xd;w;s;eq1;eq2;ed1;ed2;psid;psiq;Ef;Pm;Vavrm;Vavrr;Vavrf;Vavrref;tgovg;tgovm;Tmech;f;dpg;qplt;vg];
            %         maxTcoeff=alpha;
            %         fluctuations=calcFluctuationHE(heCoeffsToCheck,maxTcoeff);
            %         flucTag=getFluctuationTag(fluctuations,maxTcoeff);
            
            heCoeffsToCheckA=[V2a;xda;wa;sa;eq1a;eq2a;ed1a;ed2a;psida;psiqa;Efa;Pma;Vavrma;Vavrra;Vavrfa;Vavrrefa;tgovga;tgovma;Tmecha;fa;dpga;qplta;vga];
            heCoeffsToCheckB=[V2b;xdb;wb;sb;eq1b;eq2b;ed1b;ed2b;psidb;psiqb;Efb;Pmb;Vavrmb;Vavrrb;Vavrfb;Vavrrefb;tgovgb;tgovmb;Tmechb;fb;dpgb;qpltb;vgb];
            heCoeffsToCheckB=[ones(size(heCoeffsToCheckB,1),1),heCoeffsToCheckB];
            flucTag=getFluctuationTagPade(real(heCoeffsToCheckA),real(heCoeffsToCheckB),0.99*alpha,0.001);
            if flucTag && max(abs(isAgcActiveByIsland.*fc))<1e-4&&allowExit1
                addLog('Dyn model likely reached steady state!(PS)','INFO');
                exitFlag=1;
            else
                exitFlag=0;
            end
            
            if ~isempty(xw) && exitFlag==0
                %             heCoeffsToCheckw=[xd;xw;s;eq1;eq2;ed1;ed2;psid;psiq;Ef;Pm;Vavrm;Vavrr;Vavrf;Vavrref;tgovg;tgovm;Tmech;f;dpg;qplt;vg];
                %             maxTcoeffw=alpha;
                %             fluctuationsw=calcFluctuationHE(heCoeffsToCheckw,maxTcoeffw);
                %             flucTagw=getFluctuationTag(fluctuationsw,maxTcoeffw);
                heCoeffsToCheck=[V2;xw;s;eq1;eq2;ed1;ed2;psid;psiq;Ef;Vavrm;Vavrr;Vavrf;Vavrref;qplt;vg];
                heCoeffsToCheckAw=[V2a;xwa;sa;eq1a;eq2a;ed1a;ed2a;psida;psiqa;Efa;Vavrma;Vavrra;Vavrfa;Vavrrefa;qplta;vga];
                heCoeffsToCheckBw=[V2b;xwb;sb;eq1b;eq2b;ed1b;ed2b;psidb;psiqb;Efb;Vavrmb;Vavrrb;Vavrfb;Vavrrefb;qpltb;vgb];
                heCoeffsToCheckBw=[ones(size(heCoeffsToCheckBw,1),1),heCoeffsToCheckBw];
                fThreshold=0.0001;
                flucTagw=getFluctuationTagPade(real(heCoeffsToCheckAw),real(heCoeffsToCheckBw),0.99*alpha,fThreshold);
                [fluctRateMod,flucTagMod]=getFluctuationTagPadeMod(real(heCoeffsToCheckAw),real(heCoeffsToCheckBw),0.99*alpha,fThreshold);
                [ubc,lbc]=findPolynomialBounds(real(heCoeffsToCheck(:,2:end)),0.99*alpha);
                fluctRatePs=min(abs([ubc,lbc]),[],2);
                if flucTagw
                    addLog('Inter-rotor transients approximately fades away!(PS)','INFO');
                    exitFlag=2;
                else
                    exitFlag=0;
                end
            end
        end
%         heCoeffsToChecka=[xda;wa;sa;eq1a;eq2a;ed1a;ed2a;psida;psiqa;Efa;Pma;Vavrma;Vavrra;Vavrfa;Vavrrefa;tgovga;tgovma;Tmecha;fa;qplta;vga];
%         heCoeffsToCheckb=[xdb;wb;sb;eq1b;eq2b;ed1b;ed2b;psidb;psiqb;Efb;Pmb;Vavrmb;Vavrrb;Vavrfb;Vavrrefb;tgovgb;tgovmb;Tmechb;fb;qpltb;vgb];
%         fluctuationsPade=calcFluctuationHEpade(heCoeffsToChecka,heCoeffsToCheckb,maxTcoeff);
%         flucTagPade=getFluctuationTag(fluctuationsPade,maxTcoeff); 
%         if flucTagPade
%             addLog('Dyn model likely reached steady state!(Pade)','INFO');
%         end 
        
        noMoveCnt=0;
        
        tt=[dAlpha:dAlpha:(alpha-alphaTol/10),alpha]; 
        vcx=polyvalVec(Va(:,end:-1:1),tt)./polyvalVec([Vb(:,end:-1:1),ones(size(Vb,1),1)],tt);
        qcx=polyvalVec(Qa(:,end:-1:1),tt)./polyvalVec([Qb(:,end:-1:1),ones(size(Qb,1),1)],tt);
        qcx(busType==0,:)=0;
        scx=polyvalVec(sa(:,end:-1:1),tt)./polyvalVec([sb(:,end:-1:1),ones(size(sb,1),1)],tt);
        dcx=polyvalVec(da(:,end:-1:1),tt)./polyvalVec([db(:,end:-1:1),ones(size(db,1),1)],tt);
        wcx=polyvalVec(wa(:,end:-1:1),tt)./polyvalVec([wb(:,end:-1:1),ones(size(wb,1),1)],tt);
        eq1cx=polyvalVec(eq1a(:,end:-1:1),tt)./polyvalVec([eq1b(:,end:-1:1),ones(size(eq1b,1),1)],tt);
        eq2cx=polyvalVec(eq2a(:,end:-1:1),tt)./polyvalVec([eq2b(:,end:-1:1),ones(size(eq2b,1),1)],tt);
        ed1cx=polyvalVec(ed1a(:,end:-1:1),tt)./polyvalVec([ed1b(:,end:-1:1),ones(size(ed1b,1),1)],tt);
        ed2cx=polyvalVec(ed2a(:,end:-1:1),tt)./polyvalVec([ed2b(:,end:-1:1),ones(size(ed2b,1),1)],tt);
        psidcx=polyvalVec(psida(:,end:-1:1),tt)./polyvalVec([psidb(:,end:-1:1),ones(size(psidb,1),1)],tt);
        psiqcx=polyvalVec(psiqa(:,end:-1:1),tt)./polyvalVec([psiqb(:,end:-1:1),ones(size(psiqb,1),1)],tt);
                
        Pmcx=polyvalVec(Pma(:,end:-1:1),tt)./polyvalVec([Pmb(:,end:-1:1),ones(size(Pmb,1),1)],tt);
        Efcx=polyvalVec(Efa(:,end:-1:1),tt)./polyvalVec([Efb(:,end:-1:1),ones(size(Efb,1),1)],tt);
        Vavrmcx=polyvalVec(Vavrma(:,end:-1:1),tt)./polyvalVec([Vavrmb(:,end:-1:1),ones(size(Vavrmb,1),1)],tt);
        Vavrrcx=polyvalVec(Vavrra(:,end:-1:1),tt)./polyvalVec([Vavrrb(:,end:-1:1),ones(size(Vavrrb,1),1)],tt);
        Vavrfcx=polyvalVec(Vavrfa(:,end:-1:1),tt)./polyvalVec([Vavrfb(:,end:-1:1),ones(size(Vavrfb,1),1)],tt);
        Vavrrefcx=polyvalVec(Vavrrefa(:,end:-1:1),tt)./polyvalVec([Vavrrefb(:,end:-1:1),ones(size(Vavrrefb,1),1)],tt);
        tgovgcx=polyvalVec(tgovga(:,end:-1:1),tt)./polyvalVec([tgovgb(:,end:-1:1),ones(size(tgovgb,1),1)],tt);
        tgovmcx=polyvalVec(tgovma(:,end:-1:1),tt)./polyvalVec([tgovmb(:,end:-1:1),ones(size(tgovmb,1),1)],tt);
        Tmechcx=polyvalVec(Tmecha(:,end:-1:1),tt)./polyvalVec([Tmechb(:,end:-1:1),ones(size(Tmechb,1),1)],tt);
        
        fcx=polyvalVec(fa(:,end:-1:1),tt)./polyvalVec([fb(:,end:-1:1),ones(size(fb,1),1)],tt);
        dpgcx=polyvalVec(dpga(:,end:-1:1),tt)./polyvalVec([dpgb(:,end:-1:1),ones(size(dpgb,1),1)],tt);
        qpltcx=polyvalVec(qplta(:,end:-1:1),tt)./polyvalVec([qpltb(:,end:-1:1),ones(size(qpltb,1),1)],tt);
        vgcx=polyvalVec(vga(:,end:-1:1),tt)./polyvalVec([vgb(:,end:-1:1),ones(size(vgb,1),1)],tt);
%                
%         load('2bus_x.mat');
%         if isempty(tend);tlastend=0;else tlastend=tend(end);end
%         E=sw(1,4);
%         r=line(1,8);
%         x=line(1,9);
%         tstart=[tstart;tlastend];
%         tend=[tend;tt(end)+tlastend];
%         istart=[istart;((real(Va(2,1))-E)*(real(Va(2,1))-E)+imag(Va(2,1))*imag(Va(2,1)))/(r*r+x*x)];
%         iend=[iend;((real(vcx(2,end))-E)*(real(vcx(2,end))-E)+imag(vcx(2,end))*imag(vcx(2,end)))/(r*r+x*x)];
%         V2=[V2;V(2,:)];
%         save('2bus_x.mat','tstart','tend','istart','iend','V2');
%         
        t=[t,(tt+t(end))];
        VSol=[VSol,vcx];
        QSol=[QSol,qcx];
        sSol=[sSol,scx];
        dSol=[dSol,dcx];
        wSol=[wSol,wcx];
        eq1Sol=[eq1Sol,eq1cx];
        eq2Sol=[eq2Sol,eq2cx];
        ed1Sol=[ed1Sol,ed1cx];
        ed2Sol=[ed2Sol,ed2cx];
        psiqSol=[psiqSol,psiqcx];
        psidSol=[psidSol,psidcx];
        PmSol=[PmSol,Pmcx];
        EfSol=[EfSol,Efcx];
        VavrmSol=[VavrmSol,Vavrmcx];
        VavrrSol=[VavrrSol,Vavrrcx];
        VavrfSol=[VavrfSol,Vavrfcx];
        VavrrefSol=[VavrrefSol,Vavrrefcx];
        tgovgSol=[tgovgSol,tgovgcx];
        tgovmSol=[tgovmSol,tgovmcx];
        TmechSol=[TmechSol,Tmechcx];
        
        fSol=[fSol,fcx];
        dpgSol=[dpgSol,dpgcx];
        qpltSol=[qpltSol,qpltcx];
        vgSol=[vgSol,vgcx];
        
        V0x=vcx(:,end);
        Q0x=qcx(:,end);
        s0x=scx(:,end);
        d0x=dcx(:,end);w0x=wcx(:,end);eq10x=eq1cx(:,end);eq20x=eq2cx(:,end);ed10x=ed1cx(:,end);ed20x=ed2cx(:,end);psid0x=psidcx(:,end);psiq0x=psiqcx(:,end);
        Pm0x=Pmcx(:,end);Ef0x=Efcx(:,end);Vavrm0x=Vavrmcx(:,end);Vavrr0x=Vavrrcx(:,end);Vavrf0x=Vavrfcx(:,end);
        Vavrref0x=Vavrrefcx(:,end);tgovg0x=tgovgcx(:,end);tgovm0x=tgovmcx(:,end);Tmech0x=Tmechcx(:,end);
        f0x=fcx(:,end);dpg0x=dpgcx(:,end);qplt0x=qpltcx(:,end);vg0x=vgcx(:,end);
        
        alphaList=[alphaList,alpha];
        diffList=[diffList,max(abs(diff))];
        
        if allowModelSwitch
            if exitFlag==DynSimFlag.STEADY||exitFlag==DynSimFlag.QSS
                break;
            end
        end
        if exitFlag==DynSimFlag.DIFF_CTRL
            break;
        end
    end
end
if length(diffList)>1;diffList(1)=diffList(2);end
stateCurve=[...
VSol;QSol;sSol;dSol;wSol;eq1Sol;eq2Sol;ed1Sol;ed2Sol;psidSol;psiqSol;
PmSol;EfSol;VavrmSol;VavrrSol;VavrfSol;VavrrefSol;tgovgSol;tgovmSol;TmechSol;fSol;dpgSol;qpltSol;vgSol
];
finalAlpha=alphaConfirm;

if exist([iden,'.mat'],'file')
    delete([iden,'.mat']);
end
end
