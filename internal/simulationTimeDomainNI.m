function [stateCurve,tt,diffList,nxtDt,exitflag]=simulationTimeDomainNI(SimData,SysData,SysPara,x0)
% Numerical integration approach for dynamic simulation
%
% FUNCTION simulationTimeDomainNI
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
%   stateCurve - A list of states in the order of time
%   tt - A list of time points (starting with 0)
%   diffList - A list of errors
%   nxtDt - ending time step 
%   exitFlag - 
%       0  - Success and normally exit
%       -1 - Fail to finish (due to computation errors or failures)
%
% MODES (controlled by SimData.method)
%   SimData.method has format X.YZ
%   X - Method for solving differential equations (DE)
%       0 - HE
%       1 - Modified Euler
%       2 - RK4
%       3 - Trapezoidal
%
%   Y - Method for solving algebraic equations (AE)
%       0 - HE
%       1 - NR
%
%   Z - Variable step (only applies to X>1)
%       0 - Fixed step
%       1 - Adaptive step
%
%%
[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
[nState,idxs]...
    =getIndexDyn(SysData);
% [V0,Q0,s0,d0,w0,eq10,eq20,ed10,ed20,psid0,psiq0,Pm0,Ef0,Vavrm0,Vavrr0,Vavrf0,Vavrref0,tgovg0,tgovm0,tgovmech0]=unfoldX(x0,SysData);
% [~,~,sNew,dNew,~,eq1New,eq2New,ed1New,ed2New,psidNew,psiqNew,~,EfNew,~,~,~,~,~,~,~]=unfoldX(xNew,SysData);
[maxTime,segTime,dt,nlvl,taylorN,alphaTol,diffTol,diffTolMax,method,varOpt]=unfoldSimData(SimData);
[pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,~,~,~,~,Tmech1,Varref1,Ef1,Pm1,Eq11]=unfoldSysPara(SysPara);

if ~isempty(varOpt) && isfield(varOpt,'absT')
    useAbsT=1;
    absT=varOpt.absT;
else
    useAbsT=0;
end

[~,~,~,~,~,~,~,~,method]=unfoldSimData(SimData);
DEmethod=floor(method);
AEmethod=round(10*mod(method,1));
varStep=round(100*mod(method,0.1));

eval('setting_func_simulationTimeDomainNI');
if varStep~=0
% some temporary parameters
dt=min([maxStep,max([minStep,dt])]);
%
end

tt=0;
diffList=0;

% tt=0:dt:maxTime;
% if maxTime>tt(end)+1e-9;tt=[tt,maxTime];end
% if length(tt)<2;tt=[0,maxTime];end
% diffList=zeros(size(tt));

nbus=size(bus,1);
nline=size(line,1);
nInd=size(ind,1);
nZip=size(zip,1);
nSyn=size(syn,1);
nTg=size(tg,1);
nExc=size(exc,1);

if isempty(pqIncr);pqIncr=zeros(size(pq,1),2);end
if isempty(pvIncr);pvIncr=zeros(size(pv,1),1);end
if isempty(Rind0);Rind0=ones(nInd,1);end
if isempty(Rind1);Rind1=zeros(nInd,1);end
if isempty(Reind0);Reind0=ones(nInd,1);end
if isempty(Reind1);Reind1=zeros(nInd,1);end
if isempty(Rzip0);Rzip0=ones(nZip,1);end
if isempty(Rzip1);Rzip1=zeros(nZip,1);end
if isempty(Ytr0)||isempty(Ytr1)||isempty(Ysh0)||isempty(Ysh1)
    [~,Ytr0x,Ysh0x,~,~,~,~]=getYMatrix(nbus,line);
    if isempty(Ytr0);Ytr0=Ytr0x;end
    if isempty(Ytr1);Ytr1=0*Ytr0;end
    if isempty(Ysh0);Ysh0=Ysh0x;end
    if isempty(Ysh1);Ysh1=0*Ysh0;end
end
if isempty(Tmech1);Tmech1=zeros(nTg,1);end
if isempty(Varref1);Varref1=zeros(nExc,1);end
if isempty(Ef1);Ef1=zeros(nSyn,1);end
if isempty(Pm1);Pm1=zeros(nSyn,1);end
if isempty(Eq11);Eq11=zeros(nSyn,1);end

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

yShunt=zeros(nbus,1);
yShunt(shunt(:,1))=shunt(:,5)+1j*shunt(:,6);
Ysh0=Ysh0+yShunt;
shunt=zeros(0,7);
if isempty(VspSq2)
    V0=x0(idxs.vIdx);
    VspSq2=[abs(V0).*abs(V0),zeros(nbus,1)];
end

segAlpha=min([maxTime,segTime]);

maxLvlHe=3;
if nlvl>maxLvlHe
    nlvl=maxLvlHe;
end
SimDatax=foldSimData(maxTime,segAlpha,dt,nlvl,taylorN,alphaTol,diffTol,diffTolMax,method);

stateCurve=zeros(size(x0,1),1);
stateCurve(:,1)=x0;

iden=[datestr(clock,30),'_dyn'];

%
% [t,stateCurve,finalAlpha,alphaList,diffList]=hemMachinePFmultiStageDyn(SimDatax,SysData,SysParax,x0);
V0=x0(idxs.vIdx);
SysDatax=foldSysData(bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster);
dxc=integ(SysDatax,SysPara,x0,dt)/dt;
SysParaxx=foldSysPara(pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,[],Ysh0,[],[abs(V0).*abs(V0),zeros(nbus,1)],[],[],[],[],Tmech1,Varref1,Ef1,Pm1,Eq11);
diff=checkEquationBalanceSynDyn(SysDatax,SysParaxx,x0,dxc);
diffList(1)=max(abs(diff));

it=1;
stepCnt=0;
errPtoX=zeros(0,1);
nxtDt=dt;
correctRatio=1;
while tt(end)+dt<=maxTime+stepSizeEps
    xt=stateCurve(:,end);
    nxtDt=dt;
%     dt=tt(it)-tt(it-1);
    if DEmethod==1 % Improved Euler
        SysParaOrig=foldSysPara(pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,[],[],[],[],Tmech1,Varref1,Ef1,Pm1,Eq11);
        SysDataOrig=foldSysData(bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster);
        
        [xTemp,diff,exitflag]=singleStepME(SimDatax,SysDataOrig,SysParaOrig,xt,tt(end),dt,iden,AEmethod);
        
        if exitflag~=0
            if varStep
                dt=max([minStep,dt/10]);
                [xTemp,diff,exitflag]=singleStepME(SimDatax,SysDataOrig,SysParaOrig,xt,tt(end),dt,iden,AEmethod);
                if exitflag~=0
                    break;
                end
            else
                break;                
            end
        end

        stateCurve=[stateCurve,xTemp];        
        diffList=[diffList,max(abs(diff))];        
        
        if varStep&&stepCnt>=stepCntThreshold            
            [xTempMid,~,tmpExitflag]=singleStepME(SimDatax,SysDataOrig,SysParaOrig,xt,tt(end),dt/2,iden,AEmethod);
            if tmpExitflag==0
                [xTemp1,~,tmpExitflag]=singleStepME(SimDatax,SysDataOrig,SysParaOrig,xTempMid,tt(end)+dt/2,dt/2,iden,AEmethod);
            end
            
            if tmpExitflag==0
                errX=max(abs(xTemp1-xTemp));
                maxErrX=max(abs(errX));
                errPtoX=[errPtoX;maxErrX/diffList(end)];
                stepRatio=min([max([(diffTol*mean(errPtoX)/2/maxErrX)^(1/3),0.3]),2]);
                nxtDt=correctRatio^1.5*0.95*dt*stepRatio;
                correctRatio=sqrt(correctRatio);
                disp('');
            end
            
            stepCnt=0;
        end
        
    elseif DEmethod==2% R-K 4
        SysParaOrig=foldSysPara(pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,[],[],[],[],Tmech1,Varref1,Ef1,Pm1,Eq11);
        SysDataOrig=foldSysData(bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster);
        
        [xTemp,diff,exitflag]=singleStepRK4(SimDatax,SysDataOrig,SysParaOrig,xt,tt(end),dt,iden,AEmethod);
        
        if exitflag~=0
            if varStep
                dt=max([minStep,dt/10]);
                [xTemp,diff,exitflag]=singleStepRK4(SimDatax,SysDataOrig,SysParaOrig,xt,tt(end),dt,iden,AEmethod);
                if exitflag~=0
                    break;
                end
            else
                break;                
            end
        end

        stateCurve=[stateCurve,xTemp];        
        diffList=[diffList,max(abs(diff))];        
        
        if varStep&&stepCnt>=stepCntThreshold            
            [xTempMid,~,tmpExitflag]=singleStepRK4(SimDatax,SysDataOrig,SysParaOrig,xt,tt(end),dt/2,iden,AEmethod);
            if tmpExitflag==0
                [xTemp1,~,tmpExitflag]=singleStepRK4(SimDatax,SysDataOrig,SysParaOrig,xTempMid,tt(end)+dt/2,dt/2,iden,AEmethod);
            end
            
            if tmpExitflag==0
                errX=max(abs(xTemp1-xTemp));
                maxErrX=max(abs(errX));
                errPtoX=[errPtoX;maxErrX/diffList(end)];
                stepRatio=min([max([(diffTol*mean(errPtoX)/2/maxErrX)^0.2,0.3]),2]);
                nxtDt=correctRatio^1.5*0.95*dt*stepRatio;
                correctRatio=sqrt(correctRatio);
                disp('');
            end
            
            stepCnt=0;
        end

    else %trapezoidal rule
        SysParaOrig=foldSysPara(pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,[],[],[],[],Tmech1,Varref1,Ef1,Pm1,Eq11);
        SysDataOrig=foldSysData(bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster);
        
        [xTemp,diff,exitflag]=singleStepTRAP(SimDatax,SysDataOrig,SysParaOrig,xt,tt(end),dt,iden,AEmethod);
        
        if exitflag~=0
            if varStep
                dt=max([minStep,dt/10]);
                [xTemp,diff,exitflag]=singleStepTRAP(SimDatax,SysDataOrig,SysParaOrig,xt,tt(end),dt,iden,AEmethod);
                if exitflag~=0
                    break;
                end
            else
                break;                
            end
        end

        stateCurve=[stateCurve,xTemp];        
        diffList=[diffList,max(abs(diff))];        
        
        if varStep&&stepCnt>=stepCntThreshold            
            [xTempMid,~,tmpExitflag]=singleStepTRAP(SimDatax,SysDataOrig,SysParaOrig,xt,tt(end),dt/2,iden,AEmethod);
            if tmpExitflag==0
                [xTemp1,~,tmpExitflag]=singleStepTRAP(SimDatax,SysDataOrig,SysParaOrig,xTempMid,tt(end)+dt/2,dt/2,iden,AEmethod);
            end
            
            if tmpExitflag==0
                errX=max(abs(xTemp1-xTemp));
                maxErrX=max(abs(errX));
                errPtoX=[errPtoX;maxErrX/diffList(end)];
                stepRatio=min([max([(diffTol*mean(errPtoX)/2/maxErrX)^(1/3),0.3]),2]);
                nxtDt=correctRatio^1.9*0.99*dt*stepRatio;
                correctRatio=sqrt(correctRatio);
                disp('');
            end
            
            stepCnt=0;
        end 

    end
        
    tt=[tt,tt(end)+dt];
    it=it+1;
    
    if useAbsT
        strPre=['T=',num2str(absT+tt(it),'%7.5f'),', '];
    else
        strPre=['TD time (in stage)=',num2str(tt(it),'%7.5f'),', '];
    end    
    addLog([strPre,'(s), dt=',num2str(dt,'%7.5f'),',cr=',num2str(correctRatio,'%7.5f'),',nxt=',num2str(nxtDt,'%7.5f')],'INFO');
    
    if tt(end)>=maxTime-stepSizeEps
        break;
    end
    
    if varStep
        stepCnt=stepCnt+1;
        dt=nxtDt;

        dt=min([maxStep,max([minStep,dt])]);

        correctRatio=correctRatio*dt/nxtDt;

        if mod(it,stepCntThreshold*20)==0
            correctRatio=1;
        end
    end
    
    if tt(end)+dt>=maxTime-stepSizeEps
        dt=maxTime-tt(end);
    end
end

stateCurve=stateCurve(:,1:it);
tt=tt(1:it);
diffList=diffList(1:it);

if exist([iden,'.mat'],'file')
    delete([iden,'.mat']);
end
end
