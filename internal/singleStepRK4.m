function [xTemp,diff,exitflag]=singleStepRK4(SimData,SysData,SysPara,xt,t,dt,iden,AEmethod)
% Single step of Runge Kutta 4
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
%   xt - Initial system state
%	t - Relative time compared with SysPara starting time
%	dt - time step length
%	iden - Aux identifier
%	AEmethod - Method for solving algebraic equations
%		0 - HE
%		1 - NR
%
% OUTPUT
%	xTemp - State at the new step
%	diff - Error vector
%	exitflag - 
%		0 - Normally exit
%		-1 - Fail
%


[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
[nState,idxs]...
    =getIndexDyn(SysData);
[pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,~,~,~,~,Tmech1,Varref1,Ef1,Pm1,Eq11]=unfoldSysPara(SysPara);

exitflag=0;
V0=xt(idxs.vIdx);
% S1
pqx=pq;pqx(:,[4,5])=pq(:,[4,5])+pqIncr(:,1:2)*t;if size(pqIncr,2)>=4;pqx(:,[4,5])=pqx(:,[4,5])+pqIncr(:,3:4)*t*t;end
pvx=pv;pvx(:,4)=pv(:,4)+pvIncr*t;
indx=ind;if ~isempty(ind);indx(:,15:17)=ind(:,15:17).*repmat(Rind0+Rind1*t,1,3);indx(:,13)=ind(:,13)./(Reind0+Reind1*t);end
zipx=zip;if ~isempty(zip);zipx(:,5:10)=zip(:,5:10).*repmat(Rzip0+Rzip1*t,1,6);end
SysDatax=foldSysData(bus,sw,pvx,pqx,shunt,line,indx,zipx,syn,exc,tg,agc,cac,cluster);
SysData=foldSysData(bus,sw,pvx,pqx,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster);
dx1=integ(SysDatax,SysPara,xt,dt);
xTemp=xt;
xTemp=xTemp+dx1/2;
xTemp=adjustAlgebraic(SysDatax,xTemp);

dpq=pqIncr(:,1:2)*dt/2;if size(pqIncr,2)>=4;dpq=dpq+pqIncr(:,3:4)*((t+dt/2)*(t+dt/2)-t*t);end
SysParax=foldSysPara(dpq,pvIncr*dt/2,...
    Rind0+Rind1*t,Rind1*dt/2,Reind0+Reind1*t,Reind1*dt/2,Rzip0+Rzip1*t,Rzip1*dt/2,Ytr0+Ytr1*t,Ytr1*dt/2,Ysh0+Ysh1*t,Ysh1*dt/2,[VspSq2(:,1)+VspSq2(:,2)*t,VspSq2(:,2)*dt/2],[],[],[],[],...
    Tmech1*dt/2,Varref1*dt/2,Ef1*dt/2,Pm1*dt/2,Eq11*dt/2);
SysParax.iden=iden;
if exist([iden,'.mat'],'file')
    if ~exist('p_amd','var')
        load([iden,'.mat']);
    end
    SysParax.p_amd=p_amd;
else
    SysParax.p_amd=[];
end

if AEmethod==0&&exitflag==0
    try
        [Vx,Wx,Qx,fx,tempState,finalAlpha,alphas,diffRec]=hemMachinePFmultiStageAlgebraic(SimData,SysData,SysParax,xt,xTemp);
        VSol=tempState(idxs.vIdx);
        QSol=tempState(idxs.qIdx);
        fSol=tempState(idxs.fIdx);
    catch ME
        exitflag=-1;
    end
else
    [stateNew,flag,diffRec,loop]=solveAlgebraicNR(SimData,SysData,SysParax,xt,xTemp);
    if flag~=0
        addLog(['NR does not converge!'],'WARN');
        exitflag=-1;
    end
    VSol=stateNew(idxs.vIdx);
    QSol=stateNew(idxs.qIdx);
    fSol=stateNew(idxs.fIdx);
end

xTemp(idxs.vIdx)=VSol;
xTemp(idxs.qIdx)=QSol;
xTemp(idxs.fIdx)=fSol;
xTemp=adjustAlgebraic(SysDatax,xTemp);
% S2
pqx=pq;pqx(:,[4,5])=pq(:,[4,5])+pqIncr(:,1:2)*(t+dt/2);if size(pqIncr,2)>=4;pqx(:,[4,5])=pqx(:,[4,5])+pqIncr(:,3:4)*(t+dt/2)*(t+dt/2);end
pvx=pv;pvx(:,4)=pv(:,4)+pvIncr*(t+dt/2);
indx=ind;if ~isempty(ind);indx(:,15:17)=ind(:,15:17).*repmat(Rind0+Rind1*(t+dt/2),1,3);indx(:,13)=ind(:,13)./(Reind0+Reind1*(t+dt/2));end
zipx=zip;if ~isempty(zip);zipx(:,5:10)=zip(:,5:10).*repmat(Rzip0+Rzip1*(t+dt/2),1,6);end
SysDatax=foldSysData(bus,sw,pvx,pqx,shunt,line,indx,zipx,syn,exc,tg,agc,cac,cluster);
dx2=integ(SysDatax,SysPara,xTemp,dt);

xTemp=xt;
xTemp=xTemp+dx2/2;
xTemp=adjustAlgebraic(SysDatax,xTemp);
dpq=pqIncr(:,1:2)*dt/2;if size(pqIncr,2)>=4;dpq=dpq+pqIncr(:,3:4)*((t+dt)*(t+dt)-(t+dt/2)*(t+dt/2));end
SysParax=foldSysPara(dpq,pvIncr*dt/2,...
    Rind0+Rind1*t,Rind1*dt/2,Reind0+Reind1*t,Reind1*dt/2,Rzip0+Rzip1*t,Rzip1*dt/2,Ytr0+Ytr1*t,Ytr1*dt/2,Ysh0+Ysh1*t,Ysh1*dt/2,[VspSq2(:,1)+VspSq2(:,2)*t,VspSq2(:,2)*dt/2],[],[],[],[],...
    Tmech1*dt/2,Varref1*dt/2,Ef1*dt/2,Pm1*dt/2,Eq11*dt/2);

if AEmethod==0&&exitflag==0
    try
        [Vx,Wx,Qx,fx,tempState,finalAlpha,alphas,diffRec]=hemMachinePFmultiStageAlgebraic(SimData,SysData,SysParax,xt,xTemp);
        VSol=tempState(idxs.vIdx);
        QSol=tempState(idxs.qIdx);
        fSol=tempState(idxs.fIdx);
    catch ME
        exitflag=-1;
    end
else
    [stateNew,flag,diffRec,loop]=solveAlgebraicNR(SimData,SysData,SysParax,xt,xTemp);
    if flag~=0
        addLog(['NR does not converge!'],'WARN');
        exitflag=-1;
    end
    VSol=stateNew(idxs.vIdx);
    QSol=stateNew(idxs.qIdx);
    fSol=stateNew(idxs.fIdx);
end

xTemp(idxs.vIdx)=VSol;
xTemp(idxs.qIdx)=QSol;
xTemp(idxs.fIdx)=fSol;
xTemp=adjustAlgebraic(SysDatax,xTemp);
% S3
pqx=pq;pqx(:,[4,5])=pq(:,[4,5])+pqIncr(:,1:2)*(t+dt);if size(pqIncr,2)>=4;pqx(:,[4,5])=pqx(:,[4,5])+pqIncr(:,3:4)*(t+dt)*(t+dt);end
pvx=pv;pvx(:,4)=pv(:,4)+pvIncr*(t+dt);
indx=ind;if ~isempty(ind);indx(:,15:17)=ind(:,15:17).*repmat(Rind0+Rind1*(t+dt),1,3);indx(:,13)=ind(:,13)./(Reind0+Reind1*(t+dt));end
zipx=zip;if ~isempty(zip);zipx(:,5:10)=zip(:,5:10).*repmat(Rzip0+Rzip1*(t+dt),1,6);end
SysDatax=foldSysData(bus,sw,pvx,pqx,shunt,line,indx,zipx,syn,exc,tg,agc,cac,cluster);
dx3=integ(SysDatax,SysPara,xTemp,dt);

xTemp=xt;
xTemp=xTemp+dx3/2;
xTemp=adjustAlgebraic(SysDatax,xTemp);
dpq=pqIncr(:,1:2)*dt;if size(pqIncr,2)>=4;dpq=dpq+pqIncr(:,3:4)*((t+dt)*(t+dt)-t*t);end
SysParax=foldSysPara(dpq,pvIncr*dt,...
    Rind0+Rind1*t,Rind1*dt,Reind0+Reind1*t,Reind1*dt,Rzip0+Rzip1*t,Rzip1*dt,Ytr0+Ytr1*t,Ytr1*dt,Ysh0+Ysh1*t,Ysh1*dt,[VspSq2(:,1)+VspSq2(:,2)*t,VspSq2(:,2)*dt],[],[],[],[],...
    Tmech1*dt,Varref1*dt,Ef1*dt,Pm1*dt,Eq11*dt);

if AEmethod==0&&exitflag==0
    try
        [Vx,Wx,Qx,fx,tempState,finalAlpha,alphas,diffRec]=hemMachinePFmultiStageAlgebraic(SimData,SysData,SysParax,xt,xTemp);
        VSol=tempState(idxs.vIdx);
        QSol=tempState(idxs.qIdx);
        fSol=tempState(idxs.fIdx);
    catch ME
        exitflag=-1;
    end
else
    [stateNew,flag,diffRec,loop]=solveAlgebraicNR(SimData,SysData,SysParax,xt,xTemp);
    if flag~=0
        addLog(['NR does not converge!'],'WARN');
        exitflag=-1;
    end
    VSol=stateNew(idxs.vIdx);
    QSol=stateNew(idxs.qIdx);
    fSol=stateNew(idxs.fIdx);
end

xTemp(idxs.vIdx)=VSol;
xTemp(idxs.qIdx)=QSol;
xTemp(idxs.fIdx)=fSol;
xTemp=adjustAlgebraic(SysDatax,xTemp);
% S4
pqx=pq;pqx(:,[4,5])=pq(:,[4,5])+pqIncr(:,1:2)*(t+dt);if size(pqIncr,2)>=4;pqx(:,[4,5])=pqx(:,[4,5])+pqIncr(:,3:4)*(t+dt)*(t+dt);end
pvx=pv;pvx(:,4)=pv(:,4)+pvIncr*(t+dt);
indx=ind;if ~isempty(ind);indx(:,15:17)=ind(:,15:17).*repmat(Rind0+Rind1*(t+dt),1,3);indx(:,13)=ind(:,13)./(Reind0+Reind1*(t+dt));end
zipx=zip;if ~isempty(zip);zipx(:,5:10)=zip(:,5:10).*repmat(Rzip0+Rzip1*(t+dt),1,6);end
SysDatax=foldSysData(bus,sw,pvx,pqx,shunt,line,indx,zipx,syn,exc,tg,agc,cac,cluster);
dx4=integ(SysDatax,SysPara,xTemp,dt);

xTemp=xt;
xTemp=xTemp+(dx1+2*dx2+2*dx3+dx4)/6;
xTemp=adjustAlgebraic(SysDatax,xTemp);
dpq=pqIncr(:,1:2)*dt;if size(pqIncr,2)>=4;dpq=dpq+pqIncr(:,3:4)*((t+dt)*(t+dt)-t*t);end
SysParax=foldSysPara(dpq,pvIncr*dt,...
    Rind0+Rind1*t,Rind1*dt,Reind0+Reind1*t,Reind1*dt,Rzip0+Rzip1*t,Rzip1*dt,Ytr0+Ytr1*t,Ytr1*dt,Ysh0+Ysh1*t,Ysh1*dt,[VspSq2(:,1)+VspSq2(:,2)*t,VspSq2(:,2)*dt],[],[],[],[],...
    Tmech1*dt,Varref1*dt,Ef1*dt,Pm1*dt,Eq11*dt);

if AEmethod==0&&exitflag==0
    try
        [Vx,Wx,Qx,fx,tempState,finalAlpha,alphas,diffRec]=hemMachinePFmultiStageAlgebraic(SimData,SysData,SysParax,xt,xTemp);
        VSol=tempState(idxs.vIdx);
        QSol=tempState(idxs.qIdx);
        fSol=tempState(idxs.fIdx);
    catch ME
        exitflag=-1;
    end
else
    [stateNew,flag,diffRec,loop]=solveAlgebraicNR(SimData,SysData,SysParax,xt,xTemp);
    if flag~=0
        addLog(['NR does not converge!'],'WARN');
        exitflag=-1;
    end
    VSol=stateNew(idxs.vIdx);
    QSol=stateNew(idxs.qIdx);
    fSol=stateNew(idxs.fIdx);
end

xTemp(idxs.vIdx)=VSol;
xTemp(idxs.qIdx)=QSol;
xTemp(idxs.fIdx)=fSol;
xTemp=adjustAlgebraic(SysDatax,xTemp);

%         dxc=integ(SysDatax,SysPara,xTemp,dt)/dt;
dxc=dx4/dt;
SysParaxx=foldSysPara([],[],[],[],[],[],[],[],Ytr0+Ytr1*(t+dt),[],Ysh0+Ysh1*(t+dt),[],[VspSq2(:,1)+VspSq2(:,2)*(t+dt),zeros(size(VspSq2,1),1)],[],[],[],[],Tmech1,Varref1,Ef1,Pm1,Eq11);
diff=checkEquationBalanceSynDyn(SysDatax,SysParaxx,xTemp,dxc);


end