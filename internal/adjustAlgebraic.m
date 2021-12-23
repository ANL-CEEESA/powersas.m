function xt=adjustAlgebraic(SysData,xt)
% [INTERNAL] Adjust the algebraic variables (usually after numerical integration)
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
% FUNCTION adjustAlgebraic
% CAUTION: this function is used internally. Do not call it directly unless you know what you are doing.
%
% INPUT
%   SysData - System data for simulation
%	xt - Initial state
%
% OUTPUT
%	xt - System state after the algebraic variables adjusted
%
[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
[nState,idxs]...
    =getIndexDyn(SysData);
[V,Q,s,d,w,eq1,eq2,ed1,ed2,psid,psiq,Pg,Ef,Vavrm,Vavrr,Vavrf,Vavrref,tgovg,tgovm,tgovmech,f,dpg,qplt,vg]=unfoldX(xt,SysData);
%
nbus=size(bus,1);
Vmag=abs(V);
[nIslands,islands,refs]=searchIslands(bus(:,1),line(:,1:2));

nSyn=size(syn,1);
if ~isempty(exc)
    excIdx=exc(:,1);
    VavrMax=exc(:,3);
    VavrMin=exc(:,4);
    muavr0=exc(:,5);
    Tavr1=exc(:,7);
    Tavr2=exc(:,6);
    vavrf0=exc(:,8);
    Vavr0=exc(:,9);
    Tavre=exc(:,10);
    Tavrr=exc(:,11);
    
    VmagAvr=Vmag(excIdx);
    
    Efx=Vavrf;
    tavrMaxDiff=Vavrf-VavrMax;
    tavrMinDiff=Vavrf-VavrMin;
    Efx(tavrMaxDiff>0)=VavrMax(tavrMaxDiff>0);
    Efx(tavrMinDiff<0)=VavrMin(tavrMinDiff<0);
    
    xt(idxs.efIdx(excIdx))=Efx;
end

if ~isempty(tg)
    tgIdx=tg(:,1);
    wtgref=tg(:,3);
    Rtg=tg(:,4);
    Ttgmax=tg(:,5);
    Ttgmin=tg(:,6);
    Ttg2=tg(:,7);
    Ttg1=tg(:,8);
    
    xt(idxs.tgovmIdx)=(tgovg+Ttg1./Ttg2.*(wtgref-w(tgIdx))./Rtg+tgovmech);
    
    Pmx=tgovm;
    tgovMaxDiff=tgovm(:,1)-Ttgmax;
    tgovMinDiff=tgovm(:,1)-Ttgmin;
    Pmx(tgovMaxDiff>0)=Ttgmax(tgovMaxDiff>0);
    Pmx(tgovMinDiff<0)=Ttgmin(tgovMinDiff<0);
    
    xt(idxs.pgIdx(tgIdx))=Pmx;
end

synTag=zeros(nbus,1);
synTag(syn(:,1))=1:nSyn;
numSynOnBus=accumarray(syn(:,1),1,[nbus,1]);
dpgTag=ones(nbus,1);
for islIdx=1:nIslands
    busIsland=find(islands==islIdx);
    synTagIsland=synTag(busIsland);
    wIsland=w(synTagIsland(synTagIsland~=0),1);
    if ~isempty(wIsland)
        xt(idxs.fIdx(busIsland))=mean(wIsland); % note that here the freq can be different
        dpgTag(busIsland)=0;
    end
end

end
