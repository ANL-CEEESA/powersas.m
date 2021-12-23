function [stateNew,finalAlpha,alphaList,diff,SysDataUpd,SysParaUpd,maps]=...
    simulationCutLine(busTag,SysData,SysDataBase,SysPara,SimData,maps,x0,lineCutIdx)
% Simulate cut lines
%
% FUNCTION simulationCutLine
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
%   busTag - Tags for determining the partial system that will be kept 
%   SysData - System data for simulation
%   SysDataBase - Base system data for simulation
%   SysPara - Parameters representing the events happening in the system
%   SimData - Simulation parameters
%	*Map - The mapping of current data to base data
%	inv*Map - Inverse mapping, base data to current data
%   x0 - Initial system state
%   lineCutIdx - Index of the line to cut
%
% OUTPUT
%   stateNew - Solved new state
%   finalAlpha - The ending alpha
%   alphaList - A list of alphas
%   diff - Equation mismatches
%   SysDataUpd - Updated system data 
%	*Map - The mapping of current data to base data (updated)
%	inv*Map - Inverse mapping, base data to current data (updated)
%
[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
[busBase,swBase,pvBase,pqBase,shuntBase,lineBase,indBase,zipBase,synBase,excBase,tgBase,agcBase,cacBase,clusterBase]=unfoldSysData(SysDataBase);
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,fault]=unfoldSysPara(SysPara);

busMap=maps.busMap;invBusMap=maps.invBusMap;
indMap=maps.indMap;invIndMap=maps.invIndMap;
synMap=maps.synMap;invSynMap=maps.invSynMap;
lineMap=maps.lineMap;invLineMap=maps.invLineMap;

[~,~,~,~,~,~,~,~,method]=unfoldSimData(SimData);
AEmethod=round(10*mod(method,1));
 
nline=size(line,1);

pqIncrS=zeros(size(pq,1),2);pvIncrS=zeros(size(pv,1),1);
nbus=size(bus,1);
nZip=size(zip,1);
nInd=size(ind,1);
Rzip0S=ones(nZip,1);Rzip1S=zeros(nZip,1);
Rind0S=ones(nInd,1);Rind1S=zeros(nInd,1);

lineCutIdxNew=lineMap(lineCutIdx);
lineCutIdxNew=lineCutIdxNew(lineCutIdxNew~=0);

lineTag=ones(nline,1);
lineTag(lineCutIdxNew)=0;

invLineMapNew=invLineMap(lineTag==1);
lineUpd=line(lineTag==1,:);
lineMapNew=zeros(size(lineBase,1),1);
lineMapNew(invLineMapNew)=1:length(invLineMapNew);

[~,Ytr0,Ysh0,ytrfr,ytrto,yshfr,yshto]=getYMatrix(nbus,line,fault);
[~,Ytr0New,Ysh0New,~,~,~,~]=getYMatrix(nbus,lineUpd,fault);

lineMap=lineMapNew;
invLineMap=invLineMapNew;

SysParaSim=foldSysPara(pqIncrS,pvIncrS,Rind0S,Rind1S,[],[],Rzip0S,Rzip1S,Ytr0,Ytr0New-Ytr0,Ysh0,Ysh0New-Ysh0,[],[],[],[],[]);
SysDataUpd=foldSysData(bus,sw,pv,pq,shunt,lineUpd,ind,zip,syn,exc,tg,agc,cac,cluster);

if ~isempty(busTag)&&~isempty(find(busTag==0,1))
    % Special for cutting lines, need to modify the Y matrix
    [V,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~]=unfoldX(x0,SysData);
    VCutFr=V(line(lineTag==0,1));
    VCutTo=V(line(lineTag==0,2));
    ICutFr=ytrfr(lineTag==0).*(VCutFr-VCutTo);
    ICutTo=ytrto(lineTag==0).*(VCutTo-VCutFr);    
    YshCutPlusFr=ICutFr./VCutFr;
    YshCutPlusTo=ICutTo./VCutTo;    
    YshCutFr=yshfr(lineTag==0);
    YshCutTo=yshto(lineTag==0);
    
    Ysh0Mod=accumarray(line(lineTag==0,1),YshCutPlusFr+YshCutFr,[nbus,1])+...
        accumarray(line(lineTag==0,2),YshCutPlusTo+YshCutTo,[nbus,1]);
    
    SysParaSim=foldSysPara(pqIncrS,pvIncrS,Rind0S,Rind1S,[],[],Rzip0S,Rzip1S,Ytr0New,sparse(nbus,nbus),Ysh0New+Ysh0Mod,-Ysh0Mod,[],[],[],[],[]);
    
    % Purge system
    [SysDataUpdFtr,busFtr,swFtr,pvFtr,pqFtr,shuntFtr,lineFtr,indFtr,zipFtr,synFtr,excFtr,tgFtr,...
        busMap,indMap,synMap,lineMap,invBusMap,invIndMap,invSynMap,invLineMap]=...
        filterSysData(busTag,SysDataUpd,busMap,indMap,synMap,lineMap,invBusMap,invIndMap,invSynMap,invLineMap);
    SysParaSim=filterSysPara(SysParaSim,busFtr,swFtr,pvFtr,pqFtr,shuntFtr,lineFtr,indFtr,zipFtr,synFtr,excFtr,tgFtr);
    SysPara=filterSysPara(SysPara,busFtr,swFtr,pvFtr,pqFtr,shuntFtr,lineFtr,indFtr,zipFtr,synFtr,excFtr,tgFtr);
    x0=filterX(SysDataUpd,SysDataUpdFtr,x0,busFtr,swFtr,pvFtr,pqFtr,shuntFtr,lineFtr,indFtr,zipFtr,synFtr,excFtr,tgFtr);
    SysDataUpd=SysDataUpdFtr;
end

xNew=x0;

if AEmethod==0
    [stateNew,finalAlpha,alphaList,diff]=solveAlgebraicHem(SimData,SysDataUpd,SysParaSim,x0,xNew);
else
    [stateNew,flag,diff,loop]=solveAlgebraicNR(SimData,SysDataUpd,SysParaSim,x0,xNew);
    if flag==0
        finalAlpha=1;
    else
        finalAlpha=0;
    end
    alphaList=linspace(0,finalAlpha,loop);
end

SysParaUpd=SysPara;

maps.busMap=busMap;maps.invBusMap=invBusMap;
maps.indMap=indMap;maps.invIndMap=invIndMap;
maps.synMap=synMap;maps.invSynMap=invSynMap;
maps.lineMap=lineMap;maps.invLineMap=invLineMap;

end