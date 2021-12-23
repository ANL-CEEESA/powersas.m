function [stateNew,finalAlpha,alphaList,diff,SysDataUpd,maps]=...
    simulationAddConnectionLine(SysData,SysDataBase,SysPara,SimData,maps,x0,lineAddIdx)
% Simulate add line (onto the same island)
%
% FUNCTION simulationAddConnectionLine
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
%   SysData - System data for simulation
%   SysDataBase - Base system data for simulation
%   SysPara - Parameters representing the events happening in the system
%   SimData - Simulation parameters
%	*Map - The mapping of current data to base data
%	inv*Map - Inverse mapping, base data to current data
%   x0 - Initial system state
%   lineAddIdx - Index of added line
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

lineAdd=lineBase(lineAddIdx,:);
[~,Ytr0,Ysh0,~,~,~,~]=getYMatrix(nbus,line);
lineAdd(:,1)=busMap(lineAdd(:,1));
lineAdd(:,2)=busMap(lineAdd(:,2));
lineUpd=[line;lineAdd];

lineMap(lineAddIdx)=(nline+1):(nline+size(lineAdd,1));
invLineMap=[invLineMap;lineAddIdx];

[~,Ytr0New,Ysh0New,~,~,~,~]=getYMatrix(nbus,lineUpd,fault);

SysPara=foldSysPara(pqIncrS,pvIncrS,Rind0S,Rind1S,[],[],Rzip0S,Rzip1S,Ytr0,Ytr0New-Ytr0,Ysh0,Ysh0New-Ysh0,[],[],[],[],[]);
xNew=x0;

if AEmethod==0
    [stateNew,finalAlpha,alphaList,diff]=solveAlgebraicHem(SimData,SysData,SysPara,x0,xNew);
else
    [stateNew,flag,diff,loop]=solveAlgebraicNR(SimData,SysData,SysPara,x0,xNew);
    if flag==0
        finalAlpha=1;
    else
        finalAlpha=0;
    end
    alphaList=linspace(0,finalAlpha,loop);
end

SysDataUpd=foldSysData(bus,sw,pv,pq,shunt,lineUpd,ind,zip,syn,exc,tg,agc,cac,cluster);

maps.busMap=busMap;maps.invBusMap=invBusMap;
maps.indMap=indMap;maps.invIndMap=invIndMap;
maps.synMap=synMap;maps.invSynMap=invSynMap;
maps.lineMap=lineMap;maps.invLineMap=invLineMap;

end