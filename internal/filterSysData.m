function [SysDataUpd,busFtr,swFtr,pvFtr,pqFtr,shuntFtr,lineFtr,indFtr,zipFtr,synFtr,excFtr,tgFtr,...
    busMapNew,indMapNew,synMapNew,lineMapNew,invBusMapNew,invIndMapNew,invSynMapNew,invLineMapNew]=...
    filterSysData(busTag,SysData,busMap,indMap,synMap,lineMap,invBusMap,invIndMap,invSynMap,invLineMap)
% Generate filter mappings of system data when partial system is studied and filter system data
%
% FUNCTION filterSysData
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
%   busTag - 0/1 tags of the existence of buses
%   SysData - System data
%	*Map - The mapping of current data to base data
%	inv*Map - Inverse mapping, base data to current data
%
% OUTPUT
%	SysDataUpd - Updated system data
%	*Ftr - Filter of current data to rip off inactive components
%
[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
busFtr=find(busTag==1);
swFtr=find(busTag(sw(:,1))==1);
pvFtr=find(busTag(pv(:,1))==1);
pqFtr=find(busTag(pq(:,1))==1);
shuntFtr=find(busTag(shunt(:,1))==1);
lineFtr=find(busTag(line(:,1))==1&busTag(line(:,2))==1);
indFtr=find(busTag(ind(:,1))==1);
zipFtr=find(busTag(zip(:,1))==1);
synTag=busTag(syn(:,1));
synFtr=find(synTag==1);
excFtr=find(synTag(exc(:,1))==1);
tgFtr=find(synTag(tg(:,1))==1);
agcFtr=find(busTag(agc(:,1))==1);

invBusFtr=zeros(size(bus,1),1);
invBusFtr(busFtr)=1:length(busFtr);

invSynFtr=zeros(size(syn,1),1);
invSynFtr(synFtr)=1:length(synFtr);

invBusMapNew=invBusMap(busFtr);
invIndMapNew=invIndMap(indFtr);
invSynMapNew=invSynMap(synFtr);
invLineMapNew=invLineMap(lineFtr);

busMapNew=zeros(size(busMap));
indMapNew=zeros(size(indMap));
synMapNew=zeros(size(synMap));
lineMapNew=zeros(size(lineMap));

busMapNew(invBusMapNew)=1:length(invBusMapNew);
indMapNew(invIndMapNew)=1:length(invIndMapNew);
synMapNew(invSynMapNew)=1:length(invSynMapNew);
lineMapNew(invLineMapNew)=1:length(invLineMapNew);

busNew=bus(busFtr,:);
swNew=sw(swFtr,:);
pvNew=pv(pvFtr,:);
pqNew=pq(pqFtr,:);
shuntNew=shunt(shuntFtr,:);
lineNew=line(lineFtr,:);
indNew=ind(indFtr,:);
zipNew=zip(zipFtr,:);
synNew=syn(synFtr,:);
excNew=exc(excFtr,:);
tgNew=tg(tgFtr,:);
agcNew=agc(agcFtr,:);

busNew(:,1)=invBusFtr(busNew(:,1));
swNew(:,1)=invBusFtr(swNew(:,1));
pvNew(:,1)=invBusFtr(pvNew(:,1));
pqNew(:,1)=invBusFtr(pqNew(:,1));
shuntNew(:,1)=invBusFtr(shuntNew(:,1));
lineNew(:,1:2)=invBusFtr(lineNew(:,1:2));
indNew(:,1)=invBusFtr(indNew(:,1));
zipNew(:,1)=invBusFtr(zipNew(:,1));
synNew(:,1)=invBusFtr(synNew(:,1));
excNew(:,1)=invSynFtr(excNew(:,1));
tgNew(:,1)=invSynFtr(tgNew(:,1));
agcNew(:,1)=invBusFtr(agcNew(:,1));

SysDataUpd=foldSysData(busNew,swNew,pvNew,pqNew,shuntNew,lineNew,indNew,zipNew,synNew,excNew,tgNew,agcNew,cac,cluster);
end