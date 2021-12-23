function [sysDataGross,grossMaps,sysMaps,newSysDataList,originalSw]=agglomerateSystems(SysDataList,links)
% Agglomerate a cell of system data as one
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
%	SysDataList - A cell of SysData
%	links
% 		Definition of links: cell of link
% 		Definition of link: link-to-area (usually 1), self-busNum, to-busNum, serial-R, serial X
% 		Definition of link: link-to-area (usually 1), 0          , 0        , swP, swQ
%
% OUTPUT
%	sysDataGross - Agglomerated system data
%	grossMaps - Mapping from agglomerated system to subsystem
%	sysMaps - Mapping from subsystem to agglomerated system
%	newSysDataList - Cell of subsystems with external links
%

sysMaps=cell(length(SysDataList),1);
newSysDataList=cell(length(SysDataList),1);
originalSw=cell(length(SysDataList),1);

[busMain,swMain,pvMain,pqMain,shuntMain,lineMain,indMain,zipMain,synMain,excMain,tgMain,agcMain,cacMain,clusterMain]=unfoldSysData(SysDataList{1});
auxBusMain=zeros(size(busMain,1),1);

grossMaps.busAreaMap=ones(size(busMain,1),1);
grossMaps.swAreaMap=ones(size(swMain,1),1);
grossMaps.pvAreaMap=ones(size(pvMain,1),1);
grossMaps.pqAreaMap=ones(size(pqMain,1),1);
grossMaps.shuntAreaMap=ones(size(shuntMain,1),1);
grossMaps.lineAreaMap=ones(size(lineMain,1),1);
grossMaps.indAreaMap=ones(size(indMain,1),1);
grossMaps.zipAreaMap=ones(size(zipMain,1),1);
grossMaps.synAreaMap=ones(size(synMain,1),1);
grossMaps.excAreaMap=ones(size(excMain,1),1);
grossMaps.tgAreaMap=ones(size(tgMain,1),1);
grossMaps.agcAreaMap=ones(size(agcMain,1),1);
grossMaps.cacAreaMap=ones(size(cacMain,1),1);
grossMaps.clusterAreaMap=ones(size(clusterMain,1),1);

m.busMap=1:size(busMain,1);
m.swMap=1:size(swMain,1);
m.pvMap=1:size(pvMain,1);
m.pqMap=1:size(pqMain,1);
m.lineMap=1:size(lineMain,1);
m.shuntMap=1:size(shuntMain,1);
m.indMap=1:size(indMain,1);
m.zipMap=1:size(zipMain,1);
m.synMap=1:size(synMain,1);
m.excMap=1:size(excMain,1);
m.tgMap=1:size(tgMain,1);
m.agcMap=1:size(agcMain,1);
m.cacMap=1:size(cacMain,1);
m.clusterMap=1:size(clusterMain,1);
m.extBusMap=zeros(0,4);

sysMaps{1}=m;
newSysDataList{1}=SysDataList{1};
nSys=length(SysDataList);

for iSys=2:nSys
    % Warning: currently assume all the links to area 1 (main).
    link=links{iSys};  
    iSwConfig=find(link(:,2)==0,1);
    if ~isempty(iSwConfig); sSw=link(iSwConfig,[4,5]);else sSw=[0,0];end
    
    [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysDataList{iSys});
    
    swInPq=find(pq(:,1)==sw(1,1));
    if ~isempty(swInPq)
        pq(swInPq,[4,5])=pq(swInPq,[4,5])-sSw;
    else
        pq=[pq;[sw(1,1),bus(sw(1,1),2),100.0,-sSw,2.0,0.0,0,1]];
    end
%     swInPv=find(pv(:,1)==sw(1,1));
%     if ~isempty(swInPv)
%         pv(swInPv,4)=pv(swInPv,4)+sSw(1);
%     else
%         pv=[pv;[sw(1,1),bus(sw(1,1),2),100.0,sSw(1),sw(1,4),sw(1,[6,7]),2.0,0.0,0,1]];
%     end
    originalSw{iSys}=sw(1,:);
    sw=zeros(0,13);
    linkLine=link(link(:,2)~=0|link(:,3)~=0,:);
    
    nbusx=size(bus,1);
    busesInMain=busMain(linkLine(:,3));
    auxBusMain(:)=0;auxBusMain(busesInMain)=1;
    busesInMain=find(auxBusMain==1);
    auxBusMain(busesInMain)=(nbusx+1):(nbusx+size(busesInMain,1));
    busx=[bus;[((nbusx+1):(nbusx+size(busesInMain,1)))',busMain(busesInMain,2:end)]];
    linex=[line;[auxBusMain(linkLine(:,3)),    linkLine(:,2), repmat([100.00],size(linkLine,1),1),  bus(linkLine(:,2),2) , repmat([60, 0,   0.0000],size(linkLine,1),1),  linkLine(:,[4,5]), repmat([0.0000  0.00000  0.00000 0.000    0.000    0.000  1],size(linkLine,1),1)]];
    m.extBusMap=[ones(size(busesInMain)),busesInMain,ones(size(busesInMain))*iSys,((nbusx+1):(nbusx+size(busesInMain,1)))']; % The definition of external bus map: (ext-zone, ext-bus, this-zone, this-bus)
    
    newSysDataList{iSys}=foldSysData(busx,sw,pv,pq,shunt,linex,ind,zip,syn,exc,tg,agc,cac,cluster);
    
    m.busMap=(size(busMain,1)+1):(size(busMain,1)+size(bus,1));
    m.synMap=(size(synMain,1)+1):(size(synMain,1)+size(syn,1));
    bus(:,1)=m.busMap;
    pv(:,1)=m.busMap(pv(:,1));
    pq(:,1)=m.busMap(pq(:,1));
    shunt(:,1)=m.busMap(shunt(:,1));
    line(:,1)=m.busMap(line(:,1));line(:,2)=m.busMap(line(:,2));
    ind(:,1)=m.busMap(ind(:,1));
    zip(:,1)=m.busMap(zip(:,1));
    syn(:,1)=m.busMap(syn(:,1));
    exc(:,1)=m.synMap(exc(:,1));
    tg(:,1)=m.synMap(tg(:,1));
    agc(:,1)=m.busMap(agc(:,1));
    % TODO: cac and cluster
    
    line=[line;[sysMaps{linkLine(1,1)}.busMap(linkLine(:,3))',    m.busMap(linkLine(:,2))',  repmat([100.00],size(linkLine,1),1),  bus(linkLine(:,2),2) , repmat([60, 0,   0.0000],size(linkLine,1),1),  linkLine(:,[4,5]), repmat([0.0000  0.00000  0.00000 0.000    0.000    0.000  1],size(linkLine,1),1)]];
    
    m.swMap=(size(swMain,1)+1):(size(swMain,1)+size(sw,1));
    m.pvMap=(size(pvMain,1)+1):(size(pvMain,1)+size(pv,1));
    m.pqMap=(size(pqMain,1)+1):(size(pqMain,1)+size(pq,1));
    m.lineMap=(size(lineMain,1)+1):(size(lineMain,1)+size(line,1));
    m.shuntMap=(size(shuntMain,1)+1):(size(shuntMain,1)+size(shunt,1));
    m.indMap=(size(indMain,1)+1):(size(indMain,1)+size(ind,1));
    m.zipMap=(size(zipMain,1)+1):(size(zipMain,1)+size(zip,1));
    m.synMap=(size(synMain,1)+1):(size(synMain,1)+size(syn,1));
    m.excMap=(size(excMain,1)+1):(size(excMain,1)+size(exc,1));
    m.tgMap=(size(tgMain,1)+1):(size(tgMain,1)+size(tg,1));
    m.agcMap=(size(agcMain,1)+1):(size(agcMain,1)+size(agc,1));
    m.cacMap=(size(cacMain,1)+1):(size(cacMain,1)+size(cac,1));
    m.clusterMap=(size(clusterMain,1)+1):(size(clusterMain,1)+size(cluster,1));
    
    busMain=[busMain;bus];
    pvMain=[pvMain;pv];
    pqMain=[pqMain;pq];
    shuntMain=[shuntMain;shunt];
    lineMain=[lineMain;line];
    indMain=[indMain;ind];
    zipMain=[zipMain;zip];
    synMain=[synMain;syn];    
    excMain=[excMain;exc];
    tgMain=[tgMain;tg];
    agcMain=[agcMain;agc];
    % TODO: cac and cluster
    
    grossMaps.busAreaMap=[grossMaps.busAreaMap;iSys*ones(size(bus,1),1)];
    grossMaps.swAreaMap=[grossMaps.swAreaMap;iSys*ones(size(sw,1),1)];
    grossMaps.pvAreaMap=[grossMaps.pvAreaMap;iSys*ones(size(pv,1),1)];
    grossMaps.pqAreaMap=[grossMaps.pqAreaMap;iSys*ones(size(pq,1),1)];
    grossMaps.shuntAreaMap=[grossMaps.shuntAreaMap;iSys*ones(size(shunt,1),1)];
    grossMaps.lineAreaMap=[grossMaps.lineAreaMap;iSys*ones(size(line,1),1)];
    grossMaps.indAreaMap=[grossMaps.indAreaMap;iSys*ones(size(ind,1),1)];
    grossMaps.zipAreaMap=[grossMaps.zipAreaMap;iSys*ones(size(zip,1),1)];
    grossMaps.synAreaMap=[grossMaps.synAreaMap;iSys*ones(size(syn,1),1)];
    grossMaps.excAreaMap=[grossMaps.excAreaMap;iSys*ones(size(exc,1),1)];
    grossMaps.tgAreaMap=[grossMaps.tgAreaMap;iSys*ones(size(tg,1),1)];
    grossMaps.agcAreaMap=[grossMaps.agcAreaMap;iSys*ones(size(agc,1),1)];
    grossMaps.cacAreaMap=[grossMaps.cacAreaMap;iSys*ones(size(cac,1),1)];
    grossMaps.clusterAreaMap=[grossMaps.clusterAreaMap;iSys*ones(size(cluster,1),1)];    
    
    sysMaps{iSys}=m;
end

grossMaps.invBusMap=zeros(size(grossMaps.busAreaMap));
grossMaps.invSwMap=zeros(size(grossMaps.swAreaMap));
grossMaps.invPvMap=zeros(size(grossMaps.pvAreaMap));
grossMaps.invPqMap=zeros(size(grossMaps.pqAreaMap));
grossMaps.invShuntMap=zeros(size(grossMaps.shuntAreaMap));
grossMaps.invLineMap=zeros(size(grossMaps.lineAreaMap));
grossMaps.invIndMap=zeros(size(grossMaps.indAreaMap));
grossMaps.invZipMap=zeros(size(grossMaps.zipAreaMap));
grossMaps.invSynMap=zeros(size(grossMaps.synAreaMap));
grossMaps.invExcMap=zeros(size(grossMaps.excAreaMap));
grossMaps.invTgMap=zeros(size(grossMaps.tgAreaMap));
grossMaps.invAgcMap=zeros(size(grossMaps.agcAreaMap));
grossMaps.invCacMap=zeros(size(grossMaps.cacAreaMap));
grossMaps.invClusterMap=zeros(size(grossMaps.clusterAreaMap));

for iSys=1:nSys
    grossMaps.invBusMap(sysMaps{iSys}.busMap)=1:size(sysMaps{iSys}.busMap,2);
    grossMaps.invSwMap(sysMaps{iSys}.swMap)=1:size(sysMaps{iSys}.swMap,2);
    grossMaps.invPvMap(sysMaps{iSys}.pvMap)=1:size(sysMaps{iSys}.pvMap,2);
    grossMaps.invPqMap(sysMaps{iSys}.pqMap)=1:size(sysMaps{iSys}.pqMap,2);
    grossMaps.invShuntMap(sysMaps{iSys}.shuntMap)=1:size(sysMaps{iSys}.shuntMap,2);
    grossMaps.invLineMap(sysMaps{iSys}.lineMap)=1:size(sysMaps{iSys}.lineMap,2);
    grossMaps.invIndMap(sysMaps{iSys}.indMap)=1:size(sysMaps{iSys}.indMap,2);
    grossMaps.invZipMap(sysMaps{iSys}.zipMap)=1:size(sysMaps{iSys}.zipMap,2);
    grossMaps.invSynMap(sysMaps{iSys}.synMap)=1:size(sysMaps{iSys}.synMap,2);
    grossMaps.invExcMap(sysMaps{iSys}.excMap)=1:size(sysMaps{iSys}.excMap,2);
    grossMaps.invTgMap(sysMaps{iSys}.tgMap)=1:size(sysMaps{iSys}.tgMap,2);
    grossMaps.invAgcMap(sysMaps{iSys}.agcMap)=1:size(sysMaps{iSys}.agcMap,2);
    grossMaps.invCacMap(sysMaps{iSys}.cacMap)=1:size(sysMaps{iSys}.cacMap,2);
    grossMaps.invClusterMap(sysMaps{iSys}.clusterMap)=1:size(sysMaps{iSys}.clusterMap,2);
end

sysDataGross=foldSysData(busMain,swMain,pvMain,pqMain,shuntMain,lineMain,indMain,zipMain,synMain,excMain,tgMain,agcMain,cacMain,clusterMain);

end