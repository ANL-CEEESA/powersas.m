function [fBus,fPV,fPQ,fSW,fLine,fLtc,fSupply,fDemand,fShunt,fInd,fZip,fSyn,fBusName,busNewToOld,busOldToNew,pvInd,pqInd,swInd,lineInd,ltcInd,supInd,demInd,shuntInd,indInd,zipInd,synInd]...
    =forgeLargestIsland(maxIsland,nIslands,islands,bus,pv,pq,sw,line,ltc,supply,demand,shunt,ind,zip,syn,busName,newToOld,oldToNew,...
    pvIndOld,pqIndOld,swIndOld,lineIndOld,ltcIndOld,supIndOld,demIndOld,shuntIndOld,indIndOld,zipIndOld,synIndOld)
% generate system data on the lagest island
%
% FUNCTION forgeLargestIsland
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
% assume that island are coded as 1:nIslands
% set the island with largest amount of load as main island and
% calculate power flow.

%     collectPQ=zeros(size(bus,1),size(pq,2));
%     collectPQ(pq(:,1),:)=pq;
%
% %     collectIsland=islands;
%
%     loads=collectPQ(:,4).*collectPQ(:,9);
%     loadsOnIslands=zeros(nIslands,1);
%     genAvailableIslands=zeros(nIslands,1);
%
%     genAvailableIslands(islands(pv(:,1)))=1;
%     genAvailableIslands(islands(sw(:,1)))=1;
%
%     for i=1:nIslands
%         loadsOnIslands(i)=sum(loads(islands==i));
%     end
%     loadsOnIslands=loadsOnIslands.*genAvailableIslands;
%     [maxLoad,maxIsland]=max(loadsOnIslands);

% forge data structure and re-number nodes.
busNewToOld=find(islands==maxIsland);
busOldToNew=zeros(max(bus(:,1)),1);
busOldToNew(bus(islands==maxIsland,1))=(1:length(busNewToOld))';

fBus=bus;
fBus(:,1)=busOldToNew(fBus(:,1));
fBus=fBus(fBus(:,1)~=0,:);

fPV=pv;
if ~isempty(pv)
    fPV(:,1)=busOldToNew(fPV(:,1));
    availPV=(fPV(:,1)~=0).*(fPV(:,11)~=0);
    pvInd=find(availPV~=0);
    fPV=fPV(availPV~=0,:);
else
    fPV=zeros(0,11);
    pvInd=[];
end

fPQ=pq;
if ~isempty(pq)
    fPQ(:,1)=busOldToNew(fPQ(:,1));
    pqInd=find(fPQ(:,1)~=0);
    fPQ=fPQ(fPQ(:,1)~=0,:);
else
    fPQ=zeros(0,9);
    pqInd=[];
end

fSW=sw;
fSW(:,1)=busOldToNew(fSW(:,1));
swInd=find(fSW(:,1)~=0);
fSW=fSW(fSW(:,1)~=0,:);

fLine=line;
fLine(:,1)=busOldToNew(fLine(:,1));
fLine(:,2)=busOldToNew(fLine(:,2));
lineInIsland=(fLine(:,1)~=0).*(fLine(:,2)~=0);
fLine=fLine(lineInIsland~=0,:);
lineInd=find(lineInIsland~=0);

if size(ltc,1)
    fLtc=ltc;
    fLtc(:,1)=busOldToNew(fLtc(:,1));
    fLtc(:,2)=busOldToNew(fLtc(:,2));
    ltcInIsland=(fLtc(:,1)~=0).*(fLtc(:,2)~=0);
    fLtc=fLtc(ltcInIsland~=0,:);
    ltcInd=find(ltcInIsland~=0);
else
    fLtc=zeros(0,16);
    ltcInd=[];
end

if ~isempty(supply)
    fSupply=supply;
    fSupply(:,1)=busOldToNew(fSupply(:,1));
    supInd=find(fSupply(:,1)~=0);
    fSupply=fSupply(fSupply(:,1)~=0,:);
else
    fSupply=zeros(0,20);
    supInd=[];
end

if size(demand,1)
    fDemand=demand;
    fDemand(:,1)=busOldToNew(fDemand(:,1));
    demInd=find(fDemand(:,1)~=0);
    fDemand=fDemand(fDemand(:,1)~=0,:);
else
    fDemand=zeros(0,18);
    demInd=[];
end

if size(shunt,1)
    fShunt=shunt;
    fShunt(:,1)=busOldToNew(fShunt(:,1));
    shuntInd=find(fShunt(:,1)~=0);
    fShunt=fShunt(fShunt(:,1)~=0,:);
else
    fShunt=zeros(0,7);
    shuntInd=[];
end

fInd=ind;
if ~isempty(ind)
    fInd(:,1)=busOldToNew(fInd(:,1));
    indInd=find(fInd(:,1)~=0);
    fInd=fInd(fInd(:,1)~=0,:);
else
    fInd=zeros(0,20);
    indInd=[];
end

fZip=zip;
if ~isempty(zip)
    fZip(:,1)=busOldToNew(fZip(:,1));
    zipInd=find(fZip(:,1)~=0);
    fZip=fZip(fZip(:,1)~=0,:);
else
    fZip=zeros(0,12);
    zipInd=[];
end

fSyn=syn;
if ~isempty(syn)
    fSyn(:,1)=busOldToNew(fSyn(:,1));
    synInd=find(fSyn(:,1)~=0);
    fSyn=fSyn(fSyn(:,1)~=0,:);
else
    fSyn=zeros(0,26);
    synInd=[];
end

if ~isempty(busName)
    fBusName=busName(busNewToOld);
else
    fBusName={};
end

%deal with swing bus problem

if size(fSW,1)==0&&size(fPV,1)>0
    pvline=fPV(1,:);
    fSW=[pvline(1:3),pvline(5),0,pvline(6:9),pvline(10),0,1,pvline(11)];
    pvNum=size(fPV,1);
    fPV=fPV(2:pvNum,:);
    swInd=-pvInd(1); % Negative Index.
    pvInd=pvInd(2:pvNum,:);
end

% Generate new indices.

if nargin>=30
    busNewToOld=newToOld(busNewToOld);
    busOldToNew=zeros(size(busNewToOld));
    busOldToNew(busNewToOld)=(1:size(busNewToOld,1))';
    
    pvIndTemp=pvIndOld;
    pvInd=pvIndOld(pvInd);
    pqInd=pqIndOld(pqInd);
    if ~isempty(swInd)
        if swInd(1)<0
            swInd=-pvIndTemp(-swInd(1)); % Negative Index
        else
            swInd=swIndOld(swInd);
        end
    else
        swInd=[];
    end
    lineInd=lineIndOld(lineInd);
    ltcInd=ltcIndOld(ltcInd);
    supInd=supIndOld(supInd);
    demInd=demIndOld(demInd);
    shuntInd=shuntIndOld(shuntInd);
    indInd=indIndOld(indInd);
    zipInd=zipIndOld(zipInd);
    synInd=synIndOld(synInd);
end

end