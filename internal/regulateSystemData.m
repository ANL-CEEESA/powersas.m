function [bus,pv,pq,sw,line,shunt,ind,zip,syn,exc,tg,agc,cac,cluster,pm,oldToNew,newToOld]=regulateSystemData(bus,pv,pq,sw,line,shunt,ind,zip,syn,exc,tg,agc,cac,cluster)
% Regulate input data, including checking integrity of input data, align base power, renumbering bus
%
% FUNCTION regulateSystemData
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
% INPUT (will be consolidated in a future version)
%
% OUTPUT (will be consolidated in a future version)
% 
	
if nargin<9;  syn=[];end
if nargin<10; exc=[];end
if nargin<11; tg=[]; end
if nargin<12; agc=[];end
if nargin<13; cac=[];end
if nargin<14; cluster=[];end

nbus=size(bus,1);

newToOld=bus(:,1);
oldToNew=zeros(max(newToOld),1);
oldToNew(newToOld)=1:nbus;
bus(:,1)=oldToNew(bus(:,1));
pv(:,1)=oldToNew(pv(:,1));
pq(:,1)=oldToNew(pq(:,1));
sw(:,1)=oldToNew(sw(:,1));
line(:,1:2)=oldToNew(line(:,1:2));
shunt(:,1)=oldToNew(shunt(:,1));
ind(:,1)=oldToNew(ind(:,1));
zip(:,1)=oldToNew(zip(:,1));
syn(:,1)=oldToNew(syn(:,1));

pv=pv(pv(:,end)==1,:);
pq=pq(pq(:,end)==1,:);
line=line(line(:,end)==1,:);
shunt=shunt(shunt(:,end)==1,:);
ind=ind(ind(:,end)==1,:);
zip=zip(zip(:,end)==1,:);
% syn=syn(syn(:,end)==1,:);
% exc=exc(exc(:,end)==1,:);
% tg=tg(tg(:,end)==1,:);

busShunt=zeros(nbus,2);
busPQ=zeros(nbus,2);
if isempty(sw)&&~isempty(pv)
	sw=zeros(1,13);
	sw(1,1:3)=pv(1,1:3);
	sw(1,4)=pv(1,5);
	sw(1,6:end)=[9.9  -9.9  1.2  0.8  2.324  1  1  1];
	pv=pv(2:end,:);
end

pvw=zeros(nbus,size(pv,2));
if ~isempty(pv)
    for idx=2:size(pv,2)
        pvw(:,idx)=accumarray(pv(:,1),pv(:,idx).*pv(:,end),[nbus,1]);
    end
    npv=accumarray(pv(:,1),ones(size(pv(:,1))),[nbus,1]);
    pvw(:,1)=1:nbus;
    pvw=pvw(npv>0,:);
    npv=npv(npv>0);
    pvw(:,[2,3,5,8,9,10,11])=pvw(:,[2,3,5,8,9,10,11])./repmat(npv,1,7);
    pvw=pvw(pvw(:,end)>0,:);
    pvw(:,end)=1;
    pv=pvw;
end

if ~isempty(pq)
    busPQ(:,1)=busPQ(:,1)+accumarray(pq(:,1),pq(:,4).*pq(:,9),[nbus,1]);
    busPQ(:,2)=busPQ(:,2)+accumarray(pq(:,1),pq(:,5).*pq(:,9),[nbus,1]);
end
if ~isempty(shunt)
    busShunt(:,1)=accumarray(shunt(:,1),shunt(:,5).*shunt(:,7),[nbus,1]);
    busShunt(:,2)=accumarray(shunt(:,1),shunt(:,6).*shunt(:,7),[nbus,1]);
end
if ~isempty(zip)
    busShunt(:,1)=busShunt(:,1)+accumarray(zip(:,1),zip(:,5).*zip(:,12),[nbus,1]);
    busShunt(:,2)=busShunt(:,2)+accumarray(zip(:,1),zip(:,8).*zip(:,12),[nbus,1]);
    busPQ(:,1)=busPQ(:,1)+accumarray(zip(:,1),zip(:,7).*zip(:,12),[nbus,1]);
    busPQ(:,2)=busPQ(:,2)+accumarray(zip(:,1),zip(:,10).*zip(:,12),[nbus,1]);
end

pq(:,4:5)=0;
shunt(:,5:6)=0;
zipBus=zeros(nbus,12);
zipBus(:,1)=1:nbus;
zipBus(:,2)=100;
zipBus(:,3)=bus(:,2);
zipBus(:,4)=60;
zipBus(:,[5,8])=busShunt;
zipBus(:,[7,10])=busPQ;
zipBus(:,12)=1;

if ~isempty(zip)
    zipBus(zip(:,1),[6,9])=zip(:,[6,9]);
end

zip=zipBus(sum(zipBus(:,5:10),2)~=0,:);

if ~isempty(syn)
    pg=zeros(nbus,1);
    if ~isempty(pv)
        pg(pv(:,1))=pv(:,4);
    end
    %         if ~isempty(sw)
    %             pg(sw(:,1))=sw(:,10);
    %         end
    pm=pg(syn(:,1));
    nm=accumarray(syn(:,1),ones(size(syn,1),1),[nbus,1]);
    pm=pm./nm(syn(:,1));
    pg=pg-accumarray(syn(:,1),pm,[nbus,1]);
    if ~isempty(pv)
        pv(:,4)=pg(pv(:,1));
    end
else
    pm=[];
end


% Regulate non-standard data
% Non-standard base power
sBase=100;
sw(:,6:7)=sw(:,6:7).*repmat(sw(:,2)/sBase,1,2);sw(:,2)=sBase;
pv(:,[4,6,7])=pv(:,[4,6,7]).*repmat(pv(:,2)/sBase,1,3);pv(:,2)=sBase;
pq(:,[4,5])=pq(:,[4,5]).*repmat(pq(:,2)/sBase,1,2);pq(:,2)=sBase;
shunt(:,[5,6])=shunt(:,[5,6]).*repmat(shunt(:,2)/sBase,1,2);shunt(:,2)=sBase;


end