function [SysData,maps,Ef,P0,pShare,Tm]=...
    generateBlackStartNet(SysDataBase,bsBus,bsSyn,bsInd)
% Generate system data structure and mappings from base system data
%
% FUNCTION generateBlackStartNet
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

    [busBase,swBase,pvBase,pqBase,shuntBase,lineBase,indBase,zipBase,synBase,excBase,tgBase,agcBase,cacBase,clusterBase]=unfoldSysData(SysDataBase);
    busName=[];
    nbus=size(busBase,1);
    
    bus=busBase(bsBus(:,1),:);
    busMap=zeros(nbus,1);
    busMap(bsBus(:,1))=1:size(bsBus,1);
    invBusMap=bsBus(:,1);
    bus(:,1)=busMap(bus(:,1));
    
    pv=pvBase(busMap(pvBase(:,1))~=0,:);pv(:,1)=busMap(pv(:,1));
    pq=pqBase(busMap(pqBase(:,1))~=0,:);pq(:,1)=busMap(pq(:,1));
    sw=swBase(busMap(swBase(:,1))~=0,:);sw(:,1)=busMap(sw(:,1));
    zip=zeros(size(bus,1),12);
    zip(:,1)=bus(:,1);zip(:,2)=100;zip(:,3)=bus(:,2);zip(:,4)=60;
    zip(:,5:10)=bsBus(:,2:7);
    zip(:,12)=1;
    
    shunt=zeros(0,7);
    line=zeros(0,16);
    ltc=[];
    sup=[];
    dem=[];
    
    synMap=zeros(size(synBase,1),1);
    invSynMap=bsSyn(:,1);
    synMap(bsSyn(:,1))=1:length(invSynMap);
    syn=synBase(invSynMap,:);
    syn(:,1)=busMap(synBase(bsSyn(:,1),1));
    exc=excBase;if ~isempty(excBase);exc(:,1)=synMap(excBase(:,1));exc=exc(exc(:,1)~=0,:);end
    tg=tgBase;if ~isempty(tgBase);tg(:,1)=synMap(tgBase(:,1));tg=tg(tg(:,1)~=0,:);end
    Ef=bsSyn(synMap(bsSyn(:,1))~=0,2);
    P0=bsSyn(synMap(bsSyn(:,1))~=0,3);
    pShare=bsSyn(synMap(bsSyn(:,1))~=0,4);
    
    indMap=zeros(size(indBase,1),1);
    if ~isempty(bsInd)
        invIndMap=bsInd(:,1);
        indMap(bsInd(:,1))=1:length(invIndMap);
        ind=indBase(invIndMap,:);
        ind(:,1)=busMap(indBase(bsInd(:,1),1));
        Tm=bsInd(indMap(bsInd(:,1))~=0,2);
        
        ind(:,15)=Tm;ind(:,16)=0;ind(:,17)=0;
    else
        invIndMap=[];
        indMap=[];
        ind=indBase;
        Tm=[];
    end
        
    lineMap=zeros(size(lineBase,1),1);
    invLineMap=zeros(0,1);
    
    agc=agcBase(busMap(agcBase(:,1))~=0,:);
    agc(:,1)=busMap(agc(:,1));
    cac=cacBase;
    cluster=clusterBase;
    
    SysData=foldSysData(bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster);
    
    maps.busMap=busMap;maps.invBusMap=invBusMap;
    maps.indMap=indMap;maps.invIndMap=invIndMap;
    maps.synMap=synMap;maps.invSynMap=invSynMap;
    maps.lineMap=lineMap;maps.invLineMap=invLineMap;
    
end