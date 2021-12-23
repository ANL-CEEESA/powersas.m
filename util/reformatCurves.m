function [t,stateCurve]=reformatCurves(SysDataBase,tList,stateList,SysDataList,mapsList,dt)
% Reformat the curves from raw computation results
%
% FUNCTION reformatCurves
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
%
% OUTPUT
%

if nargin<6
    dt=[];
end
[nState,idxs]...
    =getIndexDyn(SysDataBase);
[busBase,swBase,pvBase,pqBase,shuntBase,lineBase,indBase,zipBase,synBase,excBase,tgBase]=unfoldSysData(SysDataBase);

if isempty(dt)
    tListx=tList;
    tListx1=tList;
else
    tListx=cell(1,length(tList));
    tListx1=cell(1,length(tList));
    residue=0.0;
    for it=1:length(tList)
        if length(tList{it})<=1
            seg=tList{it};
            tListx{it}=seg;
            tListx1{it}=seg;
            residue=residue+seg(end)-floor((residue+seg(end))/dt)*dt;
        else
            seg=((ceil((tList{it}(1)+residue)/dt)*dt):dt:(floor((tList{it}(end)+residue)/dt)*dt))-residue;
            seg=seg(seg>=0);
            if ~isempty(seg)
                if tList{it}(end)>seg(end)
                    seg=[seg,tList{it}(end)];
                end
                if seg(1)>0
                    seg=[0,seg];
                end
            else
                seg=[0,tList{it}(end)];
            end
            tListx{it}=seg;
            tListx1{it}=seg;          
            residue=residue+seg(end)-floor((residue+seg(end))/dt)*dt;
        end
    end
end
t=horzcatTime(tListx1{:});
stateCurve=nan(nState,length(t));
startIdx=1;

synToTg=zeros(size(synBase,1),1);
synToTg(tgBase(:,1))=1:size(tgBase,1);
synToExc=zeros(size(synBase,1),1);
synToExc(excBase(:,1))=1:size(excBase,1);

for seq=1:length(tListx)
    if isempty(dt)
        state=stateList{seq};
    else
        if length(tList{seq})<=1
            state=stateList{seq};
        else
            state=interp1(tList{seq},stateList{seq}',tListx{seq})';
        end
    end
    SysData=SysDataList{seq};
    busMap=mapsList{seq}.busMap;
    invBusMap=mapsList{seq}.invBusMap;
    indMap=mapsList{seq}.indMap;
    invIndMap=mapsList{seq}.invIndMap;
    synMap=mapsList{seq}.synMap;
    invSynMap=mapsList{seq}.invSynMap;
    
    [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg]=unfoldSysData(SysData);
    [V,Q,s,d,w,eq1,eq2,ed1,ed2,psid,psiq,Pm,Ef,Vavrm,Vavrr,Vavrf,Vavrref,tgovg,tgovm,tgovmech,f,dpg,qplt,vg]=...
        unfoldX(state,SysData);
    
    iTseq=startIdx:(startIdx+length(tListx{seq})-1);
    stateCurve(idxs.vIdx(invBusMap),iTseq)=V;
    stateCurve(idxs.qIdx(invBusMap),iTseq)=Q;
    stateCurve(idxs.sIdx(invIndMap),iTseq)=s;
    stateCurve(idxs.deltaIdx(invSynMap),iTseq)=d;
    stateCurve(idxs.omegaIdx(invSynMap),iTseq)=w;
    stateCurve(idxs.eq1Idx(invSynMap),iTseq)=eq1;
    stateCurve(idxs.eq2Idx(invSynMap),iTseq)=eq2;
    stateCurve(idxs.ed1Idx(invSynMap),iTseq)=ed1;
    stateCurve(idxs.ed2Idx(invSynMap),iTseq)=ed2;
    stateCurve(idxs.psidIdx(invSynMap),iTseq)=psid;
    stateCurve(idxs.psiqIdx(invSynMap),iTseq)=psiq;
    stateCurve(idxs.pgIdx(invSynMap),iTseq)=Pm;
    stateCurve(idxs.efIdx(invSynMap),iTseq)=Ef;
    stateCurve(idxs.vavrmIdx(synToExc(invSynMap(exc(:,1)))),iTseq)=Vavrm;
    stateCurve(idxs.vavrrIdx(synToExc(invSynMap(exc(:,1)))),iTseq)=Vavrr;
    stateCurve(idxs.vavrfIdx(synToExc(invSynMap(exc(:,1)))),iTseq)=Vavrf;
    stateCurve(idxs.vavrrefIdx(synToExc(invSynMap(exc(:,1)))),iTseq)=Vavrref;
    stateCurve(idxs.tgovgIdx(synToTg(invSynMap(tg(:,1)))),iTseq)=tgovg;
    stateCurve(idxs.tgovmIdx(synToTg(invSynMap(tg(:,1)))),iTseq)=tgovm;
    stateCurve(idxs.tmechIdx(synToTg(invSynMap(tg(:,1)))),iTseq)=tgovmech;  
    stateCurve(idxs.fIdx(invBusMap),iTseq)=f;
    stateCurve(idxs.dpgIdx(invBusMap),iTseq)=dpg;  
    
    startIdx=startIdx+length(tListx{seq});
end
end