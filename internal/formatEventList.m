function [eventExt,evtDyn,evtDynInd,evtDynZip]=formatEventList(event,evtDyn,evtDynInd,evtDynZip)
% Regulate non-standard event list to standard format
%
% FUNCTION formatEventList
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
EvtType=getSimEventType();
MethodType=getSimMethodType();
SimSet=getSimSettingDefault();

if size(event,2)<6
    event(:,6)=MethodType.FULL_HE;
end
if size(event,2)<7
    event(:,7)=0;
    event(event(:,4)==EvtType.DYN_SIM,7)=SimSet.DEFAULT_DT;
end

[eSortTemp,~]=sort([event(:,2);event(:,3)]);
diffEvt=eSortTemp(2:end)-eSortTemp(1:end-1);
diffEvt=diffEvt(diffEvt>1e-10);
minDiffTime=min([diffEvt;SimSet.SEG/10]);

[~,iSort]=sort(event(:,2)+(event(:,4)==EvtType.DYN_SIM)*minDiffTime/2+(event(:,4)==EvtType.DYN_SIM).*event(:,3)/eSortTemp(end)*minDiffTime/4);
eventSort=event(iSort,:);
endTagIdx=find(eventSort(:,4)==EvtType.END,1);
startT=0.;
if ~isempty(endTagIdx)
    endT=eventSort(endTagIdx,2);
else
    endT=max(max(event(:,2:3)));
end

eventSort=eventSort(eventSort(:,2)<=endT,:);
eventSort(eventSort(:,3)>endT,3)=endT;
eventSort(eventSort(:,4)~=EvtType.DYN_SIM,3)=eventSort(eventSort(:,4)~=EvtType.DYN_SIM,2);
diffTime=eventSort(:,3)-eventSort(:,2);
eventSort=eventSort(eventSort(:,4)~=EvtType.DYN_SIM|diffTime>0,:);
eventSort(1:end-1,3)=min([eventSort(1:end-1,3),eventSort(2:end,2)],[],2);

diffTime=eventSort(:,3)-eventSort(:,2);
eventSort=eventSort(eventSort(:,4)~=EvtType.DYN_SIM|diffTime>0,:);
timespots=reshape([eventSort(:,2)';eventSort(:,3)'],2*size(eventSort,1),1);
typeTag=reshape([eventSort(:,4)';zeros(size(eventSort,1),1)'],2*size(eventSort,1),1)==EvtType.DYN_SIM;
intervalDiff=[timespots(2:end)-timespots(1:end-1);0];

timespotsNz=[timespots(intervalDiff~=0);endT];
typeTagNz=typeTag(intervalDiff~=0);
stTime=timespotsNz(1:end-1);stTime=stTime(typeTagNz==0);
endTime=timespotsNz(2:end);endTime=endTime(typeTagNz==0);

[sstTime,iSort]=sort(stTime);
sendTime=endTime(iSort);
stTime=zeros(0,1);endTime=zeros(0,1);
if ~isempty(sstTime)
    for i=1:size(sstTime,1)
        if sendTime(i)-sstTime(i)>SimSet.SEG
            interval=sstTime(i):SimSet.SEG:sendTime(i);
            if interval(end)<sendTime(i)-1e-9
                interval=[interval,sendTime(i)];
            end
            stTime=[stTime;interval(1:end-1)'];
            endTime=[endTime;interval(2:end)'];
        else
            stTime=[stTime;sstTime(i)];
            endTime=[endTime;sendTime(i)];
        end
    end
end
eventAdd=[zeros(size(stTime,1),1),stTime,endTime,EvtType.DYN_SIM*ones(size(stTime,1),1),zeros(size(stTime,1),1),MethodType.FULL_HE*ones(size(stTime,1),1),SimSet.DEFAULT_DT*ones(size(stTime,1),1)];

eventExt=[eventSort;eventAdd];
[~,iSort]=sort(eventExt(:,2)+(eventExt(:,4)==EvtType.DYN_SIM)*minDiffTime/2+(eventExt(:,4)==EvtType.DYN_SIM).*eventExt(:,3)/eSortTemp(end)*minDiffTime/4);
eventExt=eventExt(iSort,:);

eventExt(:,1)=1:size(eventExt,1);
if eventExt(end,4)~=EvtType.END
    eventExt=[eventExt;[eventExt(end,1)+1,max(eventExt(end,2:3)),0,EvtType.END,0,MethodType.FULL_HE,0.0]];
end

idxDynEvt=find(eventExt(:,4)==EvtType.DYN_SIM&eventExt(:,5)~=0);
for i=1:size(idxDynEvt,1)
    item=idxDynEvt(i);
    event=eventExt(item,:);
    sstTime=event(2);
    sendTime=event(3);
    interval=sstTime:SimSet.SEG:sendTime;
    if interval(end)<sendTime-1e-9
        interval=[interval,sendTime];
    end
    stTime=interval(1:end-1)';
    endTime=interval(2:end)';
    numNewInterval=size(stTime,1);
    eventExt(item,3)=endTime(1);
    eventAdd=repmat(eventExt(item,:),numNewInterval-1,1);
    eventIdxNew=linspace(eventExt(item,1),eventExt(item,1)+1,numNewInterval+1);
    eventAdd(:,1)=eventIdxNew(2:(end-1));
    eventAdd(:,2)=stTime(2:end);
    eventAdd(:,3)=endTime(2:end);
    for iAdd=1:(numNewInterval-1)
        evtDynAdd=evtDyn(eventExt(item+iAdd-1,5),:);
        if any(evtDynAdd(:,[6,7,8,9])>0)
            if all(evtDynAdd(:,[6,7])>0) % Has induction motor item
                dynIndAdd=evtDynInd(evtDynAdd(6):evtDynAdd(7),:);
                dynIndAdd(:,2)=dynIndAdd(:,2)./(1+(stTime(iAdd+1)-sstTime)*dynIndAdd(:,2));
                evtDynAdd(6)=size(evtDynInd,1)+1;
                evtDynAdd(7)=size(evtDynInd,1)+size(dynIndAdd,1);
                evtDynInd=[evtDynInd;dynIndAdd];
            end
            if all(evtDynAdd(:,[8,9])>0)
                dynZipAdd=evtDynZip(evtDynAdd(8):evtDynAdd(9),:);
                dynZipAdd(:,2)=dynZipAdd(:,2)./(1+(stTime(iAdd+1)-sstTime)*dynZipAdd(:,2));
                evtDynAdd(8)=size(evtDynZip,1)+1;
                evtDynAdd(9)=size(evtDynZip,1)+size(dynZipAdd,1);
                evtDynZip=[evtDynZip;dynZipAdd];
            end
            evtDynAdd(1)=size(evtDyn,1)+1;
            eventAdd(iAdd,5)=evtDynAdd(1);
            evtDyn=[evtDyn;evtDynAdd];
        end
        eventExt=[eventExt(1:(item+iAdd-1),:);eventAdd(iAdd,:);eventExt((item+iAdd):end,:)];
    end
    idxDynEvt((i+1):end)=idxDynEvt((i+1):end)+numNewInterval-1;
end

eventExt(:,1)=1:size(eventExt,1);
end