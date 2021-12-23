function res=runDynamicSimulationExec(caseName,SysDataBase,simSettings,options,snapshot)
% Main function for running dynamic simulation
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
% FUNCTION runDynamicSimulationExec
%
% INPUT
%   caseName - Name of the studied case
%   SysDataBase - Studied system data struct
%   simSettings - Settings
%   options - Other options
%   snapshot - snap shot containing initial state
%
% OUTPUT
%   res struct:
%       flag - 
%       msg - Message returned
%       t - time sequence
%       stateCurve - list of system states corresponding to t
%       timestamp - the timestamp of the analysis
%       SysDataBase - Studied system data struct
%       eventList - list of system events
%       snapshot - the snapshot at the end of the analysis
%       caseName - Name of the studied case
%       simSettings - Settings
%

flag=0;
res=[];

if nargin<4
    options=regulateOptions();
else
    options=regulateOptions(options);
end
if options.dbstop
    dbstop if error
end

if isfield(options,'timestamp')
    timestamp=options.timestamp;
else
    timestamp=datestr(clock,'yyyymmddTHHMMSSFFF');
end

addLog('================= START ==================','INFO');
% addLog('--------- WARNING ---------','WARN');
% addLog(disclaimerStr,'WARN');

% restSettingFile='restoration_plan_test_039';
comment=options.comment;

if isfield(simSettings,'nlvl');nlvl=simSettings.nlvl;else nlvl=options.nlvl;end
if isfield(simSettings,'segAlpha');segAlpha=simSettings.segAlpha;else segAlpha=options.segAlpha;end
if isfield(simSettings,'dAlpha');dAlpha=simSettings.dAlpha;else dAlpha=options.dAlpha;end
if isfield(simSettings,'taylorN');taylorN=simSettings.taylorN;else taylorN=options.taylorN;end
if isfield(simSettings,'alphaTol');alphaTol=simSettings.alphaTol;else alphaTol=options.alphaTol;end
if isfield(simSettings,'diffTol');diffTol=simSettings.diffTol;else diffTol=options.diffTol;end
if isfield(simSettings,'diffTolMax');diffTolMax=simSettings.diffTolMax;else diffTolMax=options.diffTolMax;end
if isfield(simSettings,'method');method=simSettings.method;else method=options.method;end
if isfield(simSettings,'diffTolCtrl');diffTolCtrl=simSettings.diffTolCtrl;else diffTolCtrl=options.diffTolCtrl;end

eventList=simSettings.eventList;
bsSyn=simSettings.bsSyn;
bsBus=simSettings.bsBus;
bsInd=simSettings.bsInd;
Efstd=simSettings.Efstd;
evtLine=simSettings.evtLine;
evtLineSpec=simSettings.evtLineSpec;
evtZip=simSettings.evtZip;
evtZipSpec=simSettings.evtZipSpec;
evtZipSpec2=simSettings.evtZipSpec2;
evtInd=simSettings.evtInd;
evtIndSpec=simSettings.evtIndSpec;
evtSyn=simSettings.evtSyn;
evtSynSpec=simSettings.evtSynSpec;
evtFault=simSettings.evtFault;
evtFaultSpec=simSettings.evtFaultSpec;
evtDyn=simSettings.evtDyn;
evtDynPQ=simSettings.evtDynPQ;
evtDynPV=simSettings.evtDynPV;
evtDynInd=simSettings.evtDynInd;
evtDynZip=simSettings.evtDynZip;
evtDynSh=simSettings.evtDynSh;
evtDynZipRamp=simSettings.evtDynZipRamp;
evtDynTmech=simSettings.evtDynTmech;
evtDynPm=simSettings.evtDynPm;
evtDynEf=simSettings.evtDynEf;
evtDynVref=simSettings.evtDynVref;
evtDynEq1=simSettings.evtDynEq1;

if exist('envDependentScript','file')    
    envDependentScript;
end
[eventList,evtDyn,evtDynInd,evtDynZip]=formatEventList(eventList,evtDyn,evtDynInd,evtDynZip);
EvtType=getSimEventType();
MethodType=getSimMethodType();
SimSet=getSimSettingDefault();
DynSimFlag=getDynSimFlags();

if ~options.inheritParaSettings
    nlvl = options.nlvl;
    taylorN =  options.taylorN;
    alphaTol =options.alphaTol;
    diffTol =options.diffTol;
    diffTolMax =options.diffTolMax;
    method = options.method;
    diffTolCtrl =options.diffTolCtrl;
end

tTol=options.tTol;
diffTolOrig=diffTol;
allowReverse=options.allowReverse;
allowSteadyDynSwitch=options.allowSteadyDynSwitch;
useDiffCtrl=options.useDiffCtrl;
useNoMoveCtrl=options.useNoMoveCtrl;
nPool=options.nPool;

switchMaskCondition='event(2)>10 && ~(eventList(item+1,4)==EvtType.DYN_SIM&&eventList(item+1,5)~=0&&eventList(item+1,3)-eventList(item+1,2)<10)';
if useDiffCtrl&&(~exist('diffTolCtrl','var')||isempty(diffTolCtrl));useDiffCtrl=0;end

usePar=license('test','Distrib_Computing_Toolbox');
if usePar;delete(gcp('nocreate'));end
if usePar&&exist('nPool','var')&&nPool>1
    parpool(nPool);
    pctRunOnAll warning('off','MATLAB:singularMatrix');
    pctRunOnAll warning('off','MATLAB:nearlySingularMatrix');
    pctRunOnAll warning('off','MATLAB:illConditionedMatrix');
else
    warning('off','MATLAB:singularMatrix');
    warning('off','MATLAB:nearlySingularMatrix');
    warning('off','MATLAB:illConditionedMatrix');
end

hotStart=options.hotStart;
if hotStart && nargin>=5 && isstruct(snapshot) && isfield(snapshot,'SysData') && isfield(snapshot,'state')
    
else    
    hotStart=0;
end

%% Temporary code here
% eventList(:,6)=0.0;
% eventList(:,7)=0.1;
% Efstd=1.2;
% evtDynZip(:,2)=0.2;
%%
startTag=tic;
initSuccess=1;
if ~hotStart
    [busBase,swBase,pvBase,pqBase,shuntBase,lineBase,indBase,zipBase,synBase,excBase,tgBase,agcBase,cacBase,clusterBase]=unfoldSysData(SysDataBase);
    
    if eventList(1,5)==0 % Follow black start variables
        % TODO: AGC, CAC and Cluster
        nbusBase=size(busBase,1);
        synTag=zeros(nbusBase,1);
        if ~isempty(synBase);synTag(synBase(:,1))=1;end
        swBase=swBase(synTag(swBase(:,1))==0,:);
        if ~isempty(pvBase);pvBase=pvBase(synTag(pvBase(:,1))==0,:);end
        
        [busBase,pvBase,pqBase,swBase,lineBase,shuntBase,indBase,zipBase,synBase,excBase,tgBase,agcBase,cacBase,clusterBase,pm]=...
            regulateSystemData(busBase,pvBase,pqBase,swBase,lineBase,shuntBase,indBase,zipBase,synBase,excBase,tgBase,agcBase,cacBase,clusterBase);
        SysDataBase=foldSysData(busBase,swBase,pvBase,pqBase,shuntBase,lineBase,indBase,zipBase,synBase,excBase,tgBase,agcBase,cacBase,clusterBase);
        [SysData,maps,Efstd,P0,pShare,Tm]=...
            generateBlackStartNet(SysDataBase,bsBus,bsSyn,bsInd);
        
        [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
        SimDataStatic=foldSimData([],segAlpha,dAlpha,nlvl,taylorN,alphaTol,diffTol,diffTolMax,[]);
        SysPara=foldSysPara([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],(Efstd-1).*ones(size(syn,1),1),[],[]);
        
        [SysData,x0,finalAlpha,alphaList,diff]=calculateInitialState(SysData,SimDataStatic,SysPara);
    else % Full system steady-state
        
        [busBase,pvBase,pqBase,swBase,lineBase,shuntBase,indBase,zipBase,synBase,excBase,tgBase,agcBase,cacBase,clusterBase,pm]=...
            regulateSystemData(busBase,pvBase,pqBase,swBase,lineBase,shuntBase,indBase,zipBase,synBase,excBase,tgBase,agcBase,cacBase,clusterBase);
        bus=busBase;sw=swBase;pv=pvBase;pq=pqBase;shunt=shuntBase;
        line=lineBase;ind=indBase;zip=zipBase;syn=synBase;exc=excBase;tg=tgBase;
        agc=agcBase;cac=cacBase;cluster=clusterBase;
        maps.busMap=(1:size(busBase,1))';maps.invBusMap=maps.busMap;
        maps.indMap=(1:size(indBase,1))';maps.invIndMap=maps.indMap;
        maps.synMap=(1:size(synBase,1))';maps.invSynMap=maps.synMap;
        maps.lineMap=(1:size(lineBase,1))';maps.invLineMap=maps.lineMap;
        SimDataStatic=foldSimData([],segAlpha,[],nlvl,taylorN,alphaTol,diffTol,diffTolMax,[]);
        SysData=foldSysData(bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster);
        SysDataBase=foldSysData(busBase,swBase,pvBase,pqBase,shuntBase,lineBase,indBase,zipBase,synBase,excBase,tgBase,agcBase,cacBase,clusterBase);
        if size(Efstd,1)==1;Efstd=Efstd*ones(size(syn,1),1);end
        if isempty(Efstd);Ef1=[];else Ef1=(Efstd-1).*ones(size(syn,1),1);end
        SysPara=foldSysPara([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],Ef1,[],[]);
        
        [SysData,x0,finalAlpha,alphaList,diff]=calculateInitialState(SysData,SimDataStatic,SysPara);
    end
    
    if finalAlpha<1
        initSuccess=0;
        msg='Fail to get initial steady-state point. Abort.';
        flag=-1;
    end
    fault=zeros(0,4);
    
    snapshot.SysData=SysData;
    snapshot.state=x0;
    snapshot.t=0;
    snapshot.maps=maps;
    snapshot.fault=fault;
    snapshot.diff=diff;
    snapshot.step=0;
    SnapshotLilst={snapshot};
    
    SysDataList={SysData};
    stateList={x0};
    tList={0};
    mapsList={maps};
    faultList={fault};
    diffList={diff};
    stepsList={0};
    currState=x0;
    startIdx=1;
else
    [busBase,swBase,pvBase,pqBase,shuntBase,lineBase,indBase,zipBase,synBase,excBase,tgBase,agcBase,cacBase,clusterBase]=unfoldSysData(SysDataBase);
    if ~isfield(snapshot,'maps')
        maps.busMap=(1:size(busBase,1))';maps.invBusMap=maps.busMap;
        maps.indMap=(1:size(indBase,1))';maps.invIndMap=maps.indMap;
        maps.synMap=(1:size(synBase,1))';maps.invSynMap=maps.synMap;
        maps.lineMap=(1:size(lineBase,1))';maps.invLineMap=maps.lineMap;
        snapshot.maps=maps;
    end
    if ~isfield(snapshot,'fault')
        fault=zeros(0,4);
        snapshot.fault=fault;
    end
    if ~isfield(snapshot,'diff')
        diff=0;
        snapshot.diff=diff;
    end
    if ~isfield(snapshot,'t')
        snapshot.t=0;
    end
    if ~isfield(snapshot,'step')
        snapshot.step=0;
    end
    
    currState=snapshot.state(:,end);
    SysData=snapshot.SysData;
    t=0;
    maps=snapshot.maps;
    fault=snapshot.fault;
    diff=snapshot.diff;
    startIdx=1;
    SnapshotLilst={snapshot};
    
    SysDataList={SysData};
    stateList={currState};
    tList={0};
    mapsList={maps};
    faultList={fault};
    diffList={diff};
    stepsList={0};

end
item=startIdx;
SysDataOrig=[];
lastSwitchT=0.0;

switchDelay=options.switchDelay;
notifyCtrlDiff=options.notifyCtrlDiff;

while initSuccess&&item<=size(eventList,1)
    event=eventList(item,:);
    eventType=event(4);
    
    if eventType==EvtType.START %start / black start
        
    elseif eventType==EvtType.ADD_LINE %add line
        pos=event(5);
        method=event(6);
        evtLineSt=evtLine(pos,2);evtLineEnd=evtLine(pos,3);
        eventSpec=evtLineSpec(evtLineSt:evtLineEnd,:);
        eventSpec=eventSpec(eventSpec(:,2)==0,:);
        lineAdd=lineBase(eventSpec(:,1),:);
        lineBusNew=[maps.busMap(lineAdd(:,1)),maps.busMap(lineAdd(:,2))];
        lineAddConnIdx=eventSpec(lineBusNew(:,1)~=0&lineBusNew(:,2)~=0,1);
        lineAddNewIdx=eventSpec((lineBusNew(:,1)==0&lineBusNew(:,2)~=0)|(lineBusNew(:,1)~=0&lineBusNew(:,2)==0),1);
        
        if diff(end)>diffTol;diffTolSim=diff(end)+1e-9;else diffTolSim=diffTol;end
        SimData=foldSimData([],1,[],nlvl,taylorN,alphaTol,diffTolSim,diffTolMax,method);
        SysPara=foldSysPara([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],fault);
        [stateNew,finalAlpha,alphaList,diff,SysData,maps]=...
            simulationAddConnectionLine(SysData,SysDataBase,SysPara,SimData,maps,currState,lineAddConnIdx);
        if finalAlpha<1-tTol
            msg=['[',num2str(event(2)),']',' Add connection line causes singularity. Exit'];
            flag=-1;
            addLog(['Add connection line causes singularity. Exit'],'INFO');
            break;
        end
        
        snapshot.SysData=SysData;
        snapshot.state=stateNew;
        snapshot.t=0;
        snapshot.maps=maps;
        snapshot.fault=fault;
        snapshot.diff=diff;
        snapshot.step=0;
        SnapshotLilst{end+1}=snapshot;
        
        SysDataList{end+1}=SysData;
        stateList{end+1}=stateNew;tList{end+1}=0;
        mapsList{end+1}=maps;
        faultList{end+1}=fault;
        diffList{end+1}=diff;stepsList{end+1}=0;       
        currState=stateNew;
        
        SysPara=foldSysPara([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],fault);
        [stateUpd,finalAlpha,alphaList,diff,SysData,maps]=...
            simulationAddNewLine(SysData,SysDataBase,SysPara,SimData,maps,stateNew,lineAddNewIdx);
        if finalAlpha<1-tTol
            msg=['[',num2str(event(2)),']',' Add new line causes singularity. Exit'];
            flag=-1;
            addLog(['Add new line causes singularity. Exit'],'INFO');
            break;
        end
        
        snapshot.SysData=SysData;
        snapshot.state=stateUpd;
        snapshot.t=0;
        snapshot.maps=maps;
        snapshot.fault=fault;
        snapshot.diff=diff;
        snapshot.step=0;
        SnapshotLilst{end+1}=snapshot;
        
        SysDataList{end+1}=SysData;
        stateList{end+1}=stateUpd;tList{end+1}=0;
        mapsList{end+1}=maps;
        faultList{end+1}=fault;
        diffList{end+1}=diff;stepsList{end+1}=0;        
        currState=stateUpd;
        addLog(['[E',num2str(event(1)),'] Added Line(s) ',num2str(eventSpec(:,1)')],'INFO');
    elseif eventType==EvtType.ADD_LOAD %add static load
        pos=event(5);
        method=event(6);
        evtZipSt=evtZip(pos,3);evtZipEnd=evtZip(pos,4);
        if evtZip(pos,2)==0
            eventSpec=evtZipSpec(evtZipSt:evtZipEnd,:);
            eventSpec=eventSpec(eventSpec(:,2)==0,:);
            zipAdd=zipBase(eventSpec(:,1),:);
        else
            zipAdd=evtZipSpec2(evtZipSt:evtZipEnd,:);
            eventSpec=zipAdd;
        end
        
        if diff(end)>diffTol;diffTolSim=diff(end)+1e-9;else diffTolSim=diffTol;end
        [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
        nbus=size(bus,1);
        busTag=ones(nbus,1);
        
        SimData=foldSimData([],1,[],nlvl,taylorN,alphaTol,diffTolSim,diffTolMax,method);
        doContinue=1;
        SysPara=foldSysPara([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],fault);
        while 1
            [stateNew,finalAlpha,alphaList,diff,SysDataTemp,SysPara,mapsTemp]=...
                simulationAddZipLoadInstant(busTag,SysData,SysDataBase,SysPara,SimData,maps,currState,zipAdd);
            if finalAlpha<1-tTol
                [flag,doContinue,willBreak,msg]=processCollaps(SysDataTemp,stateNew(:,end),event,'Add static load',flag,doContinue);
                if willBreak
                    break; 
                end
            else
                SysData=SysDataTemp;
                maps=mapsTemp;

                break;
            end
        end
        if ~doContinue
            break;
        end
        
        snapshot.SysData=SysData;
        snapshot.state=stateNew;
        snapshot.t=0;
        snapshot.maps=maps;
        snapshot.fault=fault;
        snapshot.diff=diff;
        snapshot.step=0;
        SnapshotLilst{end+1}=snapshot;
        
        SysDataList{end+1}=SysData;
        stateList{end+1}=stateNew;tList{end+1}=0;
        mapsList{end+1}=maps;
        faultList{end+1}=fault;
        diffList{end+1}=diff;stepsList{end+1}=0;
        
        currState=stateNew;
        addLog(['[E',num2str(event(1)),'] Added Zip(s) ',num2str(eventSpec(:,1)')],'INFO');
    elseif eventType==EvtType.ADD_IND %add dynamic load
        pos=event(5);
        method=event(6);
        evtIndSt=evtInd(pos,2);evtIndEnd=evtInd(pos,3);
        eventSpec=evtIndSpec(evtIndSt:evtIndEnd,:);
        eventSpecMod=eventSpec(eventSpec(:,2)==1,:);
        eventSpec=eventSpec(eventSpec(:,2)==0,:);
        TmAdd=eventSpec(:,3);
        s0Add=eventSpec(:,4);
        
        if diff(end)>diffTol;diffTolSim=diff(end)+1e-9;else diffTolSim=diffTol;end
        [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
        nbus=size(bus,1);
        busTag=ones(nbus,1);
        
        SimData=foldSimData([],1,[],nlvl,taylorN,alphaTol,diffTolSim,diffTolMax,method);
        doContinue=1;
        SysPara=foldSysPara([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],fault);               
        while 1
            [stateNew,finalAlpha,alphaList,diff,SysDataTemp,SysPara,mapsTemp]=...
                simulationAddMotorInstant(busTag,SysData,SysDataBase,SysPara,SimData,maps,currState,eventSpec(:,1),TmAdd,s0Add);            
            if finalAlpha<1-tTol
                [flag,doContinue,willBreak,msg]=processCollaps(SysDataTemp,stateNew(:,end),event,'Add induction motor',flag,doContinue);
                if willBreak
                    break; 
                end
            else
                SysData=SysDataTemp;
                maps=mapsTemp;
                break;
            end
        end
        if ~doContinue
            break;
        end
        
        if ~isempty(eventSpecMod)
            [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
            indModIdx=maps.indMap(eventSpecMod(:,1));
            eventSpecMod=eventSpecMod(indModIdx~=0,:);
            indModIdx=indModIdx(indModIdx~=0);
            
            ind(indModIdx,15)=eventSpecMod(:,3);
            ind(indModIdx,16:17)=0;
            SysData=foldSysData(bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster);
        end
        
        snapshot.SysData=SysData;
        snapshot.state=stateNew;
        snapshot.t=0;
        snapshot.maps=maps;
        snapshot.fault=fault;
        snapshot.diff=diff;
        snapshot.step=0;
        SnapshotLilst{end+1}=snapshot;
        
        SysDataList{end+1}=SysData;
        stateList{end+1}=stateNew;tList{end+1}=0;
        mapsList{end+1}=maps;
        faultList{end+1}=fault;
        diffList{end+1}=diff;stepsList{end+1}=0;
        
        currState=stateNew;
        addLog(['[E',num2str(event(1)),'] Added Ind(s) ',num2str(eventSpec(:,1)')],'INFO');
    elseif eventType==EvtType.ADD_SYN %add syn gen
        pos=event(5);
        method=event(6);
        evtSynSt=evtSyn(pos,2);evtSynEnd=evtSyn(pos,3);
        eventSpec=evtSynSpec(evtSynSt:evtSynEnd,:);
        eventSpec=eventSpec(eventSpec(:,2)==0,:);
        d0Add=eventSpec(:,3);
        Pm0=eventSpec(:,4);
        Ef0=eventSpec(:,5);
        
        if diff(end)>diffTol;diffTolSim=diff(end)+1e-9;else diffTolSim=diffTol;end
        [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
        nbus=size(bus,1);
        busTag=ones(nbus,1);
        
        SimData=foldSimData([],1,[],nlvl,taylorN,alphaTol,diffTolSim,diffTolMax,method);
        doContinue=1;
        SysPara=foldSysPara([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],fault);                
        while 1
            [stateNew,finalAlpha,alphaList,diff,SysDataTemp,SysPara,mapsTemp]=...
                simulationAddSynInstant(busTag,SysData,SysDataBase,SysPara,SimData,maps,...
                currState,eventSpec(:,1),d0Add,Pm0,Ef0);            
            if finalAlpha<1-tTol
                [flag,doContinue,willBreak,msg]=processCollaps(SysDataTemp,stateNew(:,end),event,'Add syn gen',flag,doContinue);
                if willBreak
                    break; 
                end
            else
                SysData=SysDataTemp;
                maps=mapsTemp;
                break;
            end
        end
        if ~doContinue
            break;
        end
        
        snapshot.SysData=SysData;
        snapshot.state=stateNew;
        snapshot.t=0;
        snapshot.maps=maps;
        snapshot.fault=fault;
        snapshot.diff=diff;
        snapshot.step=0;
        SnapshotLilst{end+1}=snapshot;
        
        SysDataList{end+1}=SysData;
        stateList{end+1}=stateNew;tList{end+1}=0;
        mapsList{end+1}=maps;
        faultList{end+1}=fault;
        diffList{end+1}=diff;stepsList{end+1}=0;
        
        currState=stateNew;
        addLog(['[E',num2str(event(1)),'] Added Syn(s) ',num2str(eventSpec(:,1)')],'INFO');
    elseif eventType==EvtType.DYN_SIM %dyn simulation
        pos=event(5);
        method=event(6);
        AEmethod=round(10*mod(method,1));
        maxTime=event(3)-event(2);
        dt=event(7);
        if ~isempty(SysDataOrig); dt=min([event(7)*50,1.0]);end
        if item>1&&eventList(item-1,4)==EvtType.DYN_SIM; dt=event(7)*10;end
        if dt<SimSet.MIN_DT;dt=SimSet.MIN_DT;end
        if dt>SimSet.MAX_DT;dt=SimSet.MAX_DT;end
        if diff(end)>diffTol;diffTolSim=diff(end)+1e-9;else diffTolSim=diffTol;end
        [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
        
        nbus=size(bus,1);
        nline=size(line,1);
        nInd=size(ind,1);
        nZip=size(zip,1);
        nSyn=size(syn,1);
        nTg=size(tg,1);
        nExc=size(exc,1);
        nAgc=size(agc,1);
        nCac=size(cac,1);
        nCluster=size(cluster,1);
        
        synToTg=zeros(size(syn,1),1);
        synToTg(tg(:,1))=1:size(tg,1);
        synToExc=zeros(size(syn,1),1);
        synToExc(exc(:,1))=1:size(exc,1);
                
        switchMask=eval(switchMaskCondition);
        if pos==0
%             switchMask=switchMask;
            lastDiffCtrl=0;
            if item>1&&eventList(item-1,4)==EvtType.DYN_SIM&&eventList(item-1,6)==MethodType.RK4_NR
                lastDiffCtrl=1;
            end
            if useDiffCtrl&&AEmethod==0&&~lastDiffCtrl&&((diff(end)>diffTolCtrl&&maxTime>1)||notifyCtrlDiff)
                eventListTmp=zeros(size(eventList,1)+1,size(eventList,2));
                eventListTmp(1:item,:)=eventList(1:item,:);
                if eventList(item-1,4)==EvtType.DYN_SIM
                    eventListTmp(item,6)=MethodType.RK4_NR;
                    eventListTmp(item,7)=SimSet.DEFAULT_DT/100;
                    duration=eventListTmp(item,7)*20;
                    if duration>1;duration=1;end
                else
                    duration=maxTime/5;
                    if duration>2;duration=2;end
                end
                eventListTmp((item+1):end,:)=eventList(item:end,:);
                eventListTmp(item,3)=eventListTmp(item,2)+duration;
                eventListTmp(item+1,2)=eventListTmp(item,2)+duration;
                eventList=eventListTmp;
                event=eventList(item,:);
                method=event(6);
                AEmethod=round(10*mod(method,1));
                dt=event(7);
                maxTime=event(3)-event(2);
                diffTolSim=diffTolOrig;
            end
            varOpt.allowModelSwitch=switchMask & allowSteadyDynSwitch;
            varOpt.absT=event(2);
            varOpt.lastSwitchT=lastSwitchT;
            varOpt.switchDelay=switchDelay; 
            varOpt.useDiffCtrl=useDiffCtrl;            
            varOpt.diffTolCtrl=diffTolCtrl;
            varOpt.useNomoveCtrl=useNoMoveCtrl|(~isempty(SysDataOrig)&switchMask & allowSteadyDynSwitch);
            if ~isempty(SysDataOrig)
                SimData=foldSimData(maxTime,5,dt,nlvl,taylorN,alphaTol,diffTolSim,diffTolMax,method,varOpt);
            else
                SimData=foldSimData(maxTime,1,dt,nlvl,taylorN,alphaTol,diffTolSim,diffTolMax,method,varOpt);
            end
            SysPara=foldSysPara([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],fault);
            zipRampSpec=[];
        else
%             switchMask=switchMask;
            ePara=evtDyn(pos,:);
            pqSt=ePara(2);pqEnd=ePara(3);pvSt=ePara(4);pvEnd=ePara(5);
            indSt=ePara(6);indEnd=ePara(7);zipSt=ePara(8);zipEnd=ePara(9);
            shSt=ePara(10);shEnd=ePara(11);zipRampSt=ePara(12);zipRampEnd=ePara(13);
            TmechSt=ePara(14);TmechEnd=ePara(15);PmSt=ePara(16);PmEnd=ePara(17);
            EfSt=ePara(18);EfEnd=ePara(19);VrefSt=ePara(20);VrefEnd=ePara(21);Eq1St=ePara(22);Eq1End=ePara(23);
            if pqSt>0&&pqEnd>0;pqSpec=evtDynPQ(pqSt:pqEnd,:);pqSpec(:,1)=maps.busMap(pqSpec(:,1));pqSpec=pqSpec(pqSpec(:,1)~=0,:);else pqSpec=zeros(0,3);end
            if indSt>0&&indEnd>0;indSpec=evtDynInd(indSt:indEnd,:);indSpec(:,1)=maps.indMap(indSpec(:,1));indSpec=indSpec(indSpec(:,1)~=0,:);else indSpec=zeros(0,2);end
            if zipSt>0&&zipEnd>0;zipSpec=evtDynZip(zipSt:zipEnd,:);zipSpec(:,1)=maps.busMap(zipSpec(:,1));zipSpec=zipSpec(zipSpec(:,1)~=0,:);else zipSpec=zeros(0,2);end
            if shSt>0&&shEnd>0;shSpec=evtDynSh(shSt:shEnd,:);shSpec(:,1)=maps.busMap(shSpec(:,1));shSpec=shSpec(shSpec(:,1)~=0,:);else shSpec=zeros(0,2);end
            if zipRampSt>0&&zipRampEnd>0;zipRampSpec=evtDynZipRamp(zipRampSt:zipRampEnd,:);zipRampSpec(:,1)=maps.busMap(zipRampSpec(:,1));zipRampSpec=zipRampSpec(zipRampSpec(:,1)~=0,:);else zipRampSpec=zeros(0,12);end
            TmechSpecToPm=[];
            TmechSpecToPv=[];
            if TmechSt>0&&TmechEnd>0
                if isempty(SysDataOrig) % full dynamic model
                    TmechSpec=evtDynTmech(TmechSt:TmechEnd,:);
                    TmechSpecIdx=synToTg(maps.synMap(TmechSpec(:,1)));
                    TmechSpecToPm=TmechSpec(TmechSpecIdx==0,:);
                    TmechSpecToPm(:,1)=maps.synMap(TmechSpecToPm(:,1));
                    TmechSpecToPm=TmechSpecToPm(TmechSpecToPm(:,1)~=0,:);
                    TmechSpec=TmechSpec(TmechSpecIdx~=0,:);
                    TmechSpec(:,1)=TmechSpecIdx;
                else % QSS model
                    TmechSpecToPv=evtDynTmech(TmechSt:TmechEnd,:);
                    TmechSpecToPv(:,1)=synMapOrig(TmechSpecToPv(:,1));
                    TmechSpecToPv=TmechSpecToPv(TmechSpecToPv(:,1)~=0,:);
                    TmechSpecToPv(:,1)=SysDataOrig.syn(TmechSpecToPv(:,1),1);
                    TmechSpec=zeros(0,2);
                end
            else
                TmechSpec=zeros(0,2);
            end
            PmSpecToPv=[];
            if PmSt>0&&PmEnd>0
                if isempty(SysDataOrig) % full dynamic model
                    PmSpec=evtDynPm(PmSt:PmEnd,:);
                    PmSpec(:,1)=maps.synMap(PmSpec(:,1));
                    PmSpec=PmSpec(PmSpec(:,1)~=0,:);
                else % QSS model
                    PmSpecToPv=evtDynPm(PmSt:PmEnd,:);
                    PmSpecToPv(:,1)=synMapOrig(PmSpecToPv(:,1));
                    PmSpecToPv=PmSpecToPv(PmSpecToPv(:,1)~=0,:);
                    PmSpecToPv(:,1)=SysDataOrig.syn(PmSpecToPv(:,1),1); 
                    PmSpec=zeros(0,2);
                end
            else
                PmSpec=zeros(0,2);
            end
            if ~isempty(TmechSpecToPm)&&isempty(SysDataOrig)
                PmSpecFull=zeros(size(syn,1),2);
                PmSpecFull(:,1)=1:size(syn,1);
                PmSpecFull(PmSpec(:,1),2)=PmSpec(:,2);
                PmSpecFull(TmechSpecToPm(:,1),2)=PmSpecFull(TmechSpecToPm(:,1),2)+TmechSpecToPm(:,2);
                PmSpec=PmSpecFull(PmSpecFull(:,2)~=0,:);
            end
            
            if pvSt>0&&pvEnd>0
                pvSpec=evtDynPV(pvSt:pvEnd,:);
                pvSpec(:,1)=maps.busMap(pvSpec(:,1));
                pvSpec=pvSpec(pvSpec(:,1)~=0,:);
            else
                pvSpec=zeros(0,2);
            end
            if ~isempty(PmSpecToPv)
                pvSpecFull=zeros(size(bus,1),2);
                pvSpecFull(:,1)=1:size(bus,1);
                pvSpecFull(pvSpec(:,1),2)=pvSpec(:,2);
                pvSpecFull(PmSpecToPv(:,1),2)=PmSpecToPv(:,2);
                pvSpec=pvSpecFull(pvSpecFull(:,2)~=0,:);
            end
            if ~isempty(TmechSpecToPv)
                pvSpecFull=zeros(size(bus,1),2);
                pvSpecFull(:,1)=1:size(bus,1);
                pvSpecFull(pvSpec(:,1),2)=pvSpec(:,2);
                pvSpecFull(TmechSpecToPv(:,1),2)=TmechSpecToPv(:,2);
                pvSpec=pvSpecFull(pvSpecFull(:,2)~=0,:);
            end
            
            if EfSt>0&&EfEnd>0;EfSpec=evtDynEf(EfSt:EfEnd,:);EfSpec(:,1)=maps.synMap(EfSpec(:,1));EfSpec=EfSpec(EfSpec(:,1)~=0,:);else EfSpec=zeros(0,2);end
            if VrefSt>0&&VrefEnd>0;VrefSpec=evtDynVref(VrefSt:VrefEnd,:);VrefSpec(:,1)=synToExc(maps.synMap(VrefSpec(:,1)));VrefSpec=VrefSpec(VrefSpec(:,1)~=0,:);else VrefSpec=zeros(0,2);end
            if Eq1St>0&&Eq1End>0;Eq1Spec=evtDynEq1(Eq1St:Eq1End,:);Eq1Spec(:,1)=maps.synMap(Eq1Spec(:,1));Eq1Spec=Eq1Spec(Eq1Spec(:,1)~=0,:);else Eq1Spec=zeros(0,2);end
            
            [~,Ytr0,Ysh0,~,~,~,~]=getYMatrix(nbus,line,fault);
            
            pqBus=zeros(nbus,size(pqSpec,2)-1);pvBus=zeros(nbus,1);
            pqBus(pqSpec(:,1),:)=pqSpec(:,2:end);pvBus(pvSpec(:,1),:)=pvSpec(:,2);
            pqIncr=pqBus(pq(:,1),:);pvIncr=pvBus(pv(:,1),:);
            zipBus=zeros(nbus,1);zipBus(zipSpec(:,1))=zipSpec(:,2);
            Rzip0=ones(nZip,1);Rzip1=zipBus(zip(:,1),:);
            Rind0=ones(nInd,1);Rind1=zeros(nInd,1);Rind1(indSpec(:,1),:)=indSpec(:,2);
            Ysh1=zeros(nbus,1);Ysh1(shSpec(:,1),:)=shSpec(:,2);
            Tmech1=zeros(nTg,1);Tmech1(TmechSpec(:,1),:)=TmechSpec(:,2);
            Varref1=zeros(nExc,1);Varref1(VrefSpec(:,1),:)=VrefSpec(:,2);
            Ef1=zeros(nSyn,1);Ef1(EfSpec(:,1),:)=EfSpec(:,2);
            Pm1=zeros(nSyn,1);Pm1(PmSpec(:,1),:)=PmSpec(:,2);
            Eq11=zeros(nSyn,1);Eq11(Eq1Spec(:,1),:)=Eq1Spec(:,2);
            zipRampSpec(:,1)=maps.invBusMap(zipRampSpec(:,1));
            
            lastDiffCtrl=0;
            if item>1&&eventList(item-1,4)==EvtType.DYN_SIM&&eventList(item-1,6)==MethodType.RK4_NR
                lastDiffCtrl=1;
            end
            if useDiffCtrl&&AEmethod==0&&~lastDiffCtrl&&((diff(end)>diffTolCtrl&&maxTime>1)||notifyCtrlDiff)
                eventListTmp=zeros(size(eventList,1)+1,size(eventList,2));
                eventListTmp(1:item,:)=eventList(1:item,:);
                if eventList(item-1,4)==EvtType.DYN_SIM
                    eventListTmp(item,6)=MethodType.RK4_NR;
                    eventListTmp(item,7)=SimSet.DEFAULT_DT/100;
                    duration=eventListTmp(item,7)*20;
                    if duration>1;duration=1;end
                else
                    duration=maxTime/5;
                    if duration>2;duration=2;end
                end
                eventListTmp(item,3)=event(2)+duration;
                eventListTmp((item+1):end,:)=eventList(item:end,:);
                eventListTmp(item+1,2)=event(2)+duration;
                eventList=eventListTmp;
                
                event=eventList(item,:);
                if eventList(item+1,5)~=0
                    evtDynAdd=evtDyn(event(5),:);
                    if any(evtDynAdd(:,[6,7,8,9])>0)
                        if all(evtDynAdd(:,[6,7])>0) % Has induction motor item
                            dynIndAdd=evtDynInd(evtDynAdd(6):evtDynAdd(7),:);
                            dynIndAdd(:,2)=dynIndAdd(:,2)./(1+duration*dynIndAdd(:,2));
                            evtDynAdd(6)=size(evtDynInd,1)+1;
                            evtDynAdd(7)=size(evtDynInd,1)+size(dynIndAdd,1);
                            evtDynInd=[evtDynInd;dynIndAdd];
                        end
                        if all(evtDynAdd(:,[8,9])>0)
                            dynZipAdd=evtDynZip(evtDynAdd(8):evtDynAdd(9),:);
                            dynZipAdd(:,2)=dynZipAdd(:,2)./(1+duration*dynZipAdd(:,2));
                            evtDynAdd(8)=size(evtDynZip,1)+1;
                            evtDynAdd(9)=size(evtDynZip,1)+size(dynZipAdd,1);
                            evtDynZip=[evtDynZip;dynZipAdd];
                        end
                        evtDynAdd(:,1)=(size(evtDyn,1)+1):(size(evtDyn,1)+size(evtDynAdd,1));
                        eventList((item+1):(item+size(evtDynAdd,1)),5)=evtDynAdd(:,1);
                        evtDyn=[evtDyn;evtDynAdd];
                    end
                end
                method=event(6);
                AEmethod=round(10*mod(method,1));
                dt=event(7);
                maxTime=event(3)-event(2);
                diffTolSim=diffTolOrig;
            end
            varOpt.allowModelSwitch=switchMask & allowSteadyDynSwitch;
            varOpt.absT=event(2);
            varOpt.lastSwitchT=lastSwitchT;
            varOpt.switchDelay=switchDelay; 
            varOpt.allowExit1=0; 
            varOpt.useDiffCtrl=useDiffCtrl;
            varOpt.diffTolCtrl=diffTolCtrl;
            varOpt.useNomoveCtrl=useNoMoveCtrl|(~isempty(SysDataOrig)&switchMask & allowSteadyDynSwitch);
            if ~isempty(SysDataOrig)
                SimData=foldSimData(maxTime,5,dt,nlvl,taylorN,alphaTol,diffTolSim,diffTolMax,method,varOpt);
            else
                SimData=foldSimData(maxTime,1,dt,nlvl,taylorN,alphaTol,diffTolSim,diffTolMax,method,varOpt);
            end
            SysPara=foldSysPara(pqIncr,pvIncr,Rind0,Rind1,[],[],Rzip0,Rzip1,Ytr0,0*Ytr0,Ysh0,Ysh1,[],[],[],[],[],Tmech1,Varref1,Ef1,Pm1,Eq11,fault);
        end
        if switchMask && allowSteadyDynSwitch && pos==0 && exist('skipTag','var') && skipTag==1
            stateCurve=currState;
            t=0.0;
            finalAlpha=0.0;
            alphaList=maxTime;
            diff=0.0;
            SysDataTmp=SysData;
            exitFlag=DynSimFlag.STEADY;
        else
            skipTag=0;
            [stateCurve,t,finalAlpha,alphaList,diff,SysDataTmp,exitFlag]=generalDynSimulation(SimData,SysData,SysPara,currState,zipRampSpec,maps.busMap);
        end
        
        terminateFlag=0;
        reverseFlag=0;
        switchReady=0;
        if finalAlpha<maxTime-tTol
            if (switchMask && allowSteadyDynSwitch&& (exitFlag==DynSimFlag.QSS || (exitFlag==DynSimFlag.STEADY && pos~=0) || exitFlag==DynSimFlag.NOMOVE_CTRL))||(useDiffCtrl&&exitFlag==DynSimFlag.DIFF_CTRL)
                evtAdd=event;
                curEndTime=event(2)+finalAlpha;
                oriEndTime=event(3);
                event(3)=curEndTime;
                eventList(item,:)=event;
                if item>=size(eventList,1)
                    evtAdd(1)=event(1)+1;
                else
                    evtAdd(1)=(event(1)+eventList(item+1,1))/2;
                end
                evtAdd(2)=curEndTime;
                evtAdd(3)=oriEndTime;
                
                if evtAdd(5)~=0
                    evtDynAdd=evtDyn(evtAdd(5),:);
                    if any(evtDynAdd(:,[6,7,8,9])>0)
                        if all(evtDynAdd(:,[6,7])>0) % Has induction motor item
                            dynIndAdd=evtDynInd(evtDynAdd(6):evtDynAdd(7),:);
                            dynIndAdd(:,2)=dynIndAdd(:,2)./(1+finalAlpha*dynIndAdd(:,2));
                            evtDynAdd(6)=size(evtDynInd,1)+1;
                            evtDynAdd(7)=size(evtDynInd,1)+size(dynIndAdd,1);
                            evtDynInd=[evtDynInd;dynIndAdd];
                        end
                        if all(evtDynAdd(:,[8,9])>0)
                            dynZipAdd=evtDynZip(evtDynAdd(8):evtDynAdd(9),:);
                            dynZipAdd(:,2)=dynZipAdd(:,2)./(1+finalAlpha*dynZipAdd(:,2));
                            evtDynAdd(8)=size(evtDynZip,1)+1;
                            evtDynAdd(9)=size(evtDynZip,1)+size(dynZipAdd,1);
                            evtDynZip=[evtDynZip;dynZipAdd];
                        end
                        evtDynAdd(:,1)=(size(evtDyn,1)+1):(size(evtDyn,1)+size(evtDynAdd,1));
                        evtAdd(:,5)=evtDynAdd(:,1);
                        evtDyn=[evtDyn;evtDynAdd];
                    end
                end                
                
                eventList=[eventList(1:item,:);evtAdd;eventList((item+1):end,:)];
                if exitFlag==DynSimFlag.QSS || (exitFlag==DynSimFlag.STEADY && pos~=0)
                    switchReady=1;
                end
                if exitFlag==DynSimFlag.DIFF_CTRL
                    notifyCtrlDiff=1;
                end
            elseif allowSteadyDynSwitch&& (exitFlag==DynSimFlag.STEADY && pos==0)
                alphaList=[alphaList,maxTime-t(end)];
                diff=[diff,0.0];
                stateCurve=[stateCurve,stateCurve(:,end)];
                t=[t,maxTime];
                skipTag=1;
                addLog(['Simulation directly jumped to maxTime=',num2str(maxTime)],'INFO');
            else            
                if method==0.0||~allowReverse
                    terminateFlag=1;
                else
                    eventList(item,6)=0.0;
                    item=item-1;
                    reverseFlag=1;
                end
            end
        end
        if ~reverseFlag           
            
            SysData=SysDataTmp;
            
            snapshot.SysData=SysData;
            snapshot.state=stateCurve;
            snapshot.t=t;
            snapshot.maps=maps;
            snapshot.fault=fault;
            snapshot.diff=diff;
            snapshot.step=alphaList;
            SnapshotLilst{end+1}=snapshot;
        
            SysDataList{end+1}=SysDataTmp;
            stateList{end+1}=stateCurve;tList{end+1}=t;
            mapsList{end+1}=maps;
            faultList{end+1}=fault;
            diffList{end+1}=diff;stepsList{end+1}=alphaList;
            currState=stateCurve(:,end);
            
            if allowSteadyDynSwitch
                deltax=(stateCurve(:,end)-stateCurve(:,end-1))/(t(end)-t(end-1));
                if switchMask && isempty(SysDataOrig) && exitFlag==DynSimFlag.QSS % Dynamic model -> QSS model
                    if switchReady && eventList(item+1,4)==EvtType.DYN_SIM
                        SysDataOrig=SysData;
                        stateOrig=stateCurve(:,end);
                        mapsOrig=maps;
                        [SysData,currState]=convertSteadyDynModels(SysData,[],currState,deltax,[],SysPara,1);
                        maps.synMap=zeros(size(synBase,1),1);maps.invSynMap=zeros(0,1);
                    else
                        lastSwitchT=event(3);
                    end
                elseif (~isempty(SysDataOrig) &&...
                        (eventList(item+1,4)~=EvtType.DYN_SIM||(item<size(eventList,1)-1&&eventList(item+2,4)==EvtType.DYN_SIM&&eventList(item+2,5)~=0&&eventList(item+2,3)-eventList(item+2,2)<10)))||...
                        exitFlag==DynSimFlag.NOMOVE_CTRL % Steady-state model -> dynamic model
                    [SysData,currState]=convertSteadyDynModels(SysData,SysDataOrig,currState,deltax,stateOrig,SysPara,0);
                    maps=mapsOrig;
                    SysDataOrig=[];
                    lastSwitchT=event(3);
                end
            end
        end
        if eventList(item+1,4)~=EvtType.DYN_SIM
            skipTag=0;
            lastSwitchT=event(3);
            notifyCtrlDiff=0;
        end
        if terminateFlag
            msg=['[',num2str(event(2)+finalAlpha),']','Simulation terminated earlier. Exit'];
            flag=-1;
            addLog(['Simulation terminated earlier. Exit'],'INFO');
            break;
        end
    elseif eventType==EvtType.FAULT % fault
        pos=event(5);
        method=event(6);
        evtFaultSt=evtFault(pos,2);evtFaultEnd=evtFault(pos,3);
        eventSpec=evtFaultSpec(evtFaultSt:evtFaultEnd,:);
        
        addFault=eventSpec(eventSpec(:,5)==0,1:4);
        clearFault=eventSpec(eventSpec(:,5)~=0,1:4);
        addFault(:,1)=maps.lineMap(addFault(:,1));
        clearFault(:,1)=maps.lineMap(clearFault(:,1));
        
        faultPrev=fault;
        faultTag=zeros(size(maps.lineMap));
        faultTag(maps.invLineMap(fault(:,1)))=1;
        faultTag(maps.invLineMap(clearFault(:,1)))=0;
        fault=fault(faultTag(maps.invLineMap(fault(:,1)))==1,:);
        addFault=addFault(faultTag(maps.invLineMap(addFault(:,1)))~=1,:);
        fault=[fault;addFault];
        
        if diff(end)>diffTol;diffTolSim=diff(end)+1e-9;else diffTolSim=diffTol;end
        [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
        nbus=size(bus,1);
        busTag=ones(nbus,1);
        
        SimData=foldSimData([],1,[],nlvl,taylorN,alphaTol,diffTolSim,diffTolMax,method);
        doContinue=1;
        while 1
            SysPara=foldSysPara([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],fault);
            [stateNew,finalAlpha,alphaList,diff,SysDataTemp,faultPrev,fault,mapsTemp]=...
                simulationChangeFault(busTag,SysData,SysDataBase,SysPara,SimData,maps,currState,faultPrev,fault);
            if finalAlpha<1-tTol
                [flag,doContinue,willBreak,msg]=processCollaps(SysDataTemp,stateNew(:,end),event,'Changing fault status',flag,doContinue);
                if willBreak
                    break; 
                end
            else
                SysData=SysDataTemp;
                maps=mapsTemp;
                break;
            end
        end
        if ~doContinue
            break;
        end
        
        snapshot.SysData=SysData;
        snapshot.state=stateNew;
        snapshot.t=0;
        snapshot.maps=maps;
        snapshot.fault=fault;
        snapshot.diff=diff;
        snapshot.step=0;
        SnapshotLilst{end+1}=snapshot;
        
        SysDataList{end+1}=SysData;
        stateList{end+1}=stateNew;tList{end+1}=0;
        mapsList{end+1}=maps;
        faultList{end+1}=fault;
        diffList{end+1}=diff;stepsList{end+1}=0;
        
        currState=stateNew;
        addLog(['[E',num2str(event(1)),'] Applied fault(s) ',num2str(eventSpec(:,1)')],'INFO');
    elseif eventType==EvtType.CUT_LINE %
        pos=event(5);
        method=event(6);
        evtLineSt=evtLine(pos,2);evtLineEnd=evtLine(pos,3);
        eventSpec=evtLineSpec(evtLineSt:evtLineEnd,:);
        eventSpec=eventSpec(eventSpec(:,2)==1,:);
        lineCutIdx=eventSpec(:,1);
        
        if diff(end)>diffTol;diffTolSim=diff(end)+1e-9;else diffTolSim=diffTol;end
        [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
        nbus=size(bus,1);
        busTag=ones(nbus,1);
        
        SimData=foldSimData([],1,[],nlvl,taylorN,alphaTol,diffTolSim,diffTolMax,method);
        doContinue=1;
        SysPara=foldSysPara([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],fault);
        while 1
            [stateNew,finalAlpha,alphaList,diff,SysDataTemp,SysPara,mapsTemp]=...
                simulationCutLine(busTag,SysData,SysDataBase,SysPara,SimData,maps,currState,lineCutIdx);
            if finalAlpha<1-tTol
                [flag,doContinue,willBreak,msg]=processCollaps(SysDataTemp,stateNew(:,end),event,'Cut line',flag,doContinue);
                if willBreak
                    break; 
                end
            else
                SysData=SysDataTemp;
                maps=mapsTemp;
                break;
            end
        end
        if ~doContinue
            break;
        end
        snapshot.SysData=SysData;
        snapshot.state=stateNew;
        snapshot.t=0;
        snapshot.maps=maps;
        snapshot.fault=fault;
        snapshot.diff=diff;
        snapshot.step=0;
        SnapshotLilst{end+1}=snapshot;
        
        SysDataList{end+1}=SysData;
        stateList{end+1}=stateNew;tList{end+1}=0;
        mapsList{end+1}=maps;
        faultList{end+1}=fault;
        diffList{end+1}=diff;stepsList{end+1}=0;
        
        currState=stateNew;
        addLog(['[E',num2str(event(1)),'] Cut line(s) ',num2str(lineCutIdx')],'INFO');
    elseif eventType==EvtType.CUT_LOAD % Not necessary, cut load is adding negative load.
        
    elseif eventType==EvtType.CUT_IND %
        pos=event(5);
        method=event(6);
        evtIndSt=evtInd(pos,2);evtIndEnd=evtInd(pos,3);
        eventSpec=evtIndSpec(evtIndSt:evtIndEnd,:);
        eventSpecMod=eventSpec(eventSpec(:,2)==1,:);
        eventSpec=eventSpec(eventSpec(:,2)==2,:);
        
        if diff(end)>diffTol;diffTolSim=diff(end)+1e-9;else diffTolSim=diffTol;end
        [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
        nbus=size(bus,1);
        busTag=ones(nbus,1);
        
        SimData=foldSimData([],1,[],nlvl,taylorN,alphaTol,diffTolSim,diffTolMax,method);
        doContinue=1;
        SysPara=foldSysPara([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],fault);
        while 1
            [stateNew,finalAlpha,alphaList,diff,SysDataTemp,SysPara,mapsTemp]=...
                simulationCutMotorInstant(busTag,SysData,SysDataBase,SysPara,SimData,maps,currState,eventSpec(:,1));
            if finalAlpha<1-tTol
                [flag,doContinue,willBreak,msg]=processCollaps(SysDataTemp,stateNew(:,end),event,'Cut induction motor',flag,doContinue);
                if willBreak
                    break; 
                end
            else
                SysData=SysDataTemp;
                maps=mapsTemp;
                break;
            end
        end
        if ~doContinue
            break;
        end
        if ~isempty(eventSpecMod)
            [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
            indModIdx=maps.indMap(eventSpecMod(:,1));
            eventSpecMod=eventSpecMod(indModIdx~=0,:);
            indModIdx=indModIdx(indModIdx~=0);
            
            ind(indModIdx,15)=eventSpecMod(:,3);
            ind(indModIdx,16:17)=0;
            SysData=foldSysData(bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster);
        end
        
        snapshot.SysData=SysData;
        snapshot.state=stateNew;
        snapshot.t=0;
        snapshot.maps=maps;
        snapshot.fault=fault;
        snapshot.diff=diff;
        snapshot.step=0;
        SnapshotLilst{end+1}=snapshot;
        
        SysDataList{end+1}=SysData;
        stateList{end+1}=stateNew;tList{end+1}=0;
        mapsList{end+1}=maps;
        faultList{end+1}=fault;
        diffList{end+1}=diff;stepsList{end+1}=0;
        
        currState=stateNew;
        addLog(['[E',num2str(event(1)),'] Cut ind(s) ',num2str(eventSpec(:,1)')],'INFO');
    elseif eventType==EvtType.CUT_SYN %
        pos=event(5);
        method=event(6);
        evtSynSt=evtSyn(pos,2);evtSynEnd=evtSyn(pos,3);
        eventSpec=evtSynSpec(evtSynSt:evtSynEnd,:);
        eventSpec=eventSpec(eventSpec(:,2)==1,:);
        
        if diff(end)>diffTol;diffTolSim=diff(end)+1e-9;else diffTolSim=diffTol;end
        [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
        nbus=size(bus,1);
        busTag=ones(nbus,1);
        
        SimData=foldSimData([],1,[],nlvl,taylorN,alphaTol,diffTolSim,diffTolMax,method);
        doContinue=1;
        SysPara=foldSysPara([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],fault);
        while 1
            [stateNew,finalAlpha,alphaList,diff,SysDataTemp,SysPara,mapsTemp]=...
                simulationCutSynInstant(busTag,SysData,SysDataBase,SysPara,SimData,maps,...
                currState,eventSpec(:,1));
            
            if finalAlpha<1-tTol
                [flag,doContinue,willBreak,msg]=processCollaps(SysDataTemp,stateNew(:,end),event,'Cut synchronous generator',flag,doContinue);
                if willBreak
                    break; 
                end
            else
                SysData=SysDataTemp;
                maps=mapsTemp;
                break;
            end
        end
        if ~doContinue
            break;
        end
        SysDataList{end+1}=SysData;
        stateList{end+1}=stateNew;tList{end+1}=0;
        mapsList{end+1}=maps;
        faultList{end+1}=fault;
        diffList{end+1}=diff;stepsList{end+1}=0;
        
        currState=stateNew;
        addLog(['[E',num2str(event(1)),'] Cut syn(s) ',num2str(eventSpec(:,1)')],'INFO');
    elseif eventType==EvtType.ALT_SYN %
        pos=event(5);
        method=event(6);
        evtSynSt=evtSyn(pos,2);evtSynEnd=evtSyn(pos,3);
        eventSpec=evtSynSpec(evtSynSt:evtSynEnd,:);
        eventSpec=eventSpec(eventSpec(:,2)==2,:);
        if ~isempty(eventSpec)
            dSynCompact=eventSpec(:,3:9);            
            newSynIdx=maps.synMap(eventSpec(:,1));
            activeSynIdx=newSynIdx(newSynIdx~=0);
            deltaSyn=zeros(size(syn,1),size(dSynCompact,2));
            deltaSyn(activeSynIdx,:)=dSynCompact;
        else
            deltaSyn=zeros(size(syn,1),7);
        end
        
        if diff(end)>diffTol;diffTolSim=diff(end)+1e-9;else diffTolSim=diffTol;end
        [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
        nbus=size(bus,1);
        busTag=ones(nbus,1);
        
        SimData=foldSimData([],1,[],nlvl,taylorN,alphaTol,diffTolSim,diffTolMax,method);
        doContinue=1;
        SysPara=foldSysPara([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],fault);
        synNew=syn;
        synNew(:,[7:10,13:15])=synNew(:,[7:10,13:15])+deltaSyn;
        while 1
            [stateNew,finalAlpha,alphaList,diff,SysDataTemp,SysPara,mapsTemp]=...
                simulationAlterSynParameter(busTag,SysData,SysDataBase,SysPara,SimData,maps,currState,synNew);
            if finalAlpha<1-tTol
                [flag,doContinue,willBreak,msg]=processCollaps(SysDataTemp,stateNew(:,end),event,'Alter syn parameters',flag,doContinue);
                if willBreak
                    break; 
                end
            else
                SysData=SysDataTemp;
                maps=mapsTemp;
                break;
            end
        end
        if ~doContinue
            break;
        end
        
        snapshot.SysData=SysData;
        snapshot.state=stateNew;
        snapshot.t=0;
        snapshot.maps=maps;
        snapshot.fault=fault;
        snapshot.diff=diff;
        snapshot.step=0;
        SnapshotLilst{end+1}=snapshot;
        
        SysDataList{end+1}=SysData;
        stateList{end+1}=stateNew;tList{end+1}=0;
        mapsList{end+1}=maps;
        faultList{end+1}=fault;
        diffList{end+1}=diff;stepsList{end+1}=0;
        
        currState=stateNew;
        addLog(['[E',num2str(event(1)),'] Alter syn(s) ',num2str(eventSpec(:,1)')],'INFO');
    elseif eventType==EvtType.END % end
        msg=['[',num2str(event(2)),'] ','Simulation successfully finished.'];
        flag=1;
        
        snapshot.SysData=SysData;
        snapshot.state=currState;
        snapshot.t=0;
        snapshot.maps=maps;
        snapshot.fault=fault;
        snapshot.diff=diff;
        snapshot.step=0;
        break;
    else
        % Undefined event types, do nothing
    end
    item=item+1;
end

timeSpent=toc(startTag);
addLog(['Total computation time:',num2str(timeSpent),'s.'],'INFO');
% disp(tp)

[t,stateCurve]=reformatCurves(SysDataBase,tList,stateList,SysDataList,mapsList);

clearAllTempFiles();

if options.output
    save(['restoration/',caseName,'_',timestamp,'_simp.mat'],...
        'caseName','simSettings','eventList','SysDataBase','t','stateCurve','timeSpent','comment',...
        'snapshot','-v7.3');
    addLog(['Curves saved at: ','restoration/',caseName,'_',timestamp,'_simp.mat'],'INFO');
end
if options.output&&~options.simpleOutput
    save(['restoration/',caseName,'_',timestamp,'.mat'],...
        'caseName','simSettings','eventList','SysDataBase',...
        'SnapshotLilst','tList','stateList',...
        'timeSpent','comment','snapshot','-v7.3');
    addLog(['More detailed results saved at: ','restoration/',caseName,'_',timestamp,'.mat'],'INFO');
end

addLog('================== END ===================','INFO');
flushLogs();

res.flag=flag;
res.msg=msg;
res.t=t;
res.stateCurve=stateCurve;
res.simSettings=simSettings;
res.SysDataBase=SysDataBase;
res.caseName=caseName;
res.eventList=eventList;
res.timestamp=timestamp;
res.snapshot=snapshot;

end

function [flag,doContinue,willBreak,msg]=processCollaps(SysData,state,event,actionStr,flag,doContinue)
    willBreak=0;   
    msg='';

    [bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
    [V,~,s,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~]=unfoldX(state(:,end),SysData);
    [jac,jPIdx,jQIdx,jVIdx,jSIdx,busType]=getJacobianMod(size(bus,1),size(line,1),bus,sw,pv,pq,shunt,line,zip,ind,s,V,ones(size(ind,1)));
    [Vx,Dx]=eig(full(jac));
    dDx=diag(Dx);
    [~,iSmall]=min(abs(dDx));
    part=Vx(:,iSmall).*conj(Vx(:,iSmall));
    partBus=zeros(size(bus,1),1);
    partBus(busType~=2)=partBus(busType~=2)+part(jPIdx);
    partBus(busType==0)=partBus(busType==0)+part(jQIdx);
    [nIslands,islands,refs]=searchIslands(bus(:,1),line(:,1:2));
    if nIslands==1
        msg=['[',num2str(event(2)),'] ',actionStr,' causes singularity. Exit'];
        flag=-1;
        addLog([actionStr,' causes singularity. Exit'],'INFO');
        doContinue=0;
        willBreak=1;
    else
        addLog(['Part of the system collapses, purge collapsed parts.'],'INFO');
        partIsland=accumarray(islands,partBus,[nIslands,1]);
        [maxPartIsland,iMaxPart]=max(partIsland);
        if maxPartIsland>0.8
            busFilter=find(busTag==1);
            busTag(busFilter(islands==iMaxPart(1)))=0;
        else
            msg=['[',num2str(event(2)),'] ','No dominant collapse island can be purged. Exit'];
            flag=-1;
            addLog(['No dominant collapse island can be purged. Exit'],'INFO');
            doContinue=0;
            willBreak=1;
        end
    end
end

function options=regulateOptions(options)
if nargin<1
    options=[];
end
if ~isfield(options,'dbstop');options.dbstop=0;end
if ~isfield(options,'output');options.output=0;end
if ~isfield(options,'comment');options.comment='';end
if ~isfield(options,'tTol');options.tTol=1e-9;end
if ~isfield(options,'allowReverse');options.allowReverse=0;end
if ~isfield(options,'allowSteadyDynSwitch');options.allowSteadyDynSwitch=0;end
if ~isfield(options,'useDiffCtrl');options.useDiffCtrl=0;end
if ~isfield(options,'useNoMoveCtrl');options.useNoMoveCtrl=0;end
if ~isfield(options,'nPool');options.nPool=1;end
if ~isfield(options,'hotStart');options.hotStart=0;end
if ~isfield(options,'switchDelay');options.switchDelay=3.0;end
if ~isfield(options,'notifyCtrlDiff');options.notifyCtrlDiff=0;end
if ~isfield(options,'inheritParaSettings');options.inheritParaSettings=1;end
if ~isfield(options,'nlvl');options.nlvl=15;end
if ~isfield(options,'taylorN');options.taylorN=4;end
if ~isfield(options,'segAlpha');options.segAlpha=1;end
if ~isfield(options,'dAlpha');options.dAlpha=1;end
if ~isfield(options,'alphaTol');options.alphaTol=1.0000e-03;end
if ~isfield(options,'diffTol');options.diffTol=1.0000e-06;end
if ~isfield(options,'diffTolMax');options.diffTolMax=1.0000e-02;end
if ~isfield(options,'method');options.method=0;end
if ~isfield(options,'diffTolCtrl');options.diffTolCtrl=1.0000e-05;end
if ~isfield(options,'simpleOutput');options.simpleOutput=1;end
end