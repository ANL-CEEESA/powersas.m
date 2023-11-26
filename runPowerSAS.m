function res=runPowerSAS(simType,data,options,varargin)

%
% FUNCTION runPowerSAS
%
% Top-level function for calling PowerSAS
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

initpowersas();

if nargin<2
    %error
    ME=MException('PowerSAS:invalidArgument','Argument ''data'' cannot be absent.');
    throw(ME);
end

%% prepare options
if nargin<3
    options=regulateOptions();
else
    options=regulateOptions(options);
end
if options.dbstop
    dbstop if error
end

if ~exist(options.outputPath,'dir');mkdir(options.outputPath);end

timestamp=datestr(clock,'yyyymmddTHHMMSSFFF');
resultPath=[options.outputPath,'/',timestamp];
if options.output
    mkdir(resultPath);
    initLog([resultPath,'/',timestamp,'_','log.txt'],1000,'apparent');
else
    initLog([options.outputPath,'/~default_log.txt'],1000,'apparent');
    addLog(['OUTPUT IS DISABLED. Check ',options.outputPath,'/~default_log.txt for logs.'],'DEBUG');
end
addLog(['The session ID is ',timestamp,'.'],'INFO');
%khuang 919
options.timestamp=timestamp;

homePath=pwd;
%% prepare data 
caseName=[];
if isstruct(data)
    SysData=data;
elseif ischar(data)
    dataFile=strtrim(data);
    if strcmp(dataFile(end-1:end),'.m')
        dataFile=dataFile(1:end-2);
    end 
    delim=strfind(dataFile,'/');
    if ~isempty(delim)
        dataDir=dataFile(1:delim(end));
        options.dataPath=dataDir;
        dataFile=dataFile(delim(end)+1:end);
    end
    
    caseName=dataFile;
    
    try
        addLog(['CASE ',dataFile,'. Timestamp=',timestamp,'.'],'INFO');
        SysData=readDataFile(dataFile,options.dataPath,homePath);
        addLog(['Finished loading data file. '],'INFO');
    catch ME
        cd(homePath);
        addLog(['Read and parse data file error.'],'ERROR');
        rethrow(ME);
    end
else
    %error
    ME=MException('PowerSAS:invalidArgument','Argument ''data'' should be SysData struct or string.');
    throw(ME);
end

%% main body 
if strcmp(simType,'pf')
    res=calPF(SysData,options,caseName);
elseif strcmp(simType,'emt')
    res=calEMT(data);
elseif strcmp(simType,'cpf')
    res=calCPF(SysData,options,caseName,varargin{:});    
elseif strcmp(simType,'tsa')
    res=calTSA(SysData,options,caseName,varargin{:});
elseif strcmp(simType,'ctg')
    res=calCTG(SysData,options,caseName,varargin{:});    
elseif strcmp(simType,'n-1')
    resPf=calPF(SysData,options,caseName);
    options.hotStart=1;
    res=cell(size(SysData.line,1),1);
    for lineIdx=1:size(SysData.line,1)
        try
            res{lineIdx}=calCTG(SysData,options,caseName,lineIdx,resPf.snapshot); 
        catch ME
            
        end
    end
elseif strcmp(simType,'dyn')
    nVarargs=length(varargin);
    if nVarargs<1
        ME=MException('PowerSAS:invalidArgument','Dynamic simulation should have arguments specifying simulation settings.');
        throw(ME);
    else
        if nVarargs>=1
            settings=varargin{1};
        end
        if nVarargs>=2
            snapshot=varargin{2};
        else
            snapshot=[];
        end
        
        if isstruct(settings)
            simSettings=settings;
        elseif ischar(settings)
            simSettingFile=strtrim(settings);
            if strcmp(simSettingFile(end-1:end),'.m')
                simSettingFile=simSettingFile(1:end-2);
            end
            try
                addLog(['SETTINGS ',simSettingFile,'.'],'INFO');
                simSettings=readSimulationSettings(simSettingFile,options.settingPath,homePath);
            catch ME
                cd(homePath);
                addLog(['Read and parse settings file error.'],'ERROR');
                rethrow(ME);
            end
        else
            %error
            ME=MException('PowerSAS:invalidArgument','Argument ''settings'' should be SimSettings struct or string.');
            throw(ME);
        end        
        
        res=runDynamicSimulationExec(caseName,SysData,simSettings,options,snapshot);
    end
else
    
end
addLog('','INFO');
flushLogs();

end


function options=regulateOptions(options)
if nargin<1
    options=[];
end
if ~isfield(options,'dataPath');options.dataPath=[pwd,'/data'];end
if ~isfield(options,'settingPath');options.settingPath=options.dataPath;end
if ~isfield(options,'dbstop');options.dbstop=0;end
if ~isfield(options,'output');options.output=0;end
if ~isfield(options,'outputPath');options.outputPath=[pwd,'/output'];end
end
