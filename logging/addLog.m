function addLog(info,level,dispConfig)
% Add log 
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
% INPUT
%	info - The log message
%	level
%		DEBUG
%		INFO
%		WARN
%		ERROR
%	dispConfig - If 'NONE', then do not display the timestamp
%
global logBuffer logBufferPos maxLogBuffer logMode Global_Settings;

logOK=checkLogReady();
if ~logOK
    initLog([pwd,'/~default_log.txt'],1000,'apparent');
end

if exist('Global_Settings','var')&&~isempty(Global_Settings)
    if isfield(Global_Settings,'allosUserLogLevel')
        allosUserLogLevel=Global_Settings.allosUserLogLevel;
    else
        allosUserLogLevel=0;
    end
    if isfield(Global_Settings,'logLevel')
        logLevel=Global_Settings.logLevel;
    else
        logLevel=0;
    end
else
    allosUserLogLevel=0;
    logLevel='INFO';
end
levelCode=getLogLevelCode(level,allosUserLogLevel);
configLevelCode=getLogLevelCode(logLevel);

if levelCode>=configLevelCode
    if nargin>=3&&~isempty(dispConfig)&&strcmp(dispConfig,'NONE')
        msg=info;
    else
        msg=[datestr(clock,13),'[',level,']',info];
    end
    logBuffer{logBufferPos,1}=msg;
    if strcmp(logMode,'apparent')
        disp(msg);
    end
    if logBufferPos==maxLogBuffer
        flushLogs();
    else
        logBufferPos=logBufferPos+1;
    end
end

end

function levelCode=getLogLevelCode(level,allowUd)
if strcmp(level,'DEBUG');
    levelCode=0;
elseif strcmp(level,'INFO')
    levelCode=1;
elseif strcmp(level,'WARN')
    levelCode=2;
elseif strcmp(level,'ERROR')
    levelCode=3;
else
    if nargin>=2&&allowUd
        levelCode=4;
    else
        levelCode=-1;
    end
end
end