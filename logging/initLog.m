function initLog(fileName,maxLog,mode,forceRedirect)
% Initialize log system
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
%	fileName - The name of log file
%	maxLog - The volume of log buffer
%	mode
%		DEBUG
%		INFO
%		WARN
%		ERROR
%	forceRedirect - Force changing to another file
%

global logBuffer logFileName logBufferPos maxLogBuffer logMode;

if checkLogReady()
    addLog('Logging is already set but get new initial logging request. Flush old logs before transferring to new logging.','DEBUG');
end

if nargin<4
    forceRedirect=0;
end
flushLogs();
maxLogBuffer=maxLog;
logBuffer=cell(maxLogBuffer,1);
logBufferPos=1;
if forceRedirect
    logFileName=fileName;
    logId=fopen(logFileName,'a');
    fclose(logId);
else
    if ~isValidLogFile(logFileName)
        logFileName=fileName;
    end
    logId=fopen(logFileName,'a');
    fclose(logId);
end
logMode=mode;
flushLogs();
end

function valid=isValidLogFile(fileName)
valid= ~isempty(fileName)&&exist(fileName,'file');
end