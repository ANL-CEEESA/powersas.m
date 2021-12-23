function res=calCPF(SysData,options,caseName,varargin)

%
% FUNCTION calCPF
%
% Function for calculating continuation power flow
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

res=[];

simSettings=getDefaultSimSettings();

if nargin<2
    options=regulateOptions();
else
    options=regulateOptions(options);
end

if nargin<3||isempty(caseName)
    caseName=['d_~',num2str(size(SysData.bus,1))];
end

nVarargs=length(varargin);
if nVarargs<1
    ME=MException('PowerSAS:invalidArgument','CPF should have arguments specifying load increase direction.');
    throw(ME);
else
    zipRamp=[];
    pvRamp=[];
    if nVarargs>=1
        zipRamp=varargin{1};
    end
    if nVarargs>=2
        pvRamp=varargin{2};
    end
    if nVarargs>=3
        snapshot=varargin{3};
        hotStart=1;
    else
        snapshot=[];
        hotStart=0;
    end
    simSettings.eventList=[simSettings.eventList;[2   0.0000  options.simLen   50    1 0.0  0.0000]];
    simSettings.evtDynZipRamp=zipRamp;
    simSettings.evtDynPV=pvRamp;
    simSettings.evtDyn=zeros(1,23); simSettings.evtDyn(1,1)=1;
    simSettings.evtDyn(1,4)=1; simSettings.evtDyn(1,5)=size(pvRamp,1);
    simSettings.evtDyn(1,12)=1; simSettings.evtDyn(1,13)=size(zipRamp,1);
    options.hotStart=hotStart;
    
    res=runDynamicSimulationExec(caseName,SysData,simSettings,options,snapshot);
end

end


function options=regulateOptions(options)
if nargin<1
    options=[];
end

if ~isfield(options,'nlvl');options.nlvl=15;end
if ~isfield(options,'taylorN');options.taylorN=4;end
if ~isfield(options,'segAlpha');options.segAlpha=1;end
if ~isfield(options,'dAlpha');options.dAlpha=1;end
if ~isfield(options,'alphaTol');options.alphaTol=1e-3;end
if ~isfield(options,'diffTol');options.diffTol=1e-6;end
if ~isfield(options,'diffTolMax');options.diffTolMax=1e-2;end
if ~isfield(options,'method');options.method=0;end
if ~isfield(options,'simLen');options.simLen=200;end
if ~isfield(options,'diffTolCtrl');options.diffTolCtrl=1e-5;end
if ~isfield(options,'useNoMoveCtrl');options.useNoMoveCtrl=1;end
if ~isfield(options,'Efstd');options.Efstd=1.2;end
end