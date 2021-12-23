% EXAMPLE: extended-term simulation
% 
% in Software PowerSAS.m
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

clear; clc;

global Global_Settings

Global_Settings.logLevel='INFO';

% Full dynamic simulation 
tagFullDynStart=tic;
res_004_fulldyn=runPowerSAS('dyn','d_004_2a_bs_agc.m',setOptions('allowSteadyDynSwitch',0),'settings_d_004_2a_agc');
timeFullDyn=toc(tagFullDynStart);

clc;
% Hybrid simulation with dynamic-QSS switching
tagHybridStart=tic;
res_004=runPowerSAS('dyn','d_004_2a_bs_agc.m',setOptions('allowSteadyDynSwitch',1),'settings_d_004_2a_agc');
timeHybrid=toc(tagHybridStart);

clc;

disp(['Full dynamic simulation computation time:', num2str(timeFullDyn),' s.']);
disp(['Hybrid simulation computation time:', num2str(timeHybrid),' s.']);

plotCurves(1,res_004_fulldyn.t,res_004_fulldyn.stateCurve,res_004_fulldyn.SysDataBase,'v');
plotCurves(2,res_004.t,res_004.stateCurve,res_004.SysDataBase,'v');