% EXAMPLE: steady-state
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

clear;
clc;

global Global_Settings

Global_Settings.logLevel='INFO';

res_003=runPowerSAS('pf','d_003.m');
res_014_syn_ind_zip=runPowerSAS('pf','d_014_syn_ind_zip');

res_npcc140_ind_zip=runPowerSAS('pf','d_npcc140_ind_zip.m');
res_wecc=runPowerSAS('pf','d_wecc.m');
res_dcase2383wp_mod_ind_zip=runPowerSAS('pf','d_2383wp_mod_ind_zip.m');

% Uncomment the lines below to try running large cases
% res_ACTIVSg70k=runPowerSAS('pf','d_ACTIVSg70k.m',setOptions('dataPath',[pwd,'/data'],'nlvl',25,'segAlpha',0.4,'diffTol',1e-4));
% res_210k=runPowerSAS('pf','d_210k.m',setOptions('dataPath',[pwd,'/data'],'nlvl',25,'segAlpha',0.4,'diffTol',1e-4));
