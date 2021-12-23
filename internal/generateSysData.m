function SysData=generateSysData(dataFile,dataPath)
% Generate system data from data file
%
% FUNCTION generateSysData
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
%	dataFile - name of data file
%	dataPath - path of data file
%
% OUTPUT
%	SysData - system data
%

homePath=pwd;
cd(dataPath)
eval(dataFile)
if ~exist('Bus','var')||isempty(Bus)||isempty(Bus.con);busBase=zeros(0,6);else busBase=Bus.con;end
if ~exist('PV','var')||isempty(PV)||isempty(PV.con);pvBase=zeros(0,11);else pvBase=PV.con;end
if ~exist('PQ','var')||isempty(PQ)||isempty(PQ.con);pqBase=zeros(0,9);else pqBase=PQ.con;end
if ~exist('SW','var')||isempty(SW)||isempty(SW.con);swBase=zeros(0,13);else swBase=SW.con;end
if ~exist('Line','var')||isempty(Line)||isempty(Line.con);lineBase=zeros(0,16);else lineBase=Line.con;end
if ~exist('Ltc','var')||isempty(Ltc)||isempty(Ltc.con);ltcBase=zeros(0,16);else ltcBase=Ltc.con;end
if ~exist('Supply','var')||isempty(Supply)||isempty(Supply.con);supBase=zeros(0,20);else supBase=Supply.con;end
if ~exist('Demand','var')||isempty(Demand)||isempty(Demand.con);demBase=zeros(0,18);else demBase=Demand.con;end
if ~exist('Shunt','var')||isempty(Shunt)||isempty(Shunt.con);shuntBase=zeros(0,7);else shuntBase=Shunt.con;end
if ~exist('Ind','var')||isempty(Ind)||isempty(Ind.con);indBase=zeros(0,20);else indBase=Ind.con;end
if ~exist('Pl','var')||isempty(Pl)||isempty(Pl.con);zipBase=zeros(0,12);else zipBase=Pl.con;end
if ~exist('Syn','var')||isempty(Syn)||isempty(Syn.con);synBase=zeros(0,26);else synBase=Syn.con;end
if ~exist('Exc','var')||isempty(Exc)||isempty(Exc.con);excBase=zeros(0,14);else excBase=Exc.con;end
if ~exist('Tg','var')||isempty(Tg)||isempty(Tg.con);tgBase=zeros(0,8);else tgBase=Tg.con;end
if ~exist('Bus','var')||isempty(Bus)||isempty(Bus.names);busNameBase={};else busNameBase=Bus.names;end
if ~exist('Agc','var')||isempty(Agc)||isempty(Agc.con);agcBase=zeros(0,4);else agcBase=Agc.con;end
if ~exist('CAC','var')||isempty(CAC)||isempty(CAC.con);cacBase=zeros(0,10);else cacBase=CAC.con;end
if ~exist('Cluster','var')||isempty(Cluster)||isempty(Cluster.con);clusterBase=zeros(0,10);else clusterBase=Cluster.con;end
cd(homePath);

SysData=foldSysData(busBase,swBase,pvBase,pqBase,shuntBase,lineBase,indBase,zipBase,synBase,excBase,tgBase,agcBase,cacBase,clusterBase);
end