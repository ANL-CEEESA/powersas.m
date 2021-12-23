function [nState,idxs]...
    =getIndexDyn(SysData)
% Generate the indexes of state variables from system data
%
% FUNCTION getIndexDyn
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

nbus=size(SysData.bus,1);
nInd=size(SysData.ind,1);
nSyn=size(SysData.syn,1);
nExc=size(SysData.exc,1);
nTg=size(SysData.tg,1);
nCac=size(SysData.cac,1);
nCluster=size(SysData.cluster,1);

idxs.vIdx=(1:nbus)';
idxs.qIdx=(1:nbus)'+nbus;
idxs.sIdx=(1:nInd)'+2*nbus;
idxs.deltaIdx=(1:nSyn)'+2*nbus+nInd;
idxs.omegaIdx=(1:nSyn)'+2*nbus+nInd+nSyn;
idxs.eq1Idx=(1:nSyn)'+2*nbus+nInd+2*nSyn;
idxs.eq2Idx=(1:nSyn)'+2*nbus+nInd+3*nSyn;
idxs.ed1Idx=(1:nSyn)'+2*nbus+nInd+4*nSyn;
idxs.ed2Idx=(1:nSyn)'+2*nbus+nInd+5*nSyn;
idxs.psidIdx=(1:nSyn)'+2*nbus+nInd+6*nSyn;
idxs.psiqIdx=(1:nSyn)'+2*nbus+nInd+7*nSyn;
idxs.pgIdx=(1:nSyn)'+2*nbus+nInd+8*nSyn;
idxs.efIdx=(1:nSyn)'+2*nbus+nInd+9*nSyn;
idxs.vavrmIdx=(1:nExc)'+2*nbus+nInd+10*nSyn;
idxs.vavrrIdx=(1:nExc)'+2*nbus+nInd+10*nSyn+nExc;
idxs.vavrfIdx=(1:nExc)'+2*nbus+nInd+10*nSyn+2*nExc;
idxs.vavrrefIdx=(1:nExc)'+2*nbus+nInd+10*nSyn+3*nExc;
idxs.tgovgIdx=(1:nTg)'+2*nbus+nInd+10*nSyn+4*nExc;
idxs.tgovmIdx=(1:nTg)'+2*nbus+nInd+10*nSyn+4*nExc+nTg;
idxs.tmechIdx=(1:nTg)'+2*nbus+nInd+10*nSyn+4*nExc+2*nTg;

idxs.fIdx=(1:nbus)'+2*nbus+nInd+10*nSyn+4*nExc+3*nTg;
idxs.dpgIdx=(1:nbus)'+3*nbus+nInd+10*nSyn+4*nExc+3*nTg;
idxs.qpltIdx=(1:nCac)'+4*nbus+nInd+10*nSyn+4*nExc+3*nTg;
idxs.vgIdx=(1:nCluster)'+4*nbus+nInd+10*nSyn+4*nExc+3*nTg+nCac;

nState=4*nbus+nInd+10*nSyn+4*nExc+3*nTg+nCac+nCluster;
end