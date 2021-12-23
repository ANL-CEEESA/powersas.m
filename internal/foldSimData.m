function SimData=foldSimData(maxAlpha,segAlpha,dAlpha,nlvl,taylorN,alphaTol,diffTol,diffTolMax,method,varOpt)
% Fold simulation settings to SimData structure
%
% FUNCTION foldSimData
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
if nargin>=1 ; SimData.maxAlpha  =maxAlpha  ; else  SimData.maxAlpha  = []; end
if nargin>=2 ; SimData.segAlpha  =segAlpha  ; else  SimData.segAlpha  = []; end
if nargin>=3 ; SimData.dAlpha    =dAlpha    ; else  SimData.dAlpha    = []; end
if nargin>=4 ; SimData.nlvl      =nlvl      ; else  SimData.nlvl      = []; end
if nargin>=5 ; SimData.taylorN   =taylorN   ; else  SimData.taylorN   = []; end
if nargin>=6 ; SimData.alphaTol  =alphaTol  ; else  SimData.alphaTol  = []; end
if nargin>=7 ; SimData.diffTol   =diffTol   ; else  SimData.diffTol   = []; end
if nargin>=8 ; SimData.diffTolMax=diffTolMax; else  SimData.diffTolMax= []; end
if nargin>=9 ; SimData.method    =method    ; else  SimData.method    = []; end
if nargin>=10; SimData.varOpt    =varOpt    ; else  SimData.varOpt    = []; end
end