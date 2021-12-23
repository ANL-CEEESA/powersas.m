function [pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,...
    MatGV0,MatGV1,MatGRhs0,MatGRhs1,Tmech1,Varref1,Ef1,Pm1,Eq11,fault]=unfoldSysPara(SysPara)
% function for unfolding system parameters
%
% FUNCTION unfoldSysPara
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

pqIncr=SysPara.pqIncr;
pvIncr=SysPara.pvIncr;
Rind0=SysPara.Rind0;
Rind1=SysPara.Rind1;
Reind0=SysPara.Reind0;
Reind1=SysPara.Reind1;
Rzip0=SysPara.Rzip0;
Rzip1=SysPara.Rzip1;
Ytr0=SysPara.Ytr0;
Ytr1=SysPara.Ytr1;
Ysh0=SysPara.Ysh0;
Ysh1=SysPara.Ysh1;
VspSq2=SysPara.VspSq2;

MatGV0=SysPara.MatGV0;
MatGV1=SysPara.MatGV1;
MatGRhs0=SysPara.MatGRhs0;
MatGRhs1=SysPara.MatGRhs1;

if isfield(SysPara,'Tmech1');Tmech1=SysPara.Tmech1;else Tmech1=[];end
if isfield(SysPara,'Varref1');Varref1=SysPara.Varref1;else Varref1=[];end
if isfield(SysPara,'Ef1');Ef1=SysPara.Ef1;else Ef1=[];end
if isfield(SysPara,'Pm1');Pm1=SysPara.Pm1;else Pm1=[];end
if isfield(SysPara,'Eq11');Eq11=SysPara.Eq11;else Eq11=[];end
if isfield(SysPara,'fault');fault=SysPara.fault;else fault=[];end
end