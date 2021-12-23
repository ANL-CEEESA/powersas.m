function SysPara=foldSysPara(pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,...
    MatGV0,MatGV1,MatGRhs0,MatGRhs1,Tmech1,Varref1,Ef1,Pm1,Eq11,fault)
% Fold system parameter settings to SysPara structure
%
% FUNCTION foldSysPara
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
if nargin>=1 ; SysPara.pqIncr  =pqIncr  ; else SysPara.pqIncr  =[];end
if nargin>=2 ; SysPara.pvIncr  =pvIncr  ; else SysPara.pvIncr  =[];end
if nargin>=3 ; SysPara.Rind0   =Rind0   ; else SysPara.Rind0   =[];end
if nargin>=4 ; SysPara.Rind1   =Rind1   ; else SysPara.Rind1   =[];end
if nargin>=5 ; SysPara.Reind0  =Reind0  ; else SysPara.Reind0  =[];end
if nargin>=6 ; SysPara.Reind1  =Reind1  ; else SysPara.Reind1  =[];end
if nargin>=7 ; SysPara.Rzip0   =Rzip0   ; else SysPara.Rzip0   =[];end
if nargin>=8 ; SysPara.Rzip1   =Rzip1   ; else SysPara.Rzip1   =[];end
if nargin>=9 ; SysPara.Ytr0    =Ytr0    ; else SysPara.Ytr0    =[];end
if nargin>=10; SysPara.Ytr1    =Ytr1    ; else SysPara.Ytr1    =[];end
if nargin>=11; SysPara.Ysh0    =Ysh0    ; else SysPara.Ysh0    =[];end
if nargin>=12; SysPara.Ysh1    =Ysh1    ; else SysPara.Ysh1    =[];end
if nargin>=13; SysPara.VspSq2  =VspSq2  ; else SysPara.VspSq2  =[];end
if nargin>=14; SysPara.MatGV0  =MatGV0  ; else SysPara.MatGV0  =[];end
if nargin>=15; SysPara.MatGV1  =MatGV1  ; else SysPara.MatGV1  =[];end
if nargin>=16; SysPara.MatGRhs0=MatGRhs0; else SysPara.MatGRhs0=[];end
if nargin>=17; SysPara.MatGRhs1=MatGRhs1; else SysPara.MatGRhs1=[];end
if nargin>=18; SysPara.Tmech1  =Tmech1  ; else SysPara.Tmech1  =[];end
if nargin>=19; SysPara.Varref1 =Varref1 ; else SysPara.Varref1 =[];end
if nargin>=20; SysPara.Ef1     =Ef1     ; else SysPara.Ef1     =[];end
if nargin>=21; SysPara.Pm1     =Pm1     ; else SysPara.Pm1     =[];end
if nargin>=22; SysPara.Eq11    =Eq11    ; else SysPara.Eq11    =[];end
if nargin>=23; SysPara.fault   =fault   ; else SysPara.fault   =[];end
end