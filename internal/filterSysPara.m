function [SysParaUpd]=filterSysPara(SysPara,busFtr,swFtr,pvFtr,pqFtr,shuntFtr,lineFtr,indFtr,zipFtr,synFtr,excFtr,tgFtr)
% Filter system parameters using the filters generated by filterSysData()
%
% FUNCTION filterSysPara
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
%	SysPara - System parameters
%	*Ftr - Filter of current data to rip off inactive components
%
% OUTPUT
%	SysParaUpd - Updated system parameters
%

[pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,...
    MatGV0,MatGV1,MatGRhs0,MatGRhs1,Tmech1,Varref1,Ef1,Pm1,Eq11,fault]=unfoldSysPara(SysPara);

if ~isempty(pqIncr);pqIncr=pqIncr(pqFtr,:);end
if ~isempty(pvIncr);pvIncr=pvIncr(pvFtr,:);end
if ~isempty(Rind0);Rind0=Rind0(indFtr);end
if ~isempty(Rind1);Rind1=Rind1(indFtr);end
if ~isempty(Reind0);Reind0=Reind0(indFtr);end
if ~isempty(Reind1);Reind1=Reind1(indFtr);end
if ~isempty(Rzip0);Rzip0=Rzip0(zipFtr);end
if ~isempty(Rzip1);Rzip1=Rzip1(zipFtr);end
if ~isempty(Ytr0);Ytr0=Ytr0(busFtr,busFtr);end
if ~isempty(Ytr1);Ytr1=Ytr1(busFtr,busFtr);end
if ~isempty(Ysh0);Ysh0=Ysh0(busFtr);end
if ~isempty(Ysh1);Ysh1=Ysh1(busFtr);end
if ~isempty(VspSq2);VspSq2=VspSq2(busFtr,:);end
if ~isempty(MatGV0);MatGV0=MatGV0(synFtr,:);end
if ~isempty(MatGV1);MatGV1=MatGV1(synFtr,:);end
if ~isempty(MatGRhs0);MatGRhs0=MatGRhs0(synFtr,:);end
if ~isempty(MatGRhs1);MatGRhs1=MatGRhs1(synFtr,:);end
if ~isempty(Tmech1);Tmech1=Tmech1(tgFtr,:);end
if ~isempty(Varref1);Varref1=Varref1(excFtr,:);end
if ~isempty(Ef1);Ef1=Ef1(synFtr);end
if ~isempty(Pm1);Pm1=Pm1(synFtr);end
if ~isempty(Eq11);Eq11=Eq11(synFtr);end

if ~isempty(fault)
    nLineFtr=max([lineFtr;fault(:,1)]);
    lineTag=zeros(nLineFtr,1);
    lineTag(lineFtr)=1:length(lineFtr);
    faultTag=lineTag(fault(:,1));
    fault=fault(faultTag~=0);
    fault(:,1)=faultTag(faultTag~=0);
end

SysParaUpd=foldSysPara(pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr0,Ytr1,Ysh0,Ysh1,VspSq2,...
    MatGV0,MatGV1,MatGRhs0,MatGRhs1,Tmech1,Varref1,Ef1,Pm1,Eq11,fault);
end