function [YshInd0,Yshind1]=getLinearInterpolatorInd(nbus,ind,s0,sNew)
% Simplify induction motor as variant admittance model
%
% FUNCTION getLinearInterpolatorInd
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
%	nbus - number of buses in the system
%	ind - induction motor data
%	s0 - initial slip
%	sNew - updated slip
%
% OUTPUT
%	YshInd0 - equivalent admittance corresponding to slip s0
%	Yshind1 - the change of equivalent admittance from slip s0 to sNew
% 

nInd=size(ind,1);
indIdx=ind(:,1);

Z1=ind(:,7)+1j*ind(:,8);
Ze=1j*ind(:,13);
R2=ind(:,9);
X2=ind(:,10);

Zind0=Z1+Ze.*(R2+1j*X2.*s0)./(R2+(Ze+1j*X2).*s0);
Yind0=1./Zind0;

ZindNew=Z1+Ze.*(R2+1j*X2.*sNew)./(R2+(Ze+1j*X2).*sNew);
YindNew=1./ZindNew;

YshInd0=accumarray(indIdx,Yind0,[nbus,1]);
Yshind1=accumarray(indIdx,YindNew-Yind0,[nbus,1]);
end