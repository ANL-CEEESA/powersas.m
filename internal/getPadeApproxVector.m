function [xa,xb]=getPadeApproxVector(xc,m,tag)
% Calculate Pade coefficients by a vector of variables, called by getPadeApproxHelm()
%
% FUNCTION getPadeApproxVector
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
%	xc - HE coefficients in power series
%	m  - order on the numerator
%	tag
%		0 - proceed for computation
%		1 - this variable is ignored, do not compute
%
% OUTPUT
%	xa - Pade coeffcients on the numerator
%	xb - Pade coeffcients on the denominator 
% 


N=size(xc,2);
n=N-m;
nD=size(xc,1);

xa=zeros(nD,m);
xb=zeros(nD,n);

txc=xc(tag==0,:);
if size(txc,1)>0
    augTxc=[zeros(size(txc,1),n),txc];
    
    nT=2*n-1;
    ttxc=augTxc(:,(end-nT):(end-1));
    ytxc=-augTxc(:,(end-n+1):end);
    
    tb=solveToepLevinson(ttxc,ytxc);
    
    ta=augTxc(:,(n+1):(n+m));
    
    for k=1:n
        ta=ta+augTxc(:,(n+1-k):(n+m-k)).*repmat(tb(:,k),1,m);
    end
    
    xa(tag==0,:)=ta;
    xb(tag==0,:)=tb;
end
end