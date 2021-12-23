function x=solveToepLevinson(ct,y)
% Solving Toeplitz matrix equations with Levinson algorithm (vector)
%
% FUNCTION solveToepLevinson
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
%	ct - D*(2N-1)
%	y - D*N
%
% OUTPUT
%	x - solution to the equations
%

% overheadTag=tic;
D=size(ct,1); % The dimension of the variables
N=round((size(ct,2)+1)/2); % The size of the Toep matrix

f=zeros(D,N);
b=zeros(D,N);
temp=zeros(D,N);
x=zeros(D,N);
epsf=zeros(D,1);
epsb=zeros(D,1);
epsx=zeros(D,1);
alphaf=zeros(D,1);
betaf=zeros(D,1);
alphab=zeros(D,1);
betab=zeros(D,1);

% oh=toc(overheadTag);
% disp(sprintf('Overhead=%10.8f s.',oh));

% mainTag=tic;
f(:,1)=1./ct(:,N);
b(:,end)=f(:,1);
x(:,1)=y(:,1).*f(:,1);
for k=2:N
    epsf(:)=sum(f(:,1:(k-1)).*ct(:,((N+k-1):-1:(N+1))),2);
    epsb(:)=sum(b(:,(N-k+2):N).*ct(:,((N-1):-1:(N-k+1))),2);
    epsx(:)=sum(x(:,1:(k-1)).*ct(:,((N+k-1):-1:(N+1))),2);
    
    alphaf(:)=1./(1-epsf.*epsb);
    betaf(:)=-epsf.*alphaf;
    alphab(:)=-epsb.*alphaf;
    betab(:)=alphaf;
    
    temp(:,1:(k-1))=f(:,1:(k-1));
    f(:,1:k)=repmat(alphaf,1,k).*temp(:,1:k)+repmat(betaf,1,k).*b(:,(N-k+1):N);
    b(:,(N-k+1):N)=repmat(alphab,1,k).*temp(:,1:k)+repmat(betab,1,k).*b(:,(N-k+1):N);
    
    x(:,1:k)=x(:,1:k)+repmat((y(:,k)-epsx),1,k).*b(:,(N-k+1):N);
end

% mn=toc(mainTag);
% disp(sprintf('Main=%10.8f s.',mn));
% disp(sprintf('xxxxx/=%10.8f ',oh/(oh+mn)));
end