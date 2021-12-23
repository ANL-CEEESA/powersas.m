function [P,Q,Jac]=getJacobianPolar(V0,Y)
% Calculating Jacobian matrix in polar coordinate
% FUNCTION getJacobianPolar
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
% INPUT
%	V0 - voltage
%	Y - admittance matrix
% 
% OUTPUT
%	P - Net active power injection
%	Q - Net reactive power injection
%	Jac - Jacobian matrix
% 
nbus=size(V0,1);
Vm=abs(V0);
Va=angle(V0);
Vadiff=repmat(Va,1,nbus)-repmat(Va',nbus,1);
G=real(Y);B=imag(Y);
Gd=diag(G);Bd=diag(B);

cVa=cos(Vadiff);sVa=sin(Vadiff);
P=Vm.*((G.*cVa+B.*sVa)*Vm);
Q=-Vm.*((B.*cVa-G.*sVa)*Vm);

Vmult=repmat(Vm,1,nbus).*repmat(Vm',nbus,1);
H=-Vmult.*(G.*sVa-B.*cVa);
H(1:nbus+1:end)=Q+Vm.*Vm.*Bd;
N=-repmat(Vm,1,nbus).*(G.*cVa+B.*sVa);
N(1:nbus+1:end)=-P./Vm-Vm.*Gd;
J=Vmult.*(G.*cVa+B.*sVa);
J(1:nbus+1:end)=-P+Vm.*Vm.*Gd;
L=repmat(Vm,1,nbus).*(G.*sVa-B.*cVa);
L(1:nbus+1:end)=-Q./Vm+Vm.*Bd;

Jac=[H,N;J,L];
end