function [cosp,sinp,n]=getTaylorPolynomials(d0,n)
% Use Taylor expansion to approximate trigonometric functions
%
% FUNCTION getTaylorPolynomials
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

if n>4
    n=4;
end
if n<0
    n=0;
end

cosd0=cos(d0);
sind0=sin(d0);

cosp=zeros(size(d0,1),n+1);
sinp=zeros(size(d0,1),n+1);

if n>=0
    cosp(:,1)=cosd0;
    sinp(:,1)=sind0;
end
if n>=1
    cosp(:,1:2)=cosp(:,1:2)+[d0.*sind0,-sind0];
    sinp(:,1:2)=sinp(:,1:2)+[-d0.*cosd0,cosd0];    
end
if n>=2
    cosp(:,1:3)=cosp(:,1:3)+[-cosd0.*d0.*d0,2*cosd0.*d0,-cosd0]/2;
    sinp(:,1:3)=sinp(:,1:3)+[-sind0.*d0.*d0,2*sind0.*d0,-sind0]/2;    
end
if n>=3
    cosp(:,1:4)=cosp(:,1:4)+[-sind0.*d0.*d0.*d0,3*sind0.*d0.*d0,-3*sind0.*d0,sind0]/6;
    sinp(:,1:4)=sinp(:,1:4)+[cosd0.*d0.*d0.*d0,-3*cosd0.*d0.*d0,3*cosd0.*d0,-cosd0]/6;
end
if n>=4
    cosp(:,1:5)=cosp(:,1:5)+[cosd0.*d0.*d0.*d0.*d0,-4*cosd0.*d0.*d0.*d0,6*cosd0.*d0.*d0,-4*cosd0.*d0,cosd0]/24;
    sinp(:,1:5)=sinp(:,1:5)+[sind0.*d0.*d0.*d0.*d0,-4*sind0.*d0.*d0.*d0,6*sind0.*d0.*d0,-4*sind0.*d0,sind0]/24;
end
end