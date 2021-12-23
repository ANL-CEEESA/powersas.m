function [Y,Ytr,Ysh,ytrfr,ytrto,yshfr,yshto]=getYMatrix(nbus,line,fault)
% Generate admittance matrix
%
% FUNCTION getYMatrix
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
% fault has 4 columns: 
% line number, relative position, fault resistance, fault reactance

if nargin<3
    fault=[];
end
r=line(:,8);
x=line(:,9);
z=r+1j*x;
b=line(:,10);
chrg1 = 1j*0.5*line(:,16).*b;
chrg2 = 1j*0.5*line(:,16).*b;
y = line(:,16)./z;
if ~isempty(fault)
    zf1=z(fault(:,1)).*fault(:,2);
    zf2=z(fault(:,1)).*(1-fault(:,2));
    zfg=fault(:,3)+1j*fault(:,4);
    % delta-wye transform
    zden=zf1.*zf2+zf2.*zfg+zfg.*zf1;
    yffr=zf2./zden;
    yfto=zf1./zden;
    yftr=zfg./zden;
    
    rzf2=1./zf2;
    yfto(zf1==0&zfg==0)=rzf2(zf1==0&zfg==0);
    yftr(zf1==0&zfg==0)=rzf2(zf1==0&zfg==0);
    
    rzf1=1./zf1;
    yffr(zf2==0&zfg==0)=rzf1(zf2==0&zfg==0);
    yftr(zf2==0&zfg==0)=rzf1(zf2==0&zfg==0);
    
    y(fault(:,1))=yftr;
    chrg1(fault(:,1))=chrg1(fault(:,1))+yffr;
    chrg2(fault(:,1))=chrg2(fault(:,1))+yfto;    
end

line(line(:,11)==0,11)=1;
ts = line(:,11).*exp(1j*line(:,12)*pi/180);
ts2= ts.*conj(ts);

ytrfr=y./conj(ts);
ytrto=y./ts;
yshfr=(y+chrg1)./ts2-ytrfr;
yshto=y+chrg2-ytrto;

i_fr=line(:,1);
i_to=line(:,2);

Y = sparse(i_fr,i_to,-y./conj(ts),nbus,nbus) + ...
    sparse(i_to,i_fr,-y./ts,nbus,nbus) + ...
    sparse(i_fr,i_fr,(y+chrg1)./ts2,nbus,nbus)+ ...
    sparse(i_to,i_to,y+chrg2,nbus,nbus);

Ysh=sum(Y,2);
Ytr=Y-sparse(1:nbus,1:nbus,Ysh,nbus,nbus);
end