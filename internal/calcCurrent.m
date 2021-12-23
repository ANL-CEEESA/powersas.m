function [Ig,Im,Ii,Is]=calcCurrent(xt,syn,ind,zip,SVec,vIdx,sIdx,deltaIdx,omegaIdx,eq1Idx,eq2Idx,ed1Idx,ed2Idx,psidIdx,psiqIdx,efIdx)
% [INTERNAL] Calculate some injection currents
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
% FUNCTION calcCurrent
% CAUTION: this function is used internally. Do not call it directly unless you know what you are doing.
%
% INPUT
%
% OUTPUT
%	Ig - Currents of syn generators
%	Im - Current of induction motors
%	Ii - Current of the I component of ZIP loads
%	Is - Current of the PQ loads from SVec 
%

dg=xt(deltaIdx);
wg=xt(omegaIdx);
eq1=xt(eq1Idx);
eq2=xt(eq2Idx);
ed1=xt(ed1Idx);
ed2=xt(ed2Idx);
psid=xt(psidIdx);
psiq=xt(psiqIdx);
vbus=xt(vIdx);
s=xt(sIdx);
Ef=xt(efIdx);

if ~isempty(syn)
    nSyn=size(syn,1);
    synIdx=syn(:,1);
    Vgx=real(vbus(synIdx));
    Vgy=imag(vbus(synIdx));
    Vgd=Vgx.*sin(dg)-Vgy.*cos(dg);
    Vgq=Vgx.*cos(dg)+Vgy.*sin(dg);
    wgb=syn(:,4);
    model=syn(:,5);
    Xgl=syn(:,6);
    Rga=syn(:,7);
    Xgd=syn(:,8);
    Xgd1=syn(:,9);
    Xgd2=syn(:,10);
    Tgd1=syn(:,11);
    Tgd2=syn(:,12);
    Xgq=syn(:,13);
    Xgq1=syn(:,14);
    Xgq2=syn(:,15);
    Tgq1=syn(:,16);
    Tgq2=syn(:,17);
    Mg=syn(:,18);
    Dg=syn(:,19);
    TgAA=syn(:,24);
    gammad=Tgd2./Tgd1.*Xgd2./Xgd1.*(Xgd-Xgd1);
    gammaq=Tgq2./Tgq1.*Xgq2./Xgq1.*(Xgq-Xgq1);
    
    Igd=zeros(nSyn,1);
    Igq=zeros(nSyn,1);
end

if ~isempty(ind)
    nInd=size(ind,1);
    indIdx=ind(:,1);
    Zm1=ind(:,7)+1j*ind(:,8);
    Rm2=ind(:,9);
    Xm2=ind(:,10);
    Zmm=1j*ind(:,13);
    Tm0=ind(:,15)+ind(:,16)+ind(:,17);
    Tm1=-ind(:,16)-2*ind(:,17);
    Tm2=ind(:,17);
    Hm=ind(:,14);
end

if ~isempty(syn)
    for i=1:nSyn
        if model(i)==8
            Igd(i)=(eq2(i)-psid(i))/Xgd2(i);
            Igq(i)=(-ed2(i)-psiq(i))/Xgq2(i);
        elseif model(i)==6
            Igd(i)=(Rga(i)*(ed2(i)-Vgd(i))+Xgq2(i)*(eq2(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd2(i)*Xgq2(i));
            Igq(i)=(-Xgd2(i)*(ed2(i)-Vgd(i))+Rga(i)*(eq2(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd2(i)*Xgq2(i));
        elseif model(i)==5
            Igd(i)=(Rga(i)*(ed2(i)-Vgd(i))+Xgq2(i)*(eq2(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd2(i)*Xgq2(i));
            Igq(i)=(-Xgd2(i)*(ed2(i)-Vgd(i))+Rga(i)*(eq2(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd2(i)*Xgq2(i));
        elseif model(i)==4
            Igd(i)=(Rga(i)*(ed1(i)-Vgd(i))+Xgq1(i)*(eq1(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd1(i)*Xgq1(i));
            Igq(i)=(-Xgd1(i)*(ed1(i)-Vgd(i))+Rga(i)*(eq1(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd1(i)*Xgq1(i));
        elseif model(i)==3
            Igd(i)=(Rga(i)*(-Vgd(i))+Xgq(i)*(eq1(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd1(i)*Xgq(i));
            Igq(i)=(-Xgd1(i)*(-Vgd(i))+Rga(i)*(eq1(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd1(i)*Xgq(i));
        elseif model(i)==2
            Igd(i)=(Rga(i)*(-Vgd(i))+Xgq(i)*(Ef(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd1(i)*Xgq(i));
            Igq(i)=(-Xgd(i)*(-Vgd(i))+Rga(i)*(Ef(i)-Vgq(i)))/(Rga(i)*Rga(i)+Xgd(i)*Xgq(i));
        end
    end
    
    Ig=(Igd.*sin(dg)+Igq.*cos(dg))+1j*(-Igd.*cos(dg)+Igq.*sin(dg));
else
    Ig=[];
end

if ~isempty(ind)
    yind=1./(Zm1+Zmm.*(Rm2+1j*Xm2.*s)./(Rm2+s.*(Zmm+1j*Xm2)));
    Im=vbus(indIdx).*yind;
else
    Im=[];
end

if ~isempty(zip)
    Ii=(zip(:,6)-1j*zip(:,9)).*vbus(zip(:,1))./abs(vbus(zip(:,1)));
else
    Ii=[];
end

Is=conj(SVec./vbus);
end
