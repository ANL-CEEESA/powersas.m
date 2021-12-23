function [MatGV0,MatGV1,MatGRhs0,MatGRhs1]=getLinearInterpolatorSyn(...
    syn,synNew,d0,dNew,ed0,ed10,ed20,edNew,ed1New,ed2New,eq0,eq10,eq20,eqNew,eq1New,eq2New,psid0,psiq0,psidNew,psiqNew)
% Simplify synchronous generator as variant Thevenin model
%
% FUNCTION getLinearInterpolatorSyn
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
%	syn - synchronous generator data
%	synNew - updated synchronous generator data
%	*0 - initial states of synchronous generator
%	*New - updated states of synchronous generator
%
% OUTPUT
%	MatGV0 - equivalent admittance corresponding to initial states
%	MatGV1 - the change of equivalent admittance 
%	MatGRhs0 - equivalent voltage source corresponding to initial states
%	MatGRhs1 - the change of equivalent voltage source  
% 

nSyn=size(syn,1);
    
if nSyn>0
    synIdx=syn(:,1);
    synMod=syn(:,5);
    Rga=syn(:,7);
    Xgd=syn(:,8);
    Xgd1=syn(:,9);
    Xgd2=syn(:,10);
    Xgq=syn(:,13);
    Xgq1=syn(:,14);
    Xgq2=syn(:,15);
    
    if ~isempty(synNew)
        RgaNew=synNew(:,7);
        XgdNew=synNew(:,8);
        Xgd1New=synNew(:,9);
        Xgd2New=synNew(:,10);
        XgqNew=synNew(:,13);
        Xgq1New=synNew(:,14);
        Xgq2New=synNew(:,15);        
    else
        RgaNew=Rga;
        XgdNew=Xgd;
        Xgd1New=Xgd1;
        Xgd2New=Xgd2;
        XgqNew=Xgq;
        Xgq1New=Xgq1;
        Xgq2New=Xgq2;  
    end
    
    cosd0=cos(d0);
    sind0=sin(d0);
    cosdNew=cos(dNew);
    sindNew=sin(dNew);
    
    MatGV0=zeros(nSyn,4);
    MatGV1=zeros(nSyn,4);
    
    MatGRhs0=zeros(nSyn,2);
    MatGRhs1=zeros(nSyn,2);
    
    for i=1:nSyn
        Mg0=[sind0(i),-cosd0(i);cosd0(i),sind0(i)];
        MgNew=[sindNew(i),-cosdNew(i);cosdNew(i),sindNew(i)];
        if synMod(i)==6||synMod(i)==5
            Yg=[Rga(i),Xgq2(i);-Xgd2(i),Rga(i)]/(Rga(i)*Rga(i)+Xgd2(i)*Xgq2(i));
            YgNew=[RgaNew(i),Xgq2New(i);-Xgd2New(i),RgaNew(i)]/(RgaNew(i)*RgaNew(i)+Xgd2New(i)*Xgq2New(i));
            MatGRhs0(i,:)=(Mg0'*Yg*[ed20(i);eq20(i)]).';
            MatGRhs1(i,:)=(MgNew'*YgNew*[ed2New(i);eq2New(i)]).'-MatGRhs0(i,:);
        elseif synMod(i)==4
            Yg=[Rga(i),Xgq1(i);-Xgd1(i),Rga(i)]/(Rga(i)*Rga(i)+Xgd1(i)*Xgq1(i));
            YgNew=[RgaNew(i),Xgq1New(i);-Xgd1New(i),RgaNew(i)]/(RgaNew(i)*RgaNew(i)+Xgd1New(i)*Xgq1New(i));
            MatGRhs0(i,:)=(Mg0'*Yg*[ed10(i);eq10(i)]).';
            MatGRhs1(i,:)=(MgNew'*YgNew*[ed1New(i);eq1New(i)]).'-MatGRhs0(i,:);
        elseif synMod(i)==3
            Yg=[Rga(i),Xgq(i);-Xgd1(i),Rga(i)]/(Rga(i)*Rga(i)+Xgd1(i)*Xgq(i));
            YgNew=[RgaNew(i),XgqNew(i);-Xgd1New(i),RgaNew(i)]/(RgaNew(i)*RgaNew(i)+Xgd1New(i)*XgqNew(i));
            MatGRhs0(i,:)=(Mg0'*Yg*[ed0(i);eq10(i)]).';
            MatGRhs1(i,:)=(MgNew'*YgNew*[edNew(i);eq1New(i)]).'-MatGRhs0(i,:);
        elseif synMod(i)==2
            Yg=[Rga(i),Xgq(i);-Xgd(i),Rga(i)]/(Rga(i)*Rga(i)+Xgd(i)*Xgq(i));
            YgNew=[RgaNew(i),XgqNew(i);-XgdNew(i),RgaNew(i)]/(RgaNew(i)*RgaNew(i)+XgdNew(i)*XgqNew(i));
            MatGRhs0(i,:)=(Mg0'*Yg*[ed0(i);eq0(i)]).';
            MatGRhs1(i,:)=(MgNew'*YgNew*[edNew(i);eqNew(i)]).'-MatGRhs0(i,:);
        elseif synMod(i)==8
            Yg=zeros(2);
            YgNew=zeros(2);
            MatGRhs0(i,:)=(Mg0'*[(eq20(i)-psid0(i))/Xgd2(i);(-ed20(i)-psiq0(i))/Xgq2(i)]).';
            MatGRhs1(i,:)=(MgNew'*[(eq2New(i)-psidNew(i))/Xgd2New(i);(-ed2New(i)-psiqNew(i))/Xgq2New(i)]).'-MatGRhs0(i,:);
        end
        MatVg0=Mg0'*Yg*Mg0;
        MatVgNew=MgNew'*YgNew*MgNew;
        MatGV0(i,:)=[MatVg0(1,1),MatVg0(1,2),MatVg0(2,1),MatVg0(2,2)];
        MatGV1(i,:)=[MatVgNew(1,1),MatVgNew(1,2),MatVgNew(2,1),MatVgNew(2,2)]-MatGV0(i,:);
    end
else
    
    MatGV0=zeros(nSyn,4);
    MatGV1=zeros(nSyn,4);
    
    MatGRhs0=zeros(nSyn,2);
    MatGRhs1=zeros(nSyn,2);

end
end