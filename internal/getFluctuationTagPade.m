function flucTag=getFluctuationTagPade(coeffA,coeffB,maxTime,rateThreshold)
% Determine steady state from the fluctuation bounds using Pade approximation
%
% FUNCTION getFluctuationTagPade
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
%	coeffA - Pade coefficients on the numerator
%	coeffB - Pade coefficients on the denomenator
%	maxTime - maximum time considered
%	rateThreshold - the threshold of fluctuation rate
%
% OUTPUT
%	flucTag
%		1 - lower than threshold, entered steady state
%		0 - Not entered steady state yet
%

if nargin<4
    rateThreshold=0.0005;
end

NA=size(coeffA,2);
timePwrA=maxTime.^(1:(NA-1));
coeffAx=coeffA(:,2:NA);
constA=coeffA(:,1);
fluctuationPosA=sum(coeffAx.*(coeffAx>0).*repmat(timePwrA,size(coeffA,1),1),2);
fluctuationNegA=sum(-coeffAx.*(coeffAx<0).*repmat(timePwrA,size(coeffA,1),1),2);
NB=size(coeffB,2);
timePwrB=maxTime.^(1:(NB-1));
coeffBx=coeffB(:,2:NB);
constB=coeffB(:,1);
constTotal=constA./constB;
fluctuationPosB=sum(coeffBx.*(coeffBx>0).*repmat(timePwrB,size(coeffB,1),1),2);
fluctuationNegB=sum(-coeffBx.*(coeffBx<0).*repmat(timePwrB,size(coeffB,1),1),2);

denSgn=(constB+fluctuationPosB).*(constB-fluctuationNegB); % if the denomenator will change sign, then fluctuation is too big.

totalFluctuation=max(abs([(constA+fluctuationPosA)./(constB+fluctuationPosB)-constTotal,...
    (constA+fluctuationPosA)./(constB-fluctuationNegB)-constTotal,...
    (constA-fluctuationNegA)./(constB+fluctuationPosB)-constTotal,...
    (constA-fluctuationNegA)./(constB-fluctuationNegB)-constTotal]),[],2);

% Special: check if the coefficients on Num/Den can be directly eliminated
totalFlucRes=zeros(size(coeffA,1),1);
if NA>NB
    regA=coeffA(:,1:NB);
    res=coeffA(:,(NB+1):end);
    regB=coeffB;
    flucResPos=sum(res.*(res>0).*repmat(timePwrA(NB:end),size(res,1),1),2);
    flucResNeg=sum(-res.*(res<0).*repmat(timePwrA(NB:end),size(res,1),1),2);
    flucResNum=max([flucResPos,flucResNeg],[],2);
    flucRes=abs([flucResNum./(constB-fluctuationNegB),flucResNum./(constB+fluctuationPosB)]);
    flucRes(all(regB(:,2:end)>=0,2),1)=0;
    flucRes(all(regB(:,2:end)<=0,2),2)=0;
    totalFlucRes=max(flucRes,[],2);
elseif NB>NA
    regA=coeffA;
    regB=coeffB(:,1:NA);
    res=coeffB(:,(NA+1):end);
    flucResPos=sum(res.*(res>0).*repmat(timePwrB(NA:end),size(res,1),1),2);
    flucResNeg=sum(-res.*(res<0).*repmat(timePwrB(NA:end),size(res,1),1),2);
    flucResNum=max([flucResPos,flucResNeg],[],2);
    flucRes=abs([flucResNum./(constB-fluctuationNegB),flucResNum./(constB+fluctuationPosB)]);
    flucRes(all(regB(:,2:end)>=0,2),1)=0;
    flucRes(all(regB(:,2:end)<=0,2),2)=0;
    totalFlucRes=max(flucRes,[],2);
else
    regA=coeffA;
    regB=coeffB;
end

ratioAB=regA./regB;
ratioAB(isnan(ratioAB))=0;
ratioAB=ratioAB./repmat(mean(ratioAB,2),1,size(ratioAB,2));
flucRatioAB=max(ratioAB,[],2)-min(ratioAB,[],2);

eachFlucTag1=(totalFluctuation/maxTime)<rateThreshold;
eachFlucTag2=(flucRatioAB<rateThreshold)&((totalFlucRes/maxTime)<rateThreshold);
eachFlucTag=eachFlucTag1|eachFlucTag2;
flucTag=all(eachFlucTag);
end