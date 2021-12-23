function [fluctRate,flucTag]=getFluctuationTagPadeMod(coeffA,coeffB,maxTime,rateThreshold)
%
% FUNCTION getFluctuationTagPadeMod
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

if nargin<4
    rateThreshold=0.0005;
end

sz=max([size(coeffA,2);size(coeffB,2)]);
coeffAe=zeros(size(coeffA,1),sz);
coeffBe=zeros(size(coeffB,1),sz);
coeffAe(:,1:size(coeffA,2))=coeffA;
coeffBe(:,1:size(coeffB,2))=coeffB;

constAB=coeffAe(:,1)./coeffBe(:,1);
coeffCe=coeffAe-coeffBe.*repmat(constAB,1,sz);

[ubb,lbb]=findPolynomialBounds(coeffBe,maxTime);
[ubc,lbc]=findPolynomialBounds(coeffCe(:,2:end),maxTime);

denSgn=ubb.*lbb;

flucNumerator=max(abs([ubc,lbc]),[],2);
den=min(abs([ubb,lbb]),[],2);

fluctRate=flucNumerator./den;

flucTag=all((denSgn>0)&(fluctRate<rateThreshold));
end