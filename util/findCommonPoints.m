function [tt,sc1x,sc2x]=findCommonPoints(tt1,sc1,tt2,sc2,tt)
% Regulate two series with different time steps tt1/tt2 to the same time steps tt
%
% FUNCTION findCommonPoints
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
%	tt1 - Time series 1
%	sc1 - State series corresponding to time steps tt1
%	tt2 - Time series 2
%	sc2 - State series corresponding to time steps tt2
%	tt - Regulated time series
%
% OUTPUT
%	tt - Regulated time series
%	sc1x - Regulated state series 1
%	sc2x - Regulated state series 2
%
sc1x=zeros(size(sc1,1),size(tt,2));
sc2x=zeros(size(sc2,1),size(tt,2));
tth=1e-8;

for i=1:size(tt,2)
    
    idxttst=find(tt>=tt(i),1);
    idxtten=find(tt<=tt(i),1,'last');
    
    idx1st=find(tt1<=tt(i),1,'last');
    idx1en=find(tt1>=tt(i),1);
    if idx1st==idx1en
        sc1x(:,i)=sc1(:,idx1st);
    else
        if abs(tt1(idx1en)-tt1(idx1st))<tth
            if idxttst==idxtten
                sc1x(:,i)=(sc1(:,idx1st)+sc1(:,idx1en))/2;
            else
                idx1streal=min([idx1st,idx1en]);
                idx1enreal=max([idx1st,idx1en]);
                sc1x(:,i)=sc1(:,idx1streal+round((i-idxttst)*(idx1enreal-idx1streal)/(idxtten-idxttst)));
            end
        else
            sc1x(:,i)=(sc1(:,idx1st)*(tt1(idx1en)-tt(i))+sc1(:,idx1en)*(tt(i)-tt1(idx1st)))/...
                (tt1(idx1en)-tt1(idx1st));
        end
    end
    
    idx2st=find(tt2<=tt(i),1,'last');
    idx2en=find(tt2>=tt(i),1);
    if idx2st==idx2en
        sc2x(:,i)=sc2(:,idx2st);
    else
        if abs(tt2(idx2en)-tt2(idx2st))<tth
            if idxttst==idxtten
                sc2x(:,i)=(sc2(:,idx2st)+sc2(:,idx2en))/2;
            else
                idx2streal=min([idx2st,idx2en]);
                idx2enreal=max([idx2st,idx2en]);
                sc2x(:,i)=sc2(:,idx2streal+round((i-idxttst)*(idx2enreal-idx2streal)/(idxtten-idxttst)));
            end
        else
            sc2x(:,i)=(sc2(:,idx2st)*(tt2(idx2en)-tt(i))+sc2(:,idx2en)*(tt(i)-tt2(idx2st)))/...
                (tt2(idx2en)-tt2(idx2st));
        end
    end
end
end