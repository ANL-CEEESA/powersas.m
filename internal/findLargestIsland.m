function maxIsland=findLargestIsland(nIslands,islands,bus,pv,pq,sw)
% Find the biggest island 
%
% FUNCTION findLargestIsland
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
%	nIslands - The number of islands
%	islands - The indexes of islands
%
% OUTPUT
%	maxIsland - The index of the lagest island
%

collectPQ=zeros(size(bus,1),9);
if ~isempty(pq)
    collectPQ(pq(:,1),:)= pq;
end
% collectIsland=islands;
loads=collectPQ(:,4).*collectPQ(:,9);
islandScaleIdx=zeros(nIslands,1);
genAvailableIslands=zeros(nIslands,1);
isIsland=zeros(nIslands,1);

if ~isempty(pv)
    genAvailableIslands(islands(pv(pv(:,11)==1,1)))=1;
end
genAvailableIslands(islands(sw(:,1)))=1;

for i=1:nIslands
    isIsland(i)=sum(islands==i)-1;
    islandScaleIdx(i)=sum(loads(islands==i))*sum(islands==i);
end

islandScaleIdx=islandScaleIdx.*genAvailableIslands.*(isIsland~=0);
[maxLoad,maxIsland]=max(islandScaleIdx);
end