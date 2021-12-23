function [nIslands,islands,refs]=searchIslands(nodes,links)
% Search all islands in a system
%
% FUNCTION searchIslands
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
%	nodes - bus numbers
%	linkes - from and to terminals of all the branches
%
% OUTPUT
%	nIslands - number of islands
%	islands - the index of the island that the bus is in
%	refs - the indexes of reference buses
%

% islandCode=nodes;%assume that nodes are listed as 1:nBus.
% 
% nLine=size(links,1);
% for i=1:nLine
%     from=links(i,1);
%     to=links(i,2);
%     
%     fromIsland=islandCode(from);
%     toIsland=islandCode(to);
%     
%     if fromIsland<toIsland
%         islandCode(islandCode==toIsland)=fromIsland;
%     elseif fromIsland>toIsland
%         islandCode(islandCode==fromIsland)=toIsland;
%     end
% end
% 
% cntIsland=islandCode-nodes;
% nIslands=sum(cntIsland==0);
% islands=islandCode;
% islandBusNumber=find(cntIsland==0);
% 
% for i=1:nIslands
%     islands(islandCode==islandBusNumber(i))=i;
% end
% [sislands,iislands]=sort(islands);
% if size(islands,1)>size(islands,2)
%     refTag=sislands'-[0,sislands(1:(length(sislands)-1))'];
% else
%     refTag=sislands-[0,sislands(1:(length(sislands)-1))];
% end
% refs=iislands(refTag~=0);

global IS_OCTAVE;

nNodes=max([links(:,1);links(:,2);nodes]);
A=sparse([links(:,1);links(:,2)],[links(:,2);links(:,1)],[ones(size(links(:,1)));ones(size(links(:,1)))],nNodes,nNodes);

if IS_OCTAVE
    islands=conncompx(A);
else
    G=graph(A);
    [islands] = conncomp(G);
end
nIslands=max(islands);
if size(islands,1)<size(islands,2)
    islands=islands';
end
[sislandsx,iislandsx]=sort(islands);
refTagx=sislandsx-[0;sislandsx(1:(length(sislandsx)-1))];
refs=iislandsx(refTagx~=0);

end

function islands=conncompx(adj)
  
  nNodes=size(adj,1);
  
  visited=zeros(nNodes,1);
  islands=zeros(nNodes,1);
  
  is=1;
  while 1
    
    s=find(visited==0,1);
    if isempty(s)
      break;      
    end
    
    nb=s;
    nbt=nb;
    
    while 1
      sz=length(nb);
      szt=length(nbt);    
      
      nbx=mod(find(adj(:,nb(1:sz))~=0)-1,nNodes)+1;
%       nbx=zeros(0,1);
%       for i=1:sz
%         nbx=[nbx;find(adj(:,nb(i))~=0)];      
%       end
      nb=unique(nbx);
      nbt=unique([nbt;nb]);
      if length(nbt)==szt
        break;
      end  
    end
    islands(nbt)=is;
    visited(nbt)=1;
    
    is=is+1;
  end
  
end