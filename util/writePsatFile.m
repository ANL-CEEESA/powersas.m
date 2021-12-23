function writePsatFile(psatfilename,baseMVA, bus, gen, branch,area,dcline, gencost,namebus)
% Write PSAT data file from Matpower data 
%
% FUNCTION writePsatFile
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

%
% ----- write data in MATPOWER format -----
fid = fopen(psatfilename, 'w');

% print out function title, baseMVA, etc.
temp = psatfilename(1:length(psatfilename)-2);

% process gen
nbus=size(bus,1);
genAccum=zeros(nbus,size(gen,2));
for i=2:size(gen,2)
    genAccum(:,i)=accumarray(gen(:,1),gen(:,i),[nbus,1]);
end
genAccum(:,1)=(1:nbus)';
genCount=accumarray(gen(:,1),ones(size(gen,1),1),[nbus,1]);
genTag=zeros(nbus,1);
genTag(gen(:,1))=1;
genAccumx=genAccum(genTag==1,:);
genCountx=genCount(genTag==1);
genAccumx(:,6)=genAccumx(:,6)./genCountx;
genAccumx(:,8)=genAccumx(:,8)./genCountx;
genAccumx(genAccumx(:,8)>0,8)=1;
gen=genAccumx;

newToOld=bus(:,1);
oldToNew=zeros(max(newToOld),1);
oldToNew(newToOld)=1:nbus;

bus(:,1)=oldToNew(bus(:,1));
gen(:,1)=oldToNew(gen(:,1));
branch(:,1)=oldToNew(branch(:,1));
branch(:,2)=oldToNew(branch(:,2));

% print out bus matrix [OK]
fprintf(fid, '\n%%bus data');
fprintf(fid, '\nBus.con = [...');
for i = 1:size(bus, 1)
    % 		            bus  baseKv  v  theta area region
    fprintf(fid, '\n%4i %14.10g %14.10g %14.10g %4i %4i', bus(i,[1,10,8,9,7,11]));
end
fprintf(fid, '\n];\n');

% print out SW matrix [OK]
sw=bus((bus(:,2)==3),:);
swg=gen(gen(:,1)==sw(1),:);
fprintf(fid, '\n%%SW data');
fprintf(fid, '\nSW.con = [...');
% 	 bus  Sn  Vn v theta qmax qmin vmax vmin pg_guess loss_part ref u
fprintf(fid, '\n%4i %14.10g %14.10g %14.10g %14.10g %14.10g %14.10g %14.10g %14.10g %14.10g %14.10g %4i %4i', [sw(1),baseMVA,sw(10),swg([6,4,5]),sw([12,13]),0.0,1,1,swg(8)]);
fprintf(fid, '\n];\n');

% print out PV matrix [OK]
fprintf(fid, '\n%%PV data');
fprintf(fid, '\nPV.con = [...');
for i = 1:size(gen, 1)
    % 	bus  Sn  Vn Pg V qmax qmin vmax vmin loss_part u
    if gen(i,1)~=sw(1)
        fprintf(fid, '\n%4i %14.10g %14.10g %14.10g %14.10g %14.10g %14.10g %14.10g %14.10g %14.10g %4i', [gen(i,1),baseMVA,bus(gen(i,1),10),gen(i,2)/baseMVA,gen(i,6),gen(i,[4,5])/baseMVA,1.05,0.95,1,gen(i,8)]);
    end
end
fprintf(fid, '\n];\n');

% print out PQ matrix [OK]
fprintf(fid, '\n%%PQ data');
fprintf(fid, '\nPQ.con = [...');
for i = 1:size(bus, 1)
    % 	bus  Sn  Vn P Q vmax vmin z u
    if bus(i,1)~=sw(1)%&&isempty(find(gen(:,1)==bus(i,1), 1))
        fprintf(fid, '\n%4i %14.10g %14.10g %14.10g %14.10g %14.10g %14.10g %4i %4i', [bus(i,1),baseMVA,bus(i,10),bus(i,[3,4])/baseMVA,bus(i,[12,13]),1,1]);
    end
end
fprintf(fid, '\n];\n');


% print out Shunt matrix [OK]
fprintf(fid, '\n%%Shunt data');
fprintf(fid, '\nShunt.con = [...');
for i = 1:size(bus, 1)
    % 	bus  Sn  fn Gs Bs z u
    if bus(i,5)~=0||bus(i,6)~=0
        fprintf(fid, '\n%4i %14.10g %14.10g %14.10g %14.10g %4i %4i', [bus(i,1),baseMVA,bus(i,10),60.0,bus(i,[5,6])/baseMVA,1]);
    end
end
fprintf(fid, '\n];\n');

% print out Branch matrix [OK]
fprintf(fid, '\n%%Line Data');
fprintf(fid, '\nLine.con = [...');
for i = 1:size(branch, 1)
    % 	f t Sn Vn fn ln r x b k phi Im Pm Sm u
    bi=branch(i,1);
    bj=branch(i,2);
    kv=bus(bi,10);
    if branch(i,9)~=0
        bratio=bus(bi,10)/bus(bj,10);
    else
        bratio=0;
    end
    fprintf(fid, '\n%4i %4i %14.10g %14.10g %14.10g %4i %14.10g %14.10g %14.10g %14.10g %14.10g %14.10g %14.10g %14.10g %14.10g %14.10g %4i', [branch(i,[1,2]),baseMVA,kv,60.0,0.0,bratio,branch(i,[3,4,5,9,10]),branch(i,[6,6,6])/baseMVA,branch(i,11)]);
end
fprintf(fid, '\n];\n');

% print out Supply matrix [OK]
fprintf(fid, '\n%%Supply Data');
fprintf(fid, '\nSupply.con = [...');
for i = 1:size(gen, 1)
    % 	bus  Sn  * Pmax Pmin * * * * * * * * *   
    fprintf(fid, '\n%4i %14.10g %4i %14.10g %14.10g %4i %4i %14.10g %14.10g %4i %4i %4i %4i %4i %4i %14.10g %14.10g %4i', [gen(i,1),baseMVA,0,gen(i,9)/baseMVA,  0.0  , 0, 0, 9.0, 0.1, 0, 0, 0, 0, 0, 1,gen(i,[4,5])/baseMVA,gen(i,8)]);   
end
fprintf(fid, '\n];\n');

% print out hvdc matrix
% if ~isempty(dcline)
fprintf(fid, '\n%%fbus tbus Sn Vrn Vin fn V I Xtr Xti mr mi Ki Kp Rdc Ldc aMax aMin gaMax gaMin yRmax yRmin yiMax yiMin ctype io po vo u');
fprintf(fid, '\nHvdc.con = [...');
for i = 1:size(dcline, 1)
    fprintf(fid, '\n%4i   %4i  %14.10g %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g  %14.10g %4i',...
        [dcline(i, [1,2]),baseMVA,dcline(i,[8,9])*baseMVA,60.0,dcline(i,8)*baseMVA,dcline(i,4)/(dcline(i,8)*baseMVA),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,dcline(i,[11,10])/(dcline(i,8)*baseMVA),dcline(i,[11,10])/(dcline(i,8)*baseMVA),1,1,1,1,1]);
end;
fprintf(fid, '\n];	\n');
% end

% bus name
fprintf(fid, 'Bus.names = {...\n      ');
for i = 1:length(namebus(:,1))-1
  tmpName = namebus{i};
  newName = strrep(tmpName,'''',' ');%loadcase
  fprintf(fid, '''%s '';',newName);
  if rem(i,5) == 0; fprintf(fid,'\n      '); end
end
fprintf(fid, '''%s ''};\n',namebus{end});

fclose(fid);

end