function writePsatFilefromPsat(psatfilename,bus,sw,pv,pq,shunt,line,sup,busnames,actions,zip,ind,syn,exc,tg)
% function for writing PSAT format data
%
% FUNCTION yr_writePsatFilefromPsat
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

% print out bus matrix [OK]
fprintf(fid, '\n%%bus data: %d',size(bus, 1));
fprintf(fid, '\nBus.con = [...');
for i = 1:size(bus, 1)
    % 		            bus  baseKv  v  theta area region
    fprintf(fid, '\n%6i %10.4f %10.4f %10.4f %6i %6i', bus(i,:));
end
fprintf(fid, '\n];\n');

% print out SW matrix [OK]
fprintf(fid, '\n%%SW data: %d',size(sw,1));
fprintf(fid, '\nSW.con = [...');
for i = 1:size(sw,1)
    % 	 bus  Sn  Vn v theta qmax qmin vmax vmin pg_guess loss_part ref u
    fprintf(fid, '\n%6i %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %6i %6i', sw(i,:));
end
fprintf(fid, '\n];\n');

% print out PV matrix [OK]
fprintf(fid, '\n%%PV data: %d',size(pv, 1));
fprintf(fid, '\nPV.con = [...');
for i = 1:size(pv, 1)
    % 	bus  Sn  Vn Pg V qmax qmin vmax vmin loss_part u
    fprintf(fid, '\n%6i %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %6i', pv(i,:));
    
end
fprintf(fid, '\n];\n');

% print out PQ matrix [OK]
fprintf(fid, '\n%%PQ data: %d',size(pq, 1));
fprintf(fid, '\nPQ.con = [...');
for i = 1:size(pq, 1)
    % 	bus  Sn  Vn P Q vmax vmin z u
    fprintf(fid, '\n%6i %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %6i %6i', pq(i,:));
end
fprintf(fid, '\n];\n');


% print out Shunt matrix [OK]
fprintf(fid, '\n%%Shunt data: %d',size(shunt, 1));
fprintf(fid, '\nShunt.con = [...');
for i = 1:size(shunt, 1)
    % 	bus  Sn  fn Gs Bs z u
    fprintf(fid, '\n%6i %10.4f %10.4f %10.4f %10.4f %10.4f %6i', shunt(i,:));
end
fprintf(fid, '\n];\n');

% print out Branch matrix [OK]
fprintf(fid, '\n%%Line Data: %d',size(line, 1));
fprintf(fid, '\nLine.con = [...');
for i = 1:size(line, 1)
    % 	f t Sn Vn fn ln r x b k phi Im Pm Sm u
    fprintf(fid, '\n%6d %6d %9.5f %9.5f %9.5f %6d %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %6d', line(i,:));
end
fprintf(fid, '\n];\n');

% print out Supply matrix [OK]
if ~isempty(sup)
    fprintf(fid, '\n%%Supply Data: %d',size(sup, 1));
    fprintf(fid, '\nSupply.con = [...');
    for i = 1:size(sup, 1)
        % 	bus  Sn  * Pmax Pmin * * * * * * * * *
        fprintf(fid, '\n%6d %10.4f %10.4f %10.4f %10.4f %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %10.4f %10.4f %6d %6d %6d', sup(i,:));
    end
    fprintf(fid, '\n];\n');
end

fprintf(fid, '\nHvdc.con = [];\n\n');

% print out hvdc matrix
% if ~isempty(dcline)
% fprintf(fid, '\n%%fbus tbus Sn Vrn Vin fn V I Xtr Xti mr mi Ki Kp Rdc Ldc aMax aMin gaMax gaMin yRmax yRmin yiMax yiMin ctype io po vo u');
% fprintf(fid, '\nHvdc.con = [...');
% for i = 1:size(dcline, 1)
%     fprintf(fid, '\n%6d   %6d  %8.5f %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f %6d',...
%         [dcline(i, [1,2]),baseMVA,dcline(i,[8,9])*baseMVA,60.0,dcline(i,8)*baseMVA,dcline(i,4)/(dcline(i,8)*baseMVA),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,dcline(i,[11,10])/(dcline(i,8)*baseMVA),dcline(i,[11,10])/(dcline(i,8)*baseMVA),1,1,1,1,1]);
% end;
% fprintf(fid, '\n];	\n');
% end

% bus name
fprintf(fid, 'Bus.names = {...\n      ');
if ~isempty(busnames)
for i = 1:length(busnames(:,1))-1
  tmpName = busnames{i};
  newName = strrep(tmpName,'''',' ');%loadcase
  fprintf(fid, '''%s '';',newName);
  if rem(i,5) == 0; fprintf(fid,'\n      '); end
end
fprintf(fid, '''%s ''',busnames{end});
end
fprintf(fid, '};\n');

if nargin>=14&&~isempty(actions)
    fprintf(fid, '\n%%Manual actions for testing');
    fprintf(fid, '\naction = [...');
    for i = 1:size(actions, 1)
        % 	bus  Sn  * Pmax Pmin * * * * * * * * *
        fprintf(fid, '\n%6d %6d %10.4f %6d', actions(i,:));
    end    
    fprintf(fid, '\n];\n');
end

if nargin>=15&&~isempty(zip)
    fprintf(fid, '\n%%ZIP Data: %d',size(zip, 1));
    fprintf(fid, '\nPl.con = [...');
    for i = 1:size(zip, 1)
        fprintf(fid, '\n%6d %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %2u %2u', zip(i,:));
    end
    fprintf(fid, '\n];\n'); 
else
    fprintf(fid, '\nPl.con = [];\n');
end

if nargin>=16&&~isempty(ind)
    fprintf(fid, '\n%%Induction motor Data: %d',size(ind, 1));
    fprintf(fid, '\nInd.con = [...');
    for i = 1:size(ind, 1)
        fprintf(fid, '\n%6d %8.4g %8.4g %8.4g %6d %6d %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %2d %2d', ind(i,:));
    end
    fprintf(fid, '\n];\n');   
else
    fprintf(fid, '\nInd.con = [];\n');
end

if nargin>=17&&~isempty(syn)
    fprintf(fid, '\n%%Synchronous machine Data: %d',size(syn, 1));
    fprintf(fid, '\nSyn.con = [...');
    for i = 1:size(syn, 1)
        fprintf(fid, '\n%6d %8.4g %8.4g %8.4g %6d %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %6d %2u', syn(i,:));
    end
    fprintf(fid, '\n];\n');   
else
    fprintf(fid, '\nSyn.con = [];\n'); 
end

if nargin>=18&&~isempty(exc)
    fprintf(fid, '\n%%Exciter Data: %d',size(exc, 1));
    fprintf(fid, '\nExc.con = [...');
    for i = 1:size(exc, 1)
        fprintf(fid, '\n%6d %6d %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %2u', exc(i,:));
    end
    fprintf(fid, '\n];\n');  
else
    fprintf(fid, '\nExc.con = [];\n');   
end

if nargin>=19&&~isempty(tg)
    fprintf(fid, '\n%%Turbine governor Data: %d',size(tg, 1));
    fprintf(fid, '\nTg.con = [...');
    for i = 1:size(tg, 1)
        fprintf(fid, '\n%6d %6d %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %2u', tg(i,:));
    end
    fprintf(fid, '\n];\n'); 
else
    fprintf(fid, '\nTg.con = [];\n');    
end

fclose(fid);

end