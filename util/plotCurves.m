function plotCurves(figId,t,stateCurve,SysDataBase,variable,subIdx,lw,dot,fontSize,width,height,sizeRatio)
% Plot the system state trajectories
%
% FUNCTION plotCurves
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
%   figId - ID of figure
%   t - sequence of time
%	stateCurve - sequence of system states
%	SysDataBase - system data structure
%	variable - String denoting the variable to plot
%	subIdx - Pick part of variables to plot
%	lw - line width
%	dot - line type
%	fontSize - Font size of labels
%	width - Width (pixels) of figure
%	height - Height (pixels) of figure
%	sizeRatio - If specified, the width and height are determined by a ratio of the screen size
%	
[nState,idxs]...
    =getIndexDyn(SysDataBase);
if nargin<6
    subIdx=[];
end
if nargin<7||isempty(lw)
    lw=1;
end
if nargin<8||isempty(dot)
    ln='-';
else
    if dot
        ln='-o';
    else
        ln='-';
    end
end
if nargin<9||isempty(fontSize)
    fontSize=12;
end
if nargin<12||isempty(sizeRatio)
    sizeRatio=0.7;
end
screensize = get( 0, 'Screensize' );
if nargin<10||isempty(width)
    width=floor(screensize(3)*sizeRatio);
end
if nargin<11||isempty(height)
    height=floor(screensize(4)*sizeRatio);
end
mksz=lw*3;
xdata=t;
ydata=[];
ylbl='';
if strcmpi(variable,'list') || strcmpi(variable,'help')
    disp('------ Variable List ------')
    disp('v:Bus voltage magnitude(pu)')
    disp('a:Bus voltage angle(rad)')
    disp('delta:Gen rotor angle(rad)')
    disp('omega:Gen rotor speed(pu)')
    disp('s:Motor slip')
    disp('eq1:Gen Eq1(pu)')
    disp('eq2:Gen Eq2(pu)')
    disp('ed1:Gen Ed1(pu)')
    disp('ed2:Gen Ed2(pu)')
    disp('psid:Gen Psid(pu)')
    disp('psiq:Gen Psiq(pu)')
    disp('pg:Gen Pg(pu)')
    disp('ef:Gen Ef(pu)')
    disp('vavrm:AVR Variable m(pu)')
    disp('vavrr:AVR Variable r(pu)')
    disp('vavrf:AVR Variable f(pu)')
    disp('vref:AVR Reference volt(pu)')
    disp('tgovg:TG Variable g(pu)')
    disp('tgovm:TG Variable m(pu)')
    disp('tmech:TG Tmech0(pu)')
    disp('f:bus frequency deviation(pu)')
    disp('dpg:AGC power')
    disp('qplt')
    disp('vg')
    disp('--------------------------')
elseif strcmpi(variable,'v')
    if isempty(subIdx);idx=idxs.vIdx;else idx=idxs.vIdx(subIdx);end
    ydata=abs(stateCurve(idx,:));ylbl='Bus voltage magnitude(pu)';
elseif strcmpi(variable,'a')
    if isempty(subIdx);idx=idxs.vIdx;else idx=idxs.vIdx(subIdx);end
    ydata=angle(stateCurve(idx,:));ylbl='Bus voltage angle(rad)';    
elseif strcmpi(variable,'delta')
    if isempty(subIdx);idx=idxs.deltaIdx;else idx=idxs.deltaIdx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='Gen rotor angle(rad)';    
elseif strcmpi(variable,'omega')
    if isempty(subIdx);idx=idxs.omegaIdx;else idx=idxs.omegaIdx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='Gen rotor speed(pu)'; 
elseif strcmpi(variable,'s')
    if isempty(subIdx);idx=idxs.sIdx;else idx=idxs.sIdx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='Motor slip';     
elseif strcmpi(variable,'eq1')
    if isempty(subIdx);idx=idxs.eq1Idx;else idx=idxs.eq1Idx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='Gen Eq1(pu)';     
elseif strcmpi(variable,'eq2')
    if isempty(subIdx);idx=idxs.eq2Idx;else idx=idxs.eq2Idx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='Gen Eq2(pu)';  
elseif strcmpi(variable,'ed1')
    if isempty(subIdx);idx=idxs.ed1Idx;else idx=idxs.ed1Idx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='Gen Ed1(pu)';  
elseif strcmpi(variable,'ed2')
    if isempty(subIdx);idx=idxs.ed2Idx;else idx=idxs.ed2Idx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='Gen Ed2(pu)';     
elseif strcmpi(variable,'psid')
    if isempty(subIdx);idx=idxs.psidIdx;else idx=idxs.psidIdx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='Gen Psid(pu)'; 
elseif strcmpi(variable,'psiq')
    if isempty(subIdx);idx=idxs.psiqIdx;else idx=idxs.psiqIdx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='Gen Psiq(pu)'; 
elseif strcmpi(variable,'pg')
    if isempty(subIdx);idx=idxs.pgIdx;else idx=idxs.pgIdx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='Gen Pg(pu)'; 
elseif strcmpi(variable,'ef')
    if isempty(subIdx);idx=idxs.efIdx;else idx=idxs.efIdx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='Gen Ef(pu)';    
elseif strcmpi(variable,'vavrm')
    if isempty(subIdx);idx=idxs.vavrmIdx;else idx=idxs.vavrmIdx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='AVR Variable m(pu)';    
elseif strcmpi(variable,'vavrr')
    if isempty(subIdx);idx=idxs.vavrrIdx;else idx=idxs.vavrrIdx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='AVR Variable r(pu)';     
elseif strcmpi(variable,'vavrf')
    if isempty(subIdx);idx=idxs.vavrfIdx;else idx=idxs.vavrfIdx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='AVR Variable f(pu)';     
elseif strcmpi(variable,'vref')
    if isempty(subIdx);idx=idxs.vavrrefIdx;else idx=idxs.vavrrefIdx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='AVR Reference volt(pu)';   
elseif strcmpi(variable,'tgovg')
    if isempty(subIdx);idx=idxs.tgovgIdx;else idx=idxs.tgovgIdx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='TG Variable g(pu)'; 
elseif strcmpi(variable,'tgovm')
    if isempty(subIdx);idx=idxs.tgovmIdx;else idx=idxs.tgovmIdx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='TG Variable m(pu)'; 
elseif strcmpi(variable,'tmech')
    if isempty(subIdx);idx=idxs.tmechIdx;else idx=idxs.tmechIdx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='TG Tmech0(pu)'; 
elseif strcmpi(variable,'f')
    if isempty(subIdx);idx=idxs.fIdx;else idx=idxs.fIdx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='delta f(pu)'; 
elseif strcmpi(variable,'dpg')
    if isempty(subIdx);idx=idxs.dpgIdx;else idx=idxs.dpgIdx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='AGC delta pg(pu)';
elseif strcmpi(variable,'qplt')
    if isempty(subIdx);idx=idxs.qpltIdx;else idx=idxs.qpltIdx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='qplt(pu)'; 
elseif strcmpi(variable,'vg')
    if isempty(subIdx);idx=idxs.vgIdx;else idx=idxs.vgIdx(subIdx);end
    ydata=real(stateCurve(idx,:));ylbl='vg(pu)';  
end
if ~isempty(ydata)
    if ~isempty(figId)
        p=figure(figId);
    else
        p=figure();
    end
    movegui(p,'center');
    pos=get(p,'OuterPosition');
    centerx=pos(1)+pos(3)/2;
    centery=pos(2)+pos(4)/2;
    ymin=min(min(ydata));
    ymax=max(max(ydata));
    yRange=ymax-ymin;
    if yRange<1e-3
        yRange=1e-3;        
    end
    xmin=min(xdata);
    xmax=max(xdata);
    xDispRange=xmax-xmin;
    if xDispRange<1e-6
        xDispRange=1e-6; 
    else
        xDispRange=0;
    end
    if length(xdata)==1
        ln='-o';
    end
    plot(xdata,ydata,ln,'LineWidth',lw,'MarkerSize',mksz);
    xlabel('Time (s)','FontSize', fontSize),...
        ylabel(ylbl,'FontSize', fontSize),...
        axis([xmin-0.05*xDispRange,xmax+0.05*xDispRange,ymin-0.05*yRange,ymax+0.05*yRange]);
    set(p,'units','pixels','OuterPosition',[floor(centerx-width/2),floor(centery-height/2),width,height]);
end

end