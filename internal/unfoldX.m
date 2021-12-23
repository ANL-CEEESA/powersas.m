function [V,Q,s,d,w,eq1,eq2,ed1,ed2,psid,psiq,Pm,Ef,Vavrm,Vavrr,Vavrf,Vavrref,tgovg,tgovm,tgovmech,f,dpg,qplt,vg]=unfoldX(x,SysData)
% function for unfolding system states
%
% FUNCTION unfoldX
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

[~,idxs]...
    =getIndexDyn(SysData);

V=x(idxs.vIdx,:);
Q=x(idxs.qIdx,:);
s=x(idxs.sIdx,:);
d=x(idxs.deltaIdx,:);
w=x(idxs.omegaIdx,:);
eq1=x(idxs.eq1Idx,:);
eq2=x(idxs.eq2Idx,:);
ed1=x(idxs.ed1Idx,:);
ed2=x(idxs.ed2Idx,:);
psid=x(idxs.psidIdx,:);
psiq=x(idxs.psiqIdx,:);
Pm=x(idxs.pgIdx,:);
Ef=x(idxs.efIdx,:);
Vavrm=x(idxs.vavrmIdx,:);
Vavrr=x(idxs.vavrrIdx,:);
Vavrf=x(idxs.vavrfIdx,:);
Vavrref=x(idxs.vavrrefIdx,:);
tgovg=x(idxs.tgovgIdx,:);
tgovm=x(idxs.tgovmIdx,:);
tgovmech=x(idxs.tmechIdx,:);

f=x(idxs.fIdx,:);
dpg=x(idxs.dpgIdx,:);
qplt=x(idxs.qpltIdx,:);
vg=x(idxs.vgIdx,:);

end