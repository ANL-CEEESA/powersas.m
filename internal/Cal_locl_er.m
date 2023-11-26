function [LTE,LTE2,LTE3] = Cal_locl_er(x_ori,k,h)
% Calculation of local error
%
% FUNCTION Cal_locl_er
%
% Author: Kaiyang Huang <khuang@anl.gov>
%
% Copyright (C) 2023, UChicago Argonne, LLC. All rights reserved.
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

%   x_ori - Initial system state
%   xNew  - Updated system state (algebraic variables to be solved)
%   SysDataNew - (optional) updated system data
%	
% OUTPUT
%	Vx, Wx, Qx, fx - Last HE coefficients of algebraic variables
%	xx - Solved updated system state
%   finalAlpha - The ending length of this segment of simulation
%   alphas - Record of alphas
%   diffRec - A list of errors
%

    x       = x_ori(:,k-3:k);
    Xh      = x.*h.^(k-4:k-1);
    x_new   = polyvalVec(x_ori,h);
    eta     = 0.001;
    infcof  = max(abs(Xh(:,end))./(abs(x_new)+eta));
    trc     = infcof^(1/k);
    LTE     = infcof/(1-trc);

    LTE3    = max(abs(Xh(:,end-1)./(x_new+eta)));

    % relative error used to change order
    e1      = abs(Xh(:,end-1)./Xh(:,end));
    e2      = abs(Xh(:,end-2)./Xh(:,end));
    e3      = abs(Xh(:,end-3)./Xh(:,end-1));
    % delete nan induced by zero/zero
    f1      = find(Xh(:,end)==0);
    f2      = find(Xh(:,end-1)==0);
    f2      = f2(:);
    e1(f1)  = 0;
    e2(f1)  = 0;
    e3(f2)  = 0;
    LTE2    = infcof/min([norm(e1,Inf), sqrt(norm(e2,Inf)), sqrt(norm(e3,Inf))]);

end

