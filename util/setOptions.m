function options=setOptions(varargin)
% General function for setting options in key-value pairs
%
% FUNCTION setOptions
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

options=[];
if nargin<1
    options=[];
else
    idx=1;
    if isstruct(varargin{1})
        options=varargin{1};
        idx=idx+1;
    elseif ~ischar(varargin{1})
        idx=idx+1;
    end
    while idx<=nargin
        opt=varargin{idx};
        idx=idx+1;
        if idx<=nargin
            val=varargin{idx};
            idx=idx+1;
            if ischar(opt)
                try
                    options.(opt)=val;
                catch ME
                    disp(['setOptions() Argument ',num2str(idx-2),' & ',num2str(idx-1),': ',ME.message]);
                end
            else
                disp(['setOptions() Argument ',num2str(idx-2),': option name should be a string. Check input arguments.']);
            end
        else
            disp(['setOptions() Argument ',num2str(idx-1),': option and value should match in pairs. Check input arguments.']);
        end
    end
end
end