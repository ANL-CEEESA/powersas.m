function [Base, w, FtBus, time_step, time_step_on, TimeRange1, FtTime, TotalTime, TimeRange2, TimeRange3] = set_simulation_parameters()
% This function sets the simulation parameters.

%% settings 
Base   = 100;   % baseMVA  
w      = 120*pi;
FtBus  = 7;  % grounding fault bus
time_step    = 1*10^(-6); % pre and post fault
% time_step_on = 2*10^(-8); % during fault
time_step_on = 1*10^(-6); % during fault
TimeRange1 = [0,0.1];     % pre-fault simulation time (second)
FtTime     = 5/60; % Fault lasting time (second) (5/60 is five cycles)
TotalTime  = 0.5; %1;    % Total simulation time (second)
TimeRange2 = [TimeRange1(2),TimeRange1(2)+FtTime]; % fault-on simulation period
TimeRange3 = [TimeRange2(2),TotalTime]; % post-fault simulation time
