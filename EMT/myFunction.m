function Sys = myFunction(mpc)

%% Generator parameters
Sys.H    = [6.5 6.5 6.175 6.175]'*9;  % base of 900 MW to base 100 MW
Sys.Xl   = [0.2 0.2 0.2 0.2]'/9;
Sys.Xd   = [1.8 1.8 1.8 1.8]'/9;
Sys.Xq   = [1.7 1.7 1.7 1.7]'/9;
Sys.Xdp  = [0.3 0.3 0.3 0.3]'/9;
Sys.Xqp  = [0.55 0.55 0.55 0.55]'/9;
Sys.Xdpp = [0.25 0.25 0.25 0.25]'/9;
Sys.Xqpp = [0.25 0.25 0.25 0.25]'/9;
Sys.Ra   = [0.0025 0.0025 0.0025 0.0025]'/9;
Sys.D    = [45 45 45 45]'; % 0 on Kunder book
Sys.Td0p = [8 8 8 8]';
Sys.Tq0p = [0.4 0.4 0.4 0.4]';
Sys.Td0pp= [0.03 0.03 0.03 0.03]';
Sys.Tq0pp= [0.05 0.05 0.05 0.05]';

% add by myself
Sys.L0  =  Sys.Xl*2; % L0 in VBR synchronous generator model
Sys.T_vt=  0.05; % time constant for terminal voltage measurement

%% transfer generator parameters EMT model parameters
% Paremeter_transform  % provide parameters for VBR generator model

%% contollers
% exciter
Sys.Eta=1;   Sys.Etb=10;  Sys.Ek=100;  Sys.Ete=0.1;  Sys.Emax=3;  Sys.Emin=0;
% TGov
Sys.Gr= 0.33;   Sys.Gt1=2;  Sys.Gt2=3;  Sys.Gt3=15;  Sys.Gdt=0.4;  Sys.GVmax=10;  Sys.GVmin=0;  

%% Index
Sys.w0           = 1.0;   % nominal frequency, p.u. useful in swing equation
Sys.GenIdx       = [1 2 3 4];
Sys.NonGenIdx    = [5 6 7 8 9 10 11];
Sys.LoadIdx      = [7 9];
Sys.ShuntIdx     = [7 9];
Sys.bus_number   = mpc.bus(end,1); % total number of buses

