function [Vt,dta,Id, Iq, Edp, Eqp, Edpp, Eqpp, Pm, Efd, Sys, omg, V1, efd, P1, P2, Epp] = generator_states(mpc, Sys, Base);
%GENERATOR_STATES calculates generator states for a given power system

%% generator inside states
Vtx       = mpc.bus(:,8).*cos(mpc.bus(:,9)*pi/180);
Vty       = mpc.bus(:,8).*sin(mpc.bus(:,9)*pi/180);  
Vt        = Vtx + 1i*Vty; % terminal voltage phasor value
Vt_gen    = Vt(Sys.GenIdx,:); % terminal voltage of generator bus
GenPQ     = (mpc.gen(:,2)+1i*mpc.gen(:,3))./Base;
It_gen    = conj(GenPQ./Vt_gen); % current injection

Zbranch = mpc.branch(1,3)+1i*mpc.branch(1,4);
Ibr=(Vt(1,:)-Vt(5,:))/Zbranch;

EQ        = Vt_gen + (Sys.Ra+1i*Sys.Xq).*It_gen; 
dta       = angle(EQ);   % rotor angle, angle  of Eq, not angle of E''
omg       = zeros(length(Sys.GenIdx),1);  % rotor speed

% next step, because the past calculation starts from Vt It, which are grid reference based
% thus the E' E'' are also grid state, we need to change it to d,q state
% from grid to d,q, multiply exp(-j(delta-pi/2))
Vt_gen    = Vt_gen.*exp(1i*(pi/2-dta));        % d-q based
It_gen    = It_gen.*exp(1i*(pi/2-dta)); 
Id        = real(It_gen);  
Iq        = imag(It_gen);  

% Ed''= real(E''), Eq''=imag(E'')
Ep        = Vt_gen + Sys.Ra.*It_gen -Sys.Xqp.*Iq+1i*Sys.Xdp.*Id; % d-q based
Edp       = real(Ep);  % real of E' , d-q based 
Eqp       = imag(Ep);  % imag of E' , d-q based 
Epp       = Vt_gen + Sys.Ra.*It_gen -Sys.Xqpp.*Iq+1i*Sys.Xdpp.*Id; % d-q based
Edpp      = real(Epp);  % real of E'' 
Eqpp      = imag(Epp);  % imag of E'' 

Pe        = real(GenPQ) + (abs(It_gen).^2).*Sys.Ra; % electrical power
Pm        = Pe + Sys.D.*omg./Sys.w0;      % mechanical power    

% exciter field voltage from 6th-order equation, derivative=0
Efd  = (Sys.Xd-Sys.Xdpp)./(Sys.Xdp-Sys.Xdpp).*Eqp - (Sys.Xd-Sys.Xdp)./(Sys.Xdp-Sys.Xdpp).*Eqpp;
efd  = Efd.*Sys.Rfd./Sys.Lad;  % transfer to efd
% set Vref for exciter equation
% from transfer equation, steady state Vref-Et-V1=0
Sys.Vref     = abs(Epp)+ abs(Efd)./Sys.Ek;   % abs(Vt_gen) changed to abs(Epp)
V1           = abs(Efd)./Sys.Ek;  % V1 is a mid-term in exciter
% from transfer equation,steady state Pref=P1*Gr; P1=P2; Pm=P2;
Sys.Pref  = Pm.*Sys.Gr;
P1        = Pm;
P2        = Pm;