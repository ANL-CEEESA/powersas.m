function [Vbus_3Phase branch_number c Sys Iline_3Phase ind_line1 ind_line2 Zbranch Ibranch ...
    Amp Ang Iload_3Phase ind_load theta] = fun2(Sys, mpc, dta, int_t, w, Zload)
%% current and voltage, three phase
% all the bus voltage, three phase, assume t=ini_t, and Va=sin(wt+angle)
Vbus_3Phase=zeros(3,1,Sys.bus_number);
for k=1:1:Sys.bus_number
Vbus_3Phase(:,:,k)  = [mpc.bus(k,8).*sin(w*int_t+ mpc.bus(k,9)*pi/180);...
                       mpc.bus(k,8).*sin(w*int_t+(mpc.bus(k,9)-120)*pi/180);...
                       mpc.bus(k,8).*sin(w*int_t+(mpc.bus(k,9)+120)*pi/180)];
end
           
% all the branch circuit current, three phase, assume t=0, and Ia=sin(wt+angle)
[branch_number,c] = size(mpc.branch);
Sys.branch_number=branch_number;
Iline_3Phase=zeros(3,1,branch_number);
for k=1:1:branch_number % consider all branch and load
    ind_line1=mpc.branch(k,1);
    ind_line2=mpc.branch(k,2);
    Zbranch = mpc.branch(k,3)+1i*mpc.branch(k,4);
    Ibranch= ( mpc.bus(ind_line1,8)*exp(1i*mpc.bus(ind_line1,9)*pi/180)- ...
               mpc.bus(ind_line2,8)*exp(1i*mpc.bus(ind_line2,9)*pi/180) )/Zbranch;
    Amp=abs(Ibranch);
    Ang=angle(Ibranch);
    Iline_3Phase(:,:,k)=[Amp.*sin(w*int_t+Ang);...
                         Amp.*sin(w*int_t+Ang-120*pi/180);...
                         Amp.*sin(w*int_t+Ang+120*pi/180)];
end

% all load current
Iload_3Phase=zeros(3,1,length(Sys.LoadIdx));
for k=1:1:length(Sys.LoadIdx) % consider all branch and load
    ind_load=Sys.LoadIdx(k);  
    Iload= mpc.bus(ind_load,8)*exp(1i*mpc.bus(ind_load,9)*pi/180)/Zload(k);
    Amp=abs(Iload);
    Ang=angle(Iload);
    Iload_3Phase(:,:,k)=[Amp.*sin(w*int_t+Ang);...
                         Amp.*sin(w*int_t+Ang-120*pi/180);...
                         Amp.*sin(w*int_t+Ang+120*pi/180)];
end

%% Initialization of Theta angle
theta=dta-pi+120*pi*int_t; % I found that if int_t=0, theta=dta-pi