function Sys = update_system_parameters(Sys, Laa0_pp, Lab0_pp)
% This function updates the system parameters including Emax, Emin, Gen_R and Lline_3Phase.

Sys.Emax=Sys.Emax.*Sys.Rfd./Sys.Lad;
Sys.Emin=Sys.Emin.*Sys.Rfd./Sys.Lad;

for k=1:1:length(Sys.GenIdx)
    Sys.Gen_R(:,:,k)= [Sys.Ra(k), 0, 0; 0, Sys.Ra(k), 0; 0, 0, Sys.Ra(k)];
    L_stator(:,:,k)= [Laa0_pp(k) Lab0_pp(k) Lab0_pp(k); Lab0_pp(k) Laa0_pp(k) Lab0_pp(k); Lab0_pp(k) Lab0_pp(k) Laa0_pp(k)]/120/pi;
    % because Sys.Laa0_pp is already updated
    Sys.Lline_3Phase(:,:,k)  = Sys.Lline_3Phase(:,:,k)+ L_stator(:,:,k);
end

