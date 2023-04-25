function [L_tf_positive L_tf LoadPQ Zload R_load L_load Sys LoadShunt YShunt C_Shunt CShunt] = fun1(Sys,LC_fre,mpc,Vt,Base)
%% transformer parameters, per unit
% Assume all the transformers have the same impedance
L_tf_positive=mpc.branch(Sys.GenIdx(1),4)/(LC_fre); %positive sequence leakage reactance
L_tf=[L_tf_positive  0  0;  0  L_tf_positive  0;  0  0  L_tf_positive];

%% Constant PQ load parameter, per unit
% PQ load
LoadPQ    = (mpc.bus(Sys.LoadIdx,3)+1i*mpc.bus(Sys.LoadIdx,4))./Base;
Zload     = conj((abs(Vt(Sys.LoadIdx,:)).^2)./LoadPQ); % constant Z load
R_load    = real(Zload);
L_load    = imag(Zload)/(LC_fre);
Sys.Rload_3Phase = zeros(3,3,length(R_load));
Sys.Lload_3Phase = Sys.Rload_3Phase;
for k=1:1:size(R_load)
  Sys.Rload_3Phase(:,:,k)=[R_load(k)  0  0; 0  R_load(k)  0;  0  0  R_load(k)];
  Sys.Lload_3Phase(:,:,k)=inv([L_load(k)  0  0; 0  L_load(k)  0;  0  0  L_load(k)]);
end

% fixed shunt
LoadShunt = (mpc.bus(Sys.ShuntIdx,5)-1i*mpc.bus(Sys.ShuntIdx,6))./Base; % fixed shunt load, use - because it generate reactive power
YShunt    = conj(LoadShunt./1.^2 );   % constant Z load, do not use abs(Vt(LoadIdx,:)).^2, because voltage should be 1
C_Shunt   = imag(YShunt)/(LC_fre); 
CShunt    = zeros(3,3,length(C_Shunt));
for k=1:1:length(C_Shunt)
  CShunt(:,:,k)=[C_Shunt(k)  0  0; 0  C_Shunt(k)  0;  0  0  C_Shunt(k)];
end