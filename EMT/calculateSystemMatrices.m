function [Sys] = calculateSystemMatrices(mpc, Sys, CShunt, R_line, wL_line_positive, branch_number, L_temp, C_line);

% Calculate System Matrices
% Inputs:
% - mpc: power system data
% - Sys: structure containing system matrices
% - CShunt: shunt capacitance matrix
% - L_tf: transformer leakage inductance matrix
% - R_line: transmission line resistance
% - L_line: transmission line inductance
% - wL_line_positive: positive sequence angular frequency
% - Laa0_pp: stator inductance matrix element Laa0, prime-prime frame
% - Lab0_pp: stator inductance matrix element Lab0, prime-prime frame
% Output:
% - Sys: structure containing calculated system matrices

L_tf = L_temp.L_tf;
L_line = L_temp.L_line;
Laa0_pp = L_temp.Laa0_pp;
Lab0_pp = L_temp.Lab0_pp;


Sys.Cnode_3Phase=zeros(3,3,Sys.bus_number);
% T line
for k=1:1:branch_number
    ind_line1=mpc.branch(k,1);
    ind_line2=mpc.branch(k,2);
    % if it is transmission line branch    
    if ismember(ind_line1,Sys.NonGenIdx)
    term=mpc.branch(k,4)/wL_line_positive*C_line;
    Sys.Cnode_3Phase(:,:,ind_line1)  = Sys.Cnode_3Phase(:,:,ind_line1) + term;
    Sys.Cnode_3Phase(:,:,ind_line2)  = Sys.Cnode_3Phase(:,:,ind_line2) + term;
    end   
end
% Fix shunt
for k=1:1:length(Sys.ShuntIdx)
    Sys.Cnode_3Phase(:,:,Sys.ShuntIdx(k))= Sys.Cnode_3Phase(:,:,Sys.ShuntIdx(k))+ CShunt(:,:,k);   
end
% transformed into inverse matrices
for k=1:1:Sys.bus_number
    Sys.Cnode_3Phase(:,:,k)= inv(Sys.Cnode_3Phase(:,:,k));   
end

% all the Branch R L, consider transformer, T line, and load
Sys.Rline_3Phase=zeros(3,3,branch_number);
Sys.Lline_3Phase=zeros(3,3,branch_number);
L_stator=zeros(3,3,length(Sys.GenIdx));
% T line
for k=1:1:branch_number
    ind_line1=mpc.branch(k,1);
    ind_line2=mpc.branch(k,2);
    if ismember(ind_line1,Sys.GenIdx)
    L_stator(:,:,k)= [Laa0_pp(k) Lab0_pp(k) Lab0_pp(k); Lab0_pp(k) Laa0_pp(k) Lab0_pp(k); Lab0_pp(k) Lab0_pp(k) Laa0_pp(k)]/120/pi;
    Sys.Rline_3Phase(:,:,k)  = zeros(3,3);
    Sys.Lline_3Phase(:,:,k)  = inv(L_tf+L_stator(:,:,k));
    % if it is transmission line branch 
    elseif ismember(ind_line1,Sys.NonGenIdx) 
    Sys.Rline_3Phase(:,:,k)  = mpc.branch(k,4)/wL_line_positive*R_line;
    Sys.Lline_3Phase(:,:,k)  = inv(mpc.branch(k,4)/wL_line_positive*L_line);
    end
end

for k=1:1:length(Sys.GenIdx)
    Sys.Gen_R(:,:,k)= [Sys.Ra(k), 0, 0; 0, Sys.Ra(k), 0; 0, 0, Sys.Ra(k)];
    L_stator(:,:,k)= [Laa0_pp(k) Lab0_pp(k) Lab0_pp(k); Lab0_pp(k) Laa0_pp(k) Lab0_pp(k); Lab0_pp(k) Lab0_pp(k) Laa0_pp(k)]/120/pi;
    % because Sys.Laa0_pp is already updated
    Sys.Lline_3Phase(:,:,k)  = Sys.Lline_3Phase(:,:,k)+ L_stator(:,:,k);
end