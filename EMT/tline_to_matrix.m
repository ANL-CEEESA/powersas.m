function [R_line, L_line, C_line] = tline_to_matrix(R_line_positive, wL_line_positive, B_line_positive, LC_fre)
% Unit is per unit per km, pu/km
% Assume all the lines have the same per unit parameters
   
R_line_zero      = R_line_positive*6;  % from PSCAD line parameters, we can generally get this relationship  
wL_line_zero     = wL_line_positive*3;
L_line_positive  = wL_line_positive/(LC_fre); 
L_line_zero      = wL_line_zero/(LC_fre); 
B_line_zero      = B_line_positive/1.5;
C_line_positive  = 0.5*B_line_positive/(LC_fre); % multiply 0.5 because it is divided in to two parts
C_line_zero      = 0.5*B_line_zero/(LC_fre); 

% three phase self and mutual parameters
% Z_zero=Zs+2Zm,  Z_positive=Zs-Zm
% so we get:  Z_self=1/3*(Z_zero+2*Z_positive), Z_mutual=1/3*(Z_zero-Z_positive)
Rs_line=1/3*(R_line_zero+2*R_line_positive);  Rm_line=1/3*(R_line_zero-R_line_positive); 
Ls_line=1/3*(L_line_zero+2*L_line_positive);  Lm_line=1/3*(L_line_zero-L_line_positive);          
Cs_line=1/3*(C_line_zero+2*C_line_positive);  Cm_line=1/3*(C_line_zero-C_line_positive); 

% R L C three phase matrix  
R_line=[Rs_line, Rm_line, Rm_line; Rm_line, Rs_line, Rm_line; Rm_line, Rm_line, Rs_line];
L_line=[Ls_line, Lm_line, Lm_line; Lm_line, Ls_line, Lm_line; Lm_line, Lm_line, Ls_line];
C_line=[Cs_line, Cm_line, Cm_line; Cm_line, Cs_line, Cm_line; Cm_line, Cm_line, Cs_line];

end
