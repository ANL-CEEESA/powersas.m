function [R_line_positive, wL_line_positive, B_line_positive, int_t, LC_fre, t_base] = setDefaultParameters()

R_line_positive = 0.0001; 
wL_line_positive = 0.001; 
B_line_positive = 0.00175;
int_t = 0;   % could change to any number
LC_fre = 120*pi; % divide all L by w, consider time base is 1/120/pi
t_base =1/(120*pi);
end
