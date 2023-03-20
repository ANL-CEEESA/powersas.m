casename = 'ieee14';

% 
[mpc, warnings] = psse2mpc(horzcat(casename,'.raw'),horzcat(casename,'.m'));
check = matpower2psat(horzcat(casename,'.m'), pwd);

return