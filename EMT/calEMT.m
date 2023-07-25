
function [res] = calEMT(testSystem)
%%
[orderOfDT, timeStepOfDT] = setSASParameters();

[Base, w, FtBus, time_step, time_step_on, TimeRange1, FtTime, TotalTime, TimeRange2, TimeRange3] = set_simulation_parameters();

%%
[R_line_positive, wL_line_positive, B_line_positive, int_t, LC_fre, t_base] = setDefaultParameters();

[mpc, SUCESS]     = runpf(testSystem);

Sys = myFunction(mpc);

[Sys,Laa0_pp,Lab0_pp,Laa2_pp] = transformParameters(Sys,t_base);

[Vt,dta,Id, Iq, Edp, Eqp, Edpp, Eqpp, Pm, Efd, Sys, omg, V1, efd, P1, P2, Epp] = generator_states(mpc, Sys, Base);

[Lambda_F, Lambda_D, Lambda_Q1, Lambda_Q2] = calculateLambda(Sys, Id, Iq, efd);

[R_line, L_line, C_line] = tline_to_matrix(R_line_positive, wL_line_positive, B_line_positive,LC_fre);

[L_tf_positive L_tf LoadPQ Zload R_load L_load Sys LoadShunt YShunt C_Shunt CShunt] = fun1(Sys,LC_fre,mpc,Vt,Base);

[Vbus_3Phase branch_number c Sys Iline_3Phase ind_line1 ind_line2 Zbranch Ibranch ...
    Amp Ang Iload_3Phase ind_load theta] = fun2(Sys, mpc, dta, int_t, w, Zload);

L_temp = calculate_L_temp(L_tf, L_line, Laa0_pp, Lab0_pp);

[Sys] = calculateSystemMatrices(mpc, Sys, CShunt, R_line, wL_line_positive, branch_number, L_temp, C_line);

[Vbus, Iline, Iload] = convert_to_single_vector(Vbus_3Phase, Iline_3Phase, Iload_3Phase, Sys, branch_number);

Sys = update_system_parameters(Sys, Laa0_pp, Lab0_pp);

x_ini=[dta; theta; omg; Lambda_F; Lambda_D; Lambda_Q1; Lambda_Q2; V1; efd; P1; P2;...
       Vbus(3*length(Sys.GenIdx)+1:end); Iline; Iload; abs(Epp)]; 
 

%%
options = odeset('RelTol',1e-4,'AbsTol',1e-4); % limit the error
[t_pre,x_pre]= ode23t(@funEMT_pre,TimeRange1(1):time_step:TimeRange1(2),x_ini,options,Sys);  

x_ini_on = x_pre(end,:)'; 
[t_on,x_on] = ode23t(@funEMT_on,TimeRange2(1):time_step_on:TimeRange2(2),x_ini_on,options,Sys);  

x_ini_post  = x_on(end,:)';
[t_dt,x_dt]=DTsolver(TimeRange3,timeStepOfDT,x_ini_post,orderOfDT,Sys); %x_ini_post

t_post = t_dt;
x_post = x_dt';
res.t_final   = [t_pre;  t_on(2:end);    t_post(2:end)];
res.sol_final = [x_pre;  x_on(2:end,:);  x_post(2:end,:)];   
return