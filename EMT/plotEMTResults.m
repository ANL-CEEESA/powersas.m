close all

% t_final   = [t_pre;  t_on(2:end);    t_post(2:end)];
% sol_final = [x_pre;  x_on(2:end,:);  x_post(2:end,:)];
% % t_final   = [t_pre;  t_on(2:end)];
% % sol_final = [x_pre;  x_on(2:end,:)];
% % t_final   = [t_pre];
% % sol_final = [x_pre];

%% relative rotor angles
figure(1)
% EMT result
plot(t_final, (sol_final(:,1)-sol_final(:,1) )*180/pi,'r--','LineWidth',3) 
hold on
plot(t_final, (sol_final(:,2)-sol_final(:,1) )*180/pi,'g--','LineWidth',3) 
hold on
plot(t_final, (sol_final(:,3)-sol_final(:,1) )*180/pi,'b--','LineWidth',3)
hold on
plot(t_final, (sol_final(:,4)-sol_final(:,1) )*180/pi,'b--','LineWidth',3)
hold on
% ylim([-60,10])
xlabel('Time (s)')
ylabel('Relative rotor angles (degree)')
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 12);
title('Relative rotor angles')
legend('EMT G1','EMT G2','EMT G3','EMT G4')

%% frequency deviation
figure(2)
% EMT result
plot(t_final, sol_final(:,9),'r--','LineWidth',3) 
hold on
plot(t_final, sol_final(:,10),'g--','LineWidth',3) 
hold on
plot(t_final, sol_final(:,11),'b--','LineWidth',3) 
hold on
plot(t_final, sol_final(:,12),'m--','LineWidth',3) 
hold on
% ylim([-60,10])
% xlim([0,0.25])
xlabel('Time (s)')
ylabel('Frequency deviation (p.u.)')
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 12);
title('Frequency deviation')
legend('EMT G1','EMT G2','EMT G3','EMT G4')

%% Exciter field voltage
figure(3)
% EMT result
plot(t_final, sol_final(:,33),'r--','LineWidth',3) 
hold on
plot(t_final, sol_final(:,34),'g--','LineWidth',3) 
hold on
plot(t_final, sol_final(:,35),'b--','LineWidth',3) 
hold on
plot(t_final, sol_final(:,36),'m--','LineWidth',3) 
xlabel('Time (s)')
xlim([0,0.4])
ylabel('Exciter field voltage (p.u.)')
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 12);
title('Exciter field voltage')
legend('EMT G1','EMT G2','EMT G3','EMT G4')

%% P1
figure(4)
% EMT result
plot(t_final, sol_final(:,37),'r--','LineWidth',3) 
hold on
plot(t_final, sol_final(:,38),'g--','LineWidth',3) 
hold on
plot(t_final, sol_final(:,39),'b--','LineWidth',3) 
hold on
plot(t_final, sol_final(:,40),'m--','LineWidth',3) 

xlabel('Time (s)')
ylabel('P1(p.u.)')
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 12);
title('P1')
legend('EMT G1','EMT G2','EMT G3','EMT G4')

%% Voltage at bus 7
figure(5)
% EMT result
plot(t_final, sol_final(:,51),'r','LineWidth',1) 
hold on
plot(t_final, sol_final(:,52),'g','LineWidth',1) 
hold on
plot(t_final, sol_final(:,53),'b','LineWidth',1) 

xlabel('Time (s)')
ylabel('Voltage at bus 7(p.u.)')
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 12);
title('Voltage at bus 7')
legend('EMT phase A','EMT phase B','EMT phase C')

%% Voltage at bus 8
figure(6)
% EMT result
plot(t_final, sol_final(:,54),'r','LineWidth',1) 
hold on
plot(t_final, sol_final(:,55),'g','LineWidth',1) 
hold on
plot(t_final, sol_final(:,56),'b','LineWidth',1) 

xlabel('Time (s)')
ylabel('Voltage at bus 8(p.u.)')
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 12);
title('Voltage at bus 8')
legend('EMT phase A','EMT phase B','EMT phase C')

%% Voltage at bus 10
figure(7)
% EMT result
plot(t_final, sol_final(:,60),'r','LineWidth',1)
hold on
plot(t_final, sol_final(:,61),'g','LineWidth',1)
hold on
plot(t_final, sol_final(:,62),'b','LineWidth',1)

xlabel('Time (s)')
ylabel('Voltage at bus 10(p.u.)')
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 12);
title('Voltage at bus 10')
legend('EMT phase A','EMT phase B','EMT phase C')

%% Voltage at bus 6
figure(8) 
% EMT result
plot(t_final, sol_final(:,48),'r','LineWidth',1) 
hold on
plot(t_final, sol_final(:,49),'g','LineWidth',1) 
hold on
plot(t_final, sol_final(:,50),'b','LineWidth',1) 

xlabel('Time (s)')
ylabel('Voltage at bus 6(p.u.)')
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 12);
title('Voltage at bus 6')
legend('EMT phase A','EMT phase B','EMT phase C')

%% current on branch 5-6
figure(9) 
% EMT result
plot(t_final, sol_final(:,78),'r','LineWidth',1) 
hold on
plot(t_final, sol_final(:,79),'g','LineWidth',1) 
hold on
plot(t_final, sol_final(:,80),'b','LineWidth',1) 

xlabel('Time (s)')
ylabel('Current on branch 5-6(p.u.)')
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 12);
title('Current on branch 5-6')
legend('EMT phase A','EMT phase B','EMT phase C')

%% current on branch 6-7
figure(10) 
% EMT result
plot(t_final, sol_final(:,81),'r','LineWidth',1) 
hold on
plot(t_final, sol_final(:,82),'g','LineWidth',1) 
hold on
plot(t_final, sol_final(:,83),'b','LineWidth',1) 

xlabel('Time (s)')
ylabel('Current on branch 6-7(p.u.)')
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 12);
title('Current on branch 6-7')
legend('EMT phase A','EMT phase B','EMT phase C')

%% current on branch 7-8
figure(11) 
% EMT result
plot(t_final, sol_final(:,84),'r','LineWidth',1) 
hold on
plot(t_final, sol_final(:,85),'g','LineWidth',1) 
hold on
plot(t_final, sol_final(:,86),'b','LineWidth',1) 

xlabel('Time (s)')
ylabel('Current on branch 7-8(p.u.)')
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 12);
title('Current on branch 7-8')
legend('EMT phase A','EMT phase B','EMT phase C')

%% current on branch 9-10
figure(12) 
% EMT result
plot(t_final, sol_final(:,96),'r','LineWidth',1) 
hold on
plot(t_final, sol_final(:,97),'g','LineWidth',1) 
hold on
plot(t_final, sol_final(:,98),'b','LineWidth',1) 

xlabel('Time (s)')
ylabel('Current on branch 9-10(p.u.)')
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 12);
title('Current on branch 9-10')
legend('EMT phase A','EMT phase B','EMT phase C')

%% Voltage at bus 6
figure(13) 
% EMT result
plot(t_final, sol_final(:,108),'r','LineWidth',1) 
hold on
plot(t_final, sol_final(:,109),'g','LineWidth',1) 
hold on
plot(t_final, sol_final(:,110),'b','LineWidth',1) 
hold on
plot(t_final, sol_final(:,111),'b','LineWidth',1) 

xlabel('Time (s)')
ylabel('Voltage at bus 6(p.u.)')
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 12);
title('Voltage at bus 6')
legend('Gen1','Gen2','Gen3','Gen4')
