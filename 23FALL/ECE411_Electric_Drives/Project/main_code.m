%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ECE 411 Project Hyeongmeen Baik
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

addpath('fig')
J = 0.2
p = 8
V_R = 72
I_R = 85
N_R = 3700
N_nl = 4259
wm_R = 3700/60*2*pi;
wm_nl = N_nl/60*2*pi;
T_R = 14;

k_T = T_R/I_R %%%%% DC Machine Constant
% k_E = V_R/wm_R 

lambda_f_RMS = k_T*sqrt(2)/(p/2)/pi;

P_R = 500

R_ph = 12.4e-3;
L_ph = 154e-6;

NN = 0.1:0.1:3750;
wm = NN/60*2*pi;
we = wm*p/2;

X_ph = wm*L_ph;

Pm = zeros(size(wm));
for i = 1:1:length(NN)
    if NN(1,i)<3000
        Pm(1,i) = (NN(1,i)/3000)^3*5000;
    else
        Pm(1,i) = 5000;
    end
end

Pe_ph = -Pm;
Te_ph = Pe_ph./wm;

I_AC_rms = Te_ph*(2/(3*p*lambda_f_RMS));
Z_ph = (R_ph + j*X_ph);

E_rms = we*lambda_f_RMS;

Sout = 3*E_rms.*I_AC_rms;

%%%% E_rms I_rms
V_in = E_rms - I_AC_rms.*Z_ph;
V_in_rms = abs(V_in);

%%%%
I2R_ph = 3*I_AC_rms.^2*R_ph;
Pin = 3*V_in_rms.*I_AC_rms;
I_DC = Pin/V_R;


%%
close all
width = 600;
height = 200;
fig_i=1;

f(fig_i)= figure(fig_i);
plot(NN,Pm,'-k','LineWidth',1,'DisplayName','Mechanical Power Input')
title('Mechanical Power Input')
xlabel('Speed [RPM]')
ylabel('Power [W]')
ylim([0 6000])
grid on;
legend('boxoff')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

f(fig_i)= figure(fig_i);
plot(NN,Te_ph,'-k','LineWidth',1,'DisplayName','Torque')
title('Electromagnetic Torque')
xlabel('Speed [RPM]')
ylabel('Torque [N*M]')
grid on;
legend('boxoff')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

f(fig_i)= figure(fig_i);
plot(NN,I_AC_rms,'-k','LineWidth',1,'DisplayName','AC Current RMS')
title('AC Current RMS')
xlabel('Speed [RPM]')
ylabel('AC Current [A]')
grid on;
legend('boxoff','Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

figure(14)
plot(NN,E_rms)
title('E rms')

figure(15)
plot(NN,Sout)
title('Sout')

f(fig_i)= figure(fig_i);
plot(NN,V_in_rms,'-k','LineWidth',1,'DisplayName','AC Voltage RMS')
title('AC Voltage RMS')
xlabel('Speed [RPM]')
ylabel('AC Voltage [V]')
grid on;
legend('Location','southeast')
legend('boxoff')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

f(fig_i)= figure(fig_i);
plot(NN,I_DC,'-k','LineWidth',1,'DisplayName','DC Current')
title('DC Current')
xlabel('Speed [RPM]')
ylabel('DC Current[A]')
grid on;
legend('Location','northeast')
legend('boxoff')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

%%
%% Parameters for ECE 411 PRJ Sec2
% clear all
% clc

%% PMSM Machine Parameters
P   = 8;            % Number of poles
Rs   = 0.0124;       % Stator resistance per phase [Ohm]
Lq   = 154e-6;      % Stator q-axis inductance [H]
Ld   = 154e-6;      % Stator d-axis inductance [H]
lambda_f = lambda_f_RMS*sqrt(2);     % PMSM Rotor flux   [N-m/A]
RPM_max = 4259;     % Maximum Speed     [RPM]
Vstep = 1.0;        % Step voltage amplitude [V]

N_nl = 4259;
wm_nl = N_nl/60*2*pi;
we_nl = wm_nl *P/2;
E_rms_phase = we_nl*lambda_f_RMS
E_rms_ll_theoretical = sqrt(3)*E_rms_phase


%%
out_open = sim("ECE411_prj_OpenCircuitRotorSpinUp.slx");
%%
E_ll_rms_nl_simulation = max(out_open.V_LL.Data(:,1))/sqrt(2);

cum_theta = unwrap(out_open.theta_e.Data);

(cum_theta(end,1)-cum_theta(end-100,1))/(out_open.rpm.Time(end,1)-out_open.rpm.Time(end-100,1))

%%
close all
width = 600;
height = 200;
fig_i=20;
f(fig_i)= figure(fig_i);
plot(out_open.V_LL.Time,out_open.V_LL.Data(:,1),'-k','LineWidth',1,'DisplayName','V_{ab}')
hold on;
plot(out_open.V_LL.Time,out_open.V_LL.Data(:,2),'-b','LineWidth',1,'DisplayName','V_{bc}')
plot(out_open.V_LL.Time,out_open.V_LL.Data(:,3),'-r','LineWidth',1,'DisplayName','V_{ca}')
title('Line to Line voltages')
xlabel('time [t]')
ylabel('Voltage [V]')
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

f(fig_i)= figure(fig_i);
plot(out_open.rpm.Time,out_open.rpm.Data,'-k','LineWidth',1,'DisplayName','RPM')
title('RPM')
xlabel('time [t]')
ylabel('Rotor Velocity [RPM]')
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

f(fig_i)= figure(fig_i);
plot(out_open.theta_e.Time,out_open.theta_e.Data,'-k','LineWidth',1,'DisplayName','\theta_e')
title('theta_e')
xlabel('time [t]')
ylabel('Angle [rad]')
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

f(fig_i)= figure(fig_i);
plot(out_open.V_LL.Time,out_open.V_LL.Data(:,1),'-k','LineWidth',1,'DisplayName','V_{ab}')
hold on;
plot(out_open.V_LL.Time,out_open.V_LL.Data(:,2),'-b','LineWidth',1,'DisplayName','V_{bc}')
plot(out_open.V_LL.Time,out_open.V_LL.Data(:,3),'-r','LineWidth',1,'DisplayName','V_{ca}')
title('Line to Line voltages')
xlabel('time [t]')
ylabel('Voltage [V]')
xlim([1.95 2])
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

f(fig_i)= figure(fig_i);
plot(out_open.rpm.Time,out_open.rpm.Data,'-k','LineWidth',1,'DisplayName','RPM')
title('RPM')
xlabel('time [t]')
ylabel('Rotor Velocity [RPM]')
xlim([1.95 2])
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;


f(fig_i)= figure(fig_i);
plot(out_open.theta_e.Time,out_open.theta_e.Data,'-k','LineWidth',1,'DisplayName','\theta_e')
title('theta_e')
xlabel('time [t]')
ylabel('Angle [rad]')
xlim([1.95 2])
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

%% 2b2
va= 2/3;
vb= -1/3;
vc= -1/3;

vq = 2/3*(va-1/2*vb-1/2*vc);
vd = 2/3*(-sqrt(3)/2*vb+sqrt(3)/2*vc);

iqs = vq/Rs
ids = vd/Rs

ia = iqs;
ib = -1/2*iqs-sqrt(3)/2*ids;
ic = -1/2*iqs+sqrt(3)/2*ids;

torque_2b2 = 3/2*P/2*iqs*lambda_f_RMS*sqrt(2)
tau = Ld/Rs

ia*0.632
0.112136
%%


out_locked1 = sim("ECE411_prj_LockedRotorStepResponse11.slx")
%%
close all
width = 600;
height = 200;
fig_i=200;
f(fig_i)= figure(fig_i);
plot(out_locked1.Vin.Time,out_locked1.Vin.Data(:,1),'-k','LineWidth',1,'DisplayName','V_{in}')
hold on;
title('Step voltages')
xlabel('time [t]')
ylabel('Voltage [V]')
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

f(fig_i)= figure(fig_i);
plot(out_locked1.Iabc.Time,out_locked1.Iabc.Data(:,1),'-k','LineWidth',1,'DisplayName','I_{a}')
hold on;
title('A phase Current')
xlabel('time [t]')
ylabel('Current [A]')
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
% exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

f(fig_i)= figure(fig_i);
plot(out_locked1.torque.Time,out_locked1.torque.Data(:,1),'-k','LineWidth',1,'DisplayName','Torque')
hold on;
title('Torque')
xlabel('time [t]')
ylabel('Torque [Nm]')
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;
%%

%% 2b3
va= 0;
vb= 1/2;
vc= -1/2;

vq = 2/3*(va-1/2*vb-1/2*vc);
vd = 2/3*(-sqrt(3)/2*vb+sqrt(3)/2*vc);

iqs = vq/Rs
ids = vd/Rs

ia = iqs
ib = -1/2*iqs-sqrt(3)/2*ids
ic = -1/2*iqs+sqrt(3)/2*ids

torque_2b3 = 3/2*P/2*iqs*lambda_f_RMS*sqrt(2)
tau = Ld/Rs

40.3326*0.632
%%
out_locked2 = sim("ECE411_prj_LockedRotorStepResponse12.slx")
%%
close all
width = 600;
height = 200;
fig_i=2000;
f(fig_i)= figure(fig_i);
plot(out_locked2.Vin.Time,out_locked2.Vin.Data(:,1),'-k','LineWidth',1,'DisplayName','V_{in}')
hold on;
title('Step voltages')
xlabel('time [t]')
ylabel('Voltage [V]')
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

f(fig_i)= figure(fig_i);
plot(out_locked2.Iabc.Time,out_locked2.Iabc.Data(:,2),'-k','LineWidth',1,'DisplayName','I_{B}')
hold on;
title('B phase Current')
xlabel('time [t]')
ylabel('Current [A]')
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

f(fig_i)= figure(fig_i);
plot(out_locked2.torque.Time,out_locked2.torque.Data(:,1),'-k','LineWidth',1,'DisplayName','Torque')
hold on;
title('Torque')
xlabel('time [t]')
ylabel('Torque [Nm]')
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

%%
iq1 = squeeze(out_locked1.Iqd0.Data(1,1,:));
Torque1 = squeeze(out_locked1.torque.Data(:,1));
iq2 = squeeze(out_locked2.Iqd0.Data(1,1,:));
Torque2 = squeeze(out_locked2.torque.Data(:,1));

flux_hat1 = Torque1./iq1/(3*P/4);
flux_hat2 = Torque2./iq2/(3*P/4);
%% Rotor flux
fig_i = 7
f(fig_i)= figure(fig_i);
plot(out_locked1.Iqd0.Time,flux_hat1,'-r','LineWidth',1,'DisplayName','Estimate Lambda from 2a')
hold on;
% plot(out_locked2.Iqd0.Time,flux_hat2,'-b','LineWidth',1,'DisplayName','Torque')
title('Rotor Flux')
xlabel('time [t]')
ylabel('Flux [NM/A]')
grid on;
legend('Location','northeast')
ylim([0 0.05])
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

%% Sec 4
Jm = 0.2;       % Rotational inertia [Nms]

RPM_max = 3700;     % Maximum Speed     [RPM]
lambda_f_est = 1.0*lambda_f;

% Kp_ireg = 0.1416;  % qd current regulator proportional gain
% Ki_ireg = 141.6    % qd current regulator intergator gain
tau_ireg = 0.001;
fc_ireg = 1/tau_ireg;
wc_ireg = 2*pi*fc_ireg;
Kp_ireg = Lq*wc_ireg;
Ki_ireg = Rs*wc_ireg;

Vdc = 72;             % Nominal Vdc level
Kp_wreg = Jm*30;
Ki_wreg = Kp_wreg*(2*pi*1);

Ls = Lq;
t_end = 20;    % Use something like 2 seconds for Part 4a, longer for Part 4b, 

% maybe 20 sec for extra credit section
Te_cmd_max_pos = 10;  
Te_cmd_max_neg = -20;  % Limit of torque command Nm

Te_test_part4a = -15;
w_test_Part4a_RPM = 3000;  % 1000 RPM also

v_wind_ramp_m_per_s_per_s = 0.5;    % 16 sec to increase up to 8 m/s

%%
w_test_Part4a_RPM = 3000
out_4a = sim("ECE411_Project_Student_Starter_baik1.slx")
%%
out_4a_3000.Iabc = squeeze(out_4a.Iabc.Data)';
out_4a_3000.qd0 = squeeze(out_4a.qd0.Data);
out_4a_3000.Te = squeeze(out_4a.Te.Data)';
out_4a_3000.Time = squeeze(out_4a.Iabc.Time)';



%%
w_test_Part4a_RPM = 1000
out_4a = sim("ECE411_Project_Student_Starter_baik1.slx")
%%
out_4a_1000.Iabc = squeeze(out_4a.Iabc.Data)';
out_4a_1000.qd0 = squeeze(out_4a.qd0.Data);
out_4a_1000.Te = squeeze(out_4a.Te.Data)';
out_4a_1000.Time = squeeze(out_4a.Iabc.Time)';

%%
%%
close all;

fig_i = 40;
f(fig_i)= figure(fig_i);
plot(out_4a_3000.Time,out_4a_3000.Iabc(1,:),'-k','LineWidth',1,'DisplayName','I_a')
hold on;
plot(out_4a_3000.Time,out_4a_3000.Iabc(2,:),'-b','LineWidth',1,'DisplayName','I_b')
plot(out_4a_3000.Time,out_4a_3000.Iabc(3,:),'-r','LineWidth',1,'DisplayName','I_c')
title('I_{abc}')
xlabel('time [t]')
ylabel('Current [A]')
xlim([0.2 0.6])
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;


f(fig_i)= figure(fig_i);
plot(out_4a_3000.Time,out_4a_3000.qd0(1,:),'-k','LineWidth',1,'DisplayName','I_q')
hold on;
plot(out_4a_3000.Time,out_4a_3000.qd0(2,:),'-b','LineWidth',1,'DisplayName','I_d')
title('I_{qd}')
xlabel('time [t]')
ylabel('Current [A]')
xlim([0.2 0.6])
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;


f(fig_i)= figure(fig_i);
plot(out_4a_3000.Time,out_4a_3000.Te(1,:),'-k','LineWidth',1,'DisplayName','Torque')
hold on;
title('Torque')
xlabel('time [t]')
ylabel('Torque [Nm]')
xlim([0.2 0.6])
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

f(fig_i)= figure(fig_i);
plot(out_4a_1000.Time,out_4a_1000.Iabc(1,:),'-k','LineWidth',1,'DisplayName','I_a')
hold on;
plot(out_4a_1000.Time,out_4a_1000.Iabc(2,:),'-b','LineWidth',1,'DisplayName','I_b')
plot(out_4a_1000.Time,out_4a_1000.Iabc(3,:),'-r','LineWidth',1,'DisplayName','I_c')
title('I_{abc}')
xlabel('time [t]')
ylabel('Current [A]')
xlim([0.2 0.6])
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;


f(fig_i)= figure(fig_i);
plot(out_4a_1000.Time,out_4a_1000.qd0(1,:),'-k','LineWidth',1,'DisplayName','I_q')
hold on;
plot(out_4a_1000.Time,out_4a_1000.qd0(2,:),'-b','LineWidth',1,'DisplayName','I_d')
title('I_{qd}')
xlabel('time [t]')
ylabel('Current [A]')
xlim([0.2 0.6])
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

f(fig_i)= figure(fig_i);
plot(out_4a_1000.Time,out_4a_1000.Te(1,:),'-k','LineWidth',1,'DisplayName','Torque')
hold on;
title('Torque')
xlabel('time [t]')
ylabel('Torque [Nm]')
xlim([0.2 0.6])
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

%% Sec 4
Jm = 0.02;       % Rotational inertia [Nms]

RPM_max = 3700;     % Maximum Speed     [RPM]
lambda_f_est = 1.0*lambda_f;

% Kp_ireg = 0.1416;  % qd current regulator proportional gain
% Ki_ireg = 141.6    % qd current regulator intergator gain
tau_ireg = 0.001;
fc_ireg = 1/tau_ireg;
wc_ireg = 2*pi*fc_ireg;
Kp_ireg = Lq*wc_ireg;
Ki_ireg = Rs*wc_ireg;

Vdc = 72;             % Nominal Vdc level
Kp_wreg = Jm*30;
Ki_wreg = Kp_wreg*(2*pi*0.1);

Ls = Lq;
t_end = 20;    % Use something like 2 seconds for Part 4a, longer for Part 4b, 

% maybe 20 sec for extra credit section
Te_cmd_max_pos = 10;  
Te_cmd_max_neg = -20;  % Limit of torque command Nm

Te_test_part4a = -15;
w_test_Part4a_RPM = 3000;  % 1000 RPM also

v_wind_ramp_m_per_s_per_s = 0.5;    % 16 sec to increase up to 8 m/s

%% target power
Pm_500 = 500;
w_500 = (Pm_500/5000)^(1/3)*3000
v_wind_500 = w_500*8/3000

%%
t_end = 10;
v_wind = v_wind_500;
out_4b_500 = sim("ECE411_Project_Student_Starter_baik12.slx")

v_wind_5000 = 8;
v_wind = v_wind_5000;
out_4b_5000 = sim("ECE411_Project_Student_Starter_baik12.slx")
%%
out_4bb_500.Torque = squeeze(out_4b_500.Torque.Data)'; %%%% Te, Tecmd, Tm
out_4bb_500.speed = squeeze(out_4b_500.speed_rpm.Data)'; %%%% wr, wrcmd
out_4bb_500.power = squeeze(out_4b_500.Power.Data)'; %%%% Pm, Pdc
out_4bb_500.DCcurrent = squeeze(out_4b_500.Idc.Data)'; 
out_4bb_500.Time1 = squeeze(out_4b_500.Torque.Time)'; %%%% Torque, speed
out_4bb_500.Time2 = squeeze(out_4b_500.Power.Time)'; %%%% DC, power
%%
out_4bb_5000.Torque = squeeze(out_4b_5000.Torque.Data)'; %%%% Te, Tecmd, Tm
out_4bb_5000.speed = squeeze(out_4b_5000.speed_rpm.Data)'; %%%% wr, wrcmd
out_4bb_5000.power = squeeze(out_4b_5000.Power.Data)'; %%%% Pm, Pdc
out_4bb_5000.DCcurrent = squeeze(out_4b_5000.Idc.Data)';
out_4bb_5000.Time1 = squeeze(out_4b_5000.Torque.Time)'; %%%% Torque, speed
out_4bb_5000.Time2 = squeeze(out_4b_5000.Power.Time)';
%%
close all;
fig_i = 400;
width = 700;
height = 150;

f(fig_i)= figure(fig_i);
plot(out_4bb_500.Time1,out_4bb_500.Torque(1,:),'-k','LineWidth',1,'DisplayName','T_e')
hold on;
plot(out_4bb_500.Time1,out_4bb_500.Torque(2,:),'-r','LineWidth',1,'DisplayName','T_{ecmd}')
plot(out_4bb_500.Time1,out_4bb_500.Torque(3,:),'-b','LineWidth',1,'DisplayName','T_{wind}')
title('Torque Comparison')
xlabel('time [t]')
ylabel('Torque [Nm]')
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

f(fig_i)= figure(fig_i);
plot(out_4bb_500.Time1,out_4bb_500.speed(2,:),'-k','LineWidth',1,'DisplayName','Command rotor speed')
hold on;
plot(out_4bb_500.Time1,out_4bb_500.speed(1,:),'-b','LineWidth',1,'DisplayName','Actual rotor speed')
title('Speed Comparison')
xlabel('time [t]')
ylabel('Speed [RPM]')
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

f(fig_i)= figure(fig_i);
plot(out_4bb_500.Time2,out_4bb_500.power(1,:),'-k','LineWidth',1,'DisplayName','Mechanical Power')
hold on;
plot(out_4bb_500.Time2,out_4bb_500.power(2,:),'-b','LineWidth',1,'DisplayName','DC bus power')
title('Power Comparison')
xlabel('time [t]')
ylabel('Power [W]')
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

f(fig_i)= figure(fig_i);
plot(out_4bb_500.Time2,out_4bb_500.DCcurrent(1,:),'-k','LineWidth',1,'DisplayName','DC bus current')
title('DC bus current')
xlabel('time [t]')
ylabel('Current [A]')
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

%%%%% P =5000

f(fig_i)= figure(fig_i);
plot(out_4bb_5000.Time1,out_4bb_5000.Torque(1,:),'-k','LineWidth',1,'DisplayName','T_e')
hold on;
plot(out_4bb_5000.Time1,out_4bb_5000.Torque(2,:),'-r','LineWidth',1,'DisplayName','T_{ecmd}')
plot(out_4bb_5000.Time1,out_4bb_5000.Torque(3,:),'-b','LineWidth',1,'DisplayName','T_{wind}')
title('Torque Comparison')
xlabel('time [t]')
ylabel('Torque [Nm]')
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

f(fig_i)= figure(fig_i);
plot(out_4bb_5000.Time1,out_4bb_5000.speed(2,:),'-k','LineWidth',1,'DisplayName','Command rotor speed')
hold on;
plot(out_4bb_5000.Time1,out_4bb_5000.speed(1,:),'-b','LineWidth',1,'DisplayName','Actual rotor speed')
title('Speed Comparison')
xlabel('time [t]')
ylabel('Speed [RPM]')
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

f(fig_i)= figure(fig_i);
plot(out_4bb_5000.Time2,out_4bb_5000.power(1,:),'-k','LineWidth',1,'DisplayName','Mechanical Power')
hold on;
plot(out_4bb_5000.Time2,out_4bb_5000.power(2,:),'-b','LineWidth',1,'DisplayName','DC bus power')
title('Power Comparison')
xlabel('time [t]')
ylabel('Power [W]')
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;

f(fig_i)= figure(fig_i);
plot(out_4bb_5000.Time2,out_4bb_5000.DCcurrent(1,:),'-k','LinesWidth',1,'DisplayName','DC bus current')
title('DC bus current')
xlabel('time [t]')
ylabel('Current [A]')
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;
%%
t_end = 16;
out_optional = sim("ECE411_Project_Student_Starter_baik13.slx")

%%
out_opt.Torque = squeeze(out_optional.Torque.Data)'; %%%% Te, Tecmd, Tm
out_opt.speed = squeeze(out_optional.speed_rpm.Data)'; %%%% wr, wrcmd
out_opt.power = squeeze(out_optional.Power.Data)'; %%%% Pm, Pdc
out_opt.DCcurrent = squeeze(out_optional.Idc.Data)'; 

out_opt.VLL = squeeze(out_optional.VLL.Data)'; 
out_opt.Time1 = squeeze(out_optional.Torque.Time)'; %%%% Torque, speed
out_opt.Time2 = squeeze(out_optional.Power.Time)'; %%%% DC, power


%%
close all;
fig_i = 500;
width = 700;
height = 150;

f(fig_i)= figure(fig_i);
plot(out_4bb_500.speed,out_4bb_500.Torque(1,:),'-k','LineWidth',1,'DisplayName','T_e')
hold on;
title('Torque vs Speed')
xlabel('speed [RPM]')
ylabel('Torque [Nm]')
grid on;
legend('Location','northeast')
f(fig_i).Position = [100*fig_i 100*fig_i width height];
filename = sprintf('fig/fig%d.png',fig_i);
exportgraphics(f(fig_i),filename,'Resolution',1000);
fig_i = fig_i + 1;
