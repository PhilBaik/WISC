clc
clear
close all

vp_hat = 110*sqrt(2);
fs = 60;
Ts = 1/fs;
Po = 13;
Vo = 26;
Io = Po/Vo;
w = 2*pi*fs;
vd = 1;
vdce = 2*vd+Vo;

%%%%% Q1
N_max = 28/110/sqrt(2)/0.9;
N_min = 28/110/sqrt(2);
N = N_max;

vs_hat = N*vp_hat;

te = 1/w*sqrt(2*(1-28/vs_hat));
Rs = 4*te*(vs_hat-Vo-2*vd)/Ts/Io;
Io = 4*te*(vs_hat - vdce)/Rs/Ts;
C_min = Io*(Ts/2-2*te)/0.05/Vo;

R = 0.9;
Rs_2 = 4/sqrt(2)/pi*(vdce)/R/Io*(1-R)^(3/2)
%%%% Q2
%%%%%% vs
Tss = 1e-4;
tt = 0:Tss:Ts;
vs = vs_hat*cos(w*tt);
vdce_plot = (2*vd + Vo)*ones(size(vs));

%%%%%% is
is_high = (vs-vdce)/Rs;
is_high(is_high<0) = 0; 
is_low = (vs+vdce)/Rs;
is_low(is_low>0) = 0; 

is = is_high + is_low;
is_rect_high = is_high;
is_rect_high(is_rect_high>0) = (vs_hat-vdce)/Rs;
is_rect_low = is_low;
is_rect_low(is_rect_low<0) = -(vs_hat-vdce)/Rs;
is_rect = is_rect_high + is_rect_low;

%%%%%% id1
id = abs(is);
id_rect = abs(is_rect);

id1 = is;
id1_rect = is_rect;
id1(id1<0) = 0;
id1_rect(id1_rect<0) = 0;

%%%%%% vd1
vd1 = 13+vs/2-26;
vd1(is>0)=1;
vd1(is<0)=-27;

%%%%%% Capacitor
ic_rect = id_rect - Io;

%%%%%% Load
io = 1/2*ones(size(tt));

%% Q3
%%%%% a)
is_peak = max(is_rect);
Is_RMS = sqrt(2/Ts*(is_peak^2)*2*te);

Is1 = 2*sqrt(2)*is_peak/pi*sin(2*pi/Ts*te);
THD_F = sqrt(Is_RMS^2-Is1^2)/Is1;

%%%%% d)
Vs_RMS = vs_hat/sqrt(2)
Ss = Vs_RMS*Is_RMS;
Ps = 4*vs_hat*is_peak*sin(w*te)/2/pi;
PFs = Ps/Ss

%%%%% e)
ic_ppeak = max(ic_rect);
ic_npeak = min(ic_rect);
Ic_RMS = sqrt(1/(Ts/2)*(2*te*ic_ppeak^2-(Ts/2-2*te)*ic_npeak^2))

%%%%% f)
id_peak = max(id_rect)
id_npeak = min(id_rect);
Id_ave = 2/Ts*te*id_peak
%% plotting
close all;
figure_size = [500 300]

figure1 = figure(1);
figure1.Position(3:4) = figure_size;
plot(tt,vs,'DisplayName','v_s')
title('v_s vs time')
ylabel('voltage [V]')
xlabel('time [s]')
xlim([0 Ts])
grid on;
hold on;
plot(tt,vdce_plot,'DisplayName','v_{dce}')
plot(tt,-vdce_plot,'DisplayName','-v_{dce}')
scatter(te,vdce,'black','DisplayName','t_e')
legend('Location','best')

figure2=figure(2);
figure2.Position(3:4) = figure_size;
plot(tt,is,'DisplayName','i_s')
title('i_s vs time')
ylabel('current [A]')
xlabel('time [s]')
xlim([0 Ts])
grid on;
hold on;
plot(tt,is_rect,'DisplayName','i{s.approximation}')
legend('Location','best')

figure3=figure(3);
figure3.Position(3:4) = figure_size;
yyaxis left;
plot(tt,id1,'DisplayName','i_{d1}')
hold on;
plot(tt,id1_rect,'DisplayName','i_{d1.approx}')
ylabel('current [A]')
title('diode vs time')

yyaxis right;
ylabel('voltage [V]')
plot(tt,vd1,'DisplayName','v_{d1}')
xlabel('time [s]')
xlim([0 Ts])
grid on;
legend('Location','best')

figure4=figure(4);
figure4.Position(3:4) = figure_size;
plot(tt,ic_rect,'DisplayName','i_{C}')
yyaxis left;
title('Capacitor vs time')
ylabel('Current [A]')
hold on;
xlabel('time [s]')
xlim([0 Ts])
yyaxis right;
ylabel('Voltage [V]')
ylim([20 30])
grid on;
legend('Location','best')

figure5=figure(5);
figure5.Position(3:4) = figure_size;
plot(tt,io,'DisplayName','i_{o}')
yyaxis left;
title('Load vs time')
ylabel('Current [A]')
hold on;
xlabel('time [s]')
xlim([0 Ts])
yyaxis right;
ylabel('Voltage [V]')
ylim([20 30])
grid on;
legend('Location','best')
