%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Power Electronics ECE 412 Homework #2
%%%%%%% Writer: Hyeongmeen Baik
%%%%%%% Email: hyeongmeen.baik@wisc.edu
%%%%%%% Link: github.com/PhilBaik/WISC_2023Fall/blob/master/ECE412_PE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;

VF= 1
Vin = 100
D = 0.75
Vc = D*Vin - (1-D)*VF
Io = 10

tcr = 20e-9;
tvf = 30e-9;
tvr = 35e-9;
tcf = 15e-9;

Woff = 1/2*(Vin+VF)*Io*(tcf+tvr)
Won = 1/2*(Vin+VF)*Io*(tcr+tvf)
Rds_on = 0.18

Is = D*Io
Is_RMS = sqrt(D)*Io
is_peak = Io
Ps_cond = Rds_on*D*Io^2

Id = Io-D*Io
Id_RMS = sqrt(1-D)*Io
id_peak = Io
Pd_cond = VF*Io*(1-D)

%% Q3
fs_low = 1e3;
fs_high = 500e3;

Ps_swit1 = fs_low*(Woff+Won)
Ps_swit2 = fs_high*(Woff+Won)
Ps_tot1 = Ps_swit1 + Ps_cond
Ps_tot2 = Ps_swit2 + Ps_cond

RthJA = 40;
RthJC = 0.83;
RthCS = 0.24;
Rthcase = 4;
RthJA_heat_sink = RthJC+RthCS+Rthcase;


Tdiff1 = RthJA*Ps_tot1
Tdiff2 = RthJA*Ps_tot2

Tdiff1_case = RthJA_heat_sink*Ps_tot1
Tdiff2_case = RthJA_heat_sink*Ps_tot2

P_tot1 = Ps_tot1 + Pd_cond
P_tot2 = Ps_tot2 + Pd_cond

Pin = Vin*Is

Po1 = Pin-P_tot1
Po2 = Pin-P_tot2

n1 = Po1/Pin
n2 = Po2/Pin

Vo1 = Po1/Io
Vo2 = Po2/Io

%% q5
temp_diff_max = 125
Ps_max = temp_diff_max/RthJA_heat_sink
Ps_switch_max = Ps_max - Ps_cond
fs_max = Ps_switch_max/(Woff+Won)
