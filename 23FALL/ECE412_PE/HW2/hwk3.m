%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Hyeongmeen Baik
%%%% ECE412 Power Electronics
%%%% HW 3


%% Q3

clc;
clear;
Vo_min = 5
Vo_max = 20
IESR = 10
RESR = 0.05
Vin = 12

%%%% Vo_min case
D1 = (IESR*RESR+Vo_min)/(Vo_min+Vin)
Io1 = IESR-D1*IESR
Iin1 = D1*IESR
Pin1 = Vin*Iin1
Po1 = Io1*Vo_min
Rload1 = Vo_min/Io1
n1 = Po1/Pin1*100


%%%% Vo_max case
D2 = (IESR*RESR+Vo_max)/(Vo_max+Vin)
Io2 = IESR-D2*IESR
Iin2 = D2*IESR
Pin2 = Vin*Iin2
Po2 = Io2*Vo_max
Rload2 = Vo_max/Io2
n2 = Po2/Pin2*100

PESR = IESR^2*RESR


%% Q4
%%%% Vo_min case
Vs1 = (1-D1)*(Vin+Vo_min)
VsRMS1 = sqrt(1-D1)*(Vin+Vo_min)
VsPeak1 = Vin+Vo_min
Vd1 = D1*(Vo_min+Vin)
VdRMS1 = sqrt(D1)*(Vo_min+Vin)
Vdpeak1 = Vo_min+Vin

%%%% Vo_max case
Vs2 = (1-D2)*(Vin+Vo_max)
VsRMS2 = sqrt(1-D2)*(Vin+Vo_max)
VsPeak2 = Vin+Vo_max
Vd2 = D2*(Vo_max+Vin)
VdRMS2 = sqrt(D2)*(Vo_max+Vin)
Vdpeak2 = Vo_max+Vin

%% Q5
clc;
clear;
Vin = 12;
Vo = 20;
D_boost = (Vo-Vin)/Vo
D_buckboost = Vo/(Vin+Vo)

deliL = 100e-3
fs = 50e3
Ts = 1/fs

L_boost = Vin*D_boost*Ts/deliL
L_buckboost = Vin*D_buckboost*Ts/deliL

%%%% ideal case
Vsbuckboost = (1-D_buckboost)*(Vin+Vo)
VsbuckboostRMS = sqrt(1-D_buckboost)*(Vin+Vo)
VsbuckboostPeak = (Vin+Vo)
Vdbuckboost = D_buckboost*(Vo+Vin)
VdbuckboostRMS = sqrt(D_buckboost)*(Vo+Vin)
VdbuckboostPeak = Vo+Vin

%%%% ideal case
Vsboost = (1-D_boost)*(Vo)
VsboostRMS = sqrt(1-D_boost)*(Vo)
VsboostPeak = (Vo)
Vdboost = D_buckboost*(Vo-Vin)
VdboostRMS = sqrt(D_buckboost)*(Vo-Vin)
VdboostPeak = Vo-Vin

%%%% Considering ESR for boost converter

Vin = 12
Vo= 20
RESR = 0.05

syms Iin
eqn=Vin-Iin*RESR == Vo/Iin*20/5.5652
S=solve(eqn)
Iin = double(S(1,1))
D = 1-3.5938/Iin
