%%%%%%%%%%%% ECE412 HW 5
%%%%%%%%%%%% Hyeongmeen Baik

clc
clear

Lm = 150e-6
d = 0.4
Fs = 75e3
Ts = 1/Fs
Vin = 100
Vo = 10
n = 15/100
Io = 20


syms M(D)
M(D) = D*n/(1-D)
subs(M(D),D,0.4)

del_im = Vin*d*Ts/Lm
del_is = del_im/n

is_min = Io/(1-d)-del_is/2
is_mid = is_min + del_is/2
is_max = is_min + del_is

vp_low  = -1*Vo/n

im_max = n*(is_min+del_is)
im_min = n*is_min
Im = (im_max+im_min)/2
Ip = Im*d
Im/n*(1-d)



%% Q3 a
clc
clear

T1 = 10e-3
T2 = 0.1333
fp1 = 1/T1
fp2 = 1/T2
Tsample = 1e-7
VDC = 120;
Tend = 1

tt = 0:Tsample:Tend-Tsample;
m1 = 0.6
d1a = 1/2 + m1/2*cos(2*pi*fp1*tt);
d2a = 1/2 - m1/2*cos(2*pi*fp1*tt);

T_tri = 0.02e-3;
f_tri = 1/T_tri;
tri = (sawtooth(tt*2*pi*f_tri,1/2)+1)/2;

difference = d1a-tri;
difference_index = difference>0;
difference(difference_index) = 1;
difference_index = difference<=0;
difference(difference_index) = 0;
switch_signal = difference;

difference2 = d2a-tri;
difference_index = difference2>0;
difference2(difference_index) = 1;
difference_index = difference2<=0;
difference2(difference_index) = 0;
switch_signal2 = difference2;

vp1a = switch_signal*VDC;
vp2a = switch_signal2*VDC;
vla = vp1a-vp2a;
%% plotting q3a
close all;

figure(31)
subplot(2,1,1)
plot(tt,d1a,'DisplayName','duty 1')
hold on;
plot(tt,tri,'DisplayName','tri')
xlim([0 T_tri])
grid on
legend

subplot(2,1,2)
plot(tt,switch_signal,'DisplayName','Switching Signal')
xlim([0 T_tri])
grid on
legend

figure(32)
subplot(2,1,1)
plot(tt,d2a,'DisplayName','duty 2')
hold on;
plot(tt,tri,'DisplayName','tri')
xlim([0 T_tri])
grid on
legend

subplot(2,1,2)
plot(tt,switch_signal2,'DisplayName','Switching Signal2')
xlim([0 T_tri])
grid on
legend

figure(33)
plot(tt,vla,'DisplayName','vl')
hold on;
xlim([0 T_tri])
grid on
legend

%% Q3 b

m2 = 0.2
d1b = 1/2 + m2/2*cos(2*pi*fp1*tt);
d2b = 1/2 - m2/2*cos(2*pi*fp1*tt);

T_tri = 0.02e-3;
f_tri = 1/T_tri;
tri = (sawtooth(tt*2*pi*f_tri,1/2)+1)/2;

difference = d1b-tri;
difference_index = difference>0;
difference(difference_index) = 1;
difference_index = difference<=0;
difference(difference_index) = 0;
switch_signal = difference;

difference2 = d2b-tri;
difference_index = difference2>0;
difference2(difference_index) = 1;
difference_index = difference2<=0;
difference2(difference_index) = 0;
switch_signal2 = difference2;

vp1b = switch_signal*VDC;
vp2b = switch_signal2*VDC;
vlb = vp1b-vp2b;
%% plotting q3b
close all;

figure(331)
subplot(2,1,1)
plot(tt,d1b,'DisplayName','duty 1')
hold on;
plot(tt,tri,'DisplayName','tri')
xlim([0 T_tri])
grid on
legend

subplot(2,1,2)
plot(tt,switch_signal,'DisplayName','Switching Signal')
xlim([0 T_tri])
grid on
legend

figure(332)
subplot(2,1,1)
plot(tt,d2b,'DisplayName','duty 2')
hold on;
plot(tt,tri,'DisplayName','tri')
xlim([0 T_tri])
grid on
legend

subplot(2,1,2)
plot(tt,switch_signal2,'DisplayName','Switching Signal2')
xlim([0 T_tri])
grid on
legend

figure(333)
plot(tt,vlb,'DisplayName','vl')
hold on;
xlim([0 T_tri])
grid on
legend

%% Q5 a)
fp1 = 100;
V1 = 72
Vl = V1
Zload = 0.5 + j*0.1*(1e-3)*2*pi*fp1
Il = Vl/Zload;
Il_abs = abs(Il)
Il_phase = angle(Il)
Il_angle = Il_phase*180/pi;

%% Q5 b)
fp1 = 7.5019;
V1 = 24
Vl = V1
Zload = 0.5 + j*0.1*(1e-3)*2*pi*fp1
Il = Vl/Zload;
Il_abs = abs(Il)
Il_phase = angle(Il)
Il_angle = Il_phase*180/pi

%% Verification 
n = length(tt);
ff = (0:1:n-1)/Tsample/n;
fft_vl = fft(vla);
ff_shift = (-n/2:n/2-1)*(1/Tsample/n); 
fft_vl_shift = fftshift(fft_vl);

figure(51)
plot(ff_shift,abs(fft_vl_shift));

max(abs(fft_vl_shift))*2/n


