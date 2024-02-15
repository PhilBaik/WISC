%%%%% Hyeongmeen Baik HW3-1

%% Q3 )a)
clc
clear

s = 0.02

r1 = 0.11;
r2 = 0.15;
L1 = 1.6e-3;
L2 = 1.7e-3;
Lm = 38.2e-3;
Vll = 480
Vs = Vll/sqrt(3)

fe = 60;
p =8;

we = fe*2*pi;
wrm = we*(1-s)*2/p;

X1 = L1*we
X2 = L2*we
Xm = Lm*we

Zin = r1+j*X1 + j*Xm*(j*X2 + r2/s)/(j*Xm+j*X2+r2/s);

I1 = Vs/Zin
I1_abs = abs(I1)
I1_ph = angle(I1)

I2 = (Vs-I1*(r1+j*X1))/(j*X2+r2/s)
I2_abs = abs(I2)
I2_ph = angle(I2)

Po_mech = 3*r2*(1-s)/s*I2_abs^2
Te = Po_mech/wrm

pf = cos(I1_ph)

I2R1 = 3*I1_abs^2*r1
I2R2 = 3*I2_abs^2*r2

Pin = 3*real(Vs*conj(I1))

n = Po_mech/Pin*100

%%
s1 = 1
s2 = -1

% for s1 = 1

ss = -1:0.001:1;
Zin = r1+j*X1 + j*Xm*(j*X2 + r2./ss)./(j*Xm+j*X2+r2./ss);

wrm = we*(1-ss)/p*2;
I1 = Vs./Zin;
I2 = (Vs-I1*(r1+j*X1))./(j*X2+r2./ss);
I2_abs = abs(I2);

Po_mech = 3*r2*(1-ss)./ss.*I2_abs.^2;
Te = Po_mech./wrm;

figure(11)
plot(wrm,Te)

%% Q3-b
PB = 55.9275e3;
wB = 2*pi*60;
wmB = 2*wB/p;
VB = Vs ;
IB = PB/3/VB;
ZB = VB/IB;
TB = PB/wmB;

r1_pu = r1/ZB;
r2_pu = r2/ZB;
X1_pu = X1/ZB;
X2_pu = X2/ZB;
Xm_pu = Xm/ZB;
Rm = 180/ZB;

%% Q3-b
PB = 55.9275e3;
wB = 2*pi*60;
wmB = 2*wB/p;
VB = Vs ;
IB = PB/3/VB;
ZB = VB/IB;
TB = PB/wmB;

r1_pu = r1/ZB;
r2_pu = r2/ZB;
X1_pu = X1/ZB;
X2_pu = X2/ZB;
Xm_pu = Xm/ZB;
Rm = 180/ZB;

%% Q4
clc
clear
r1 = 0.02
r2 = 0.026
X1 = 0.08
X2 = 0.08
Xm = 4.1
rm = 30
Vs=1

% s = 1
s =1
Zm = j*Xm*rm/(j*Xm+rm);
Zin = r1+j*X1 + Zm*(j*X2+r2/s)/(Zm+j*X2+r2/s)
I1 = 1/Zin
I1_abs = abs(I1)
I2 = (Vs-I1*(r1+j*X1))/(j*X2+r2/s)
I2_abs = abs(I2)

Te = abs(I2)^2*r2/s/1

Is_hat = Vs/(X1 + X2)
Ts_hat = r2*Vs^2/(X1+X2)^2

%% q4b)
wr = 1.032;
we = 1;
s = 1-wr;
Te = -1.07;
Ir = sqrt(Te*s*we/r2);

Vx = Ir*(j*X2+r2/s);
Im = Vx/Zm;

Is = Ir+Im
Is_abs = abs(Is)
Is_ph = angle(Is)

Vs = Is*(r1+j*X1)+Vx
Vs_abs = abs(Vs)

P1 = real(Vs*conj(Is))
P2 = Ir^2*r2*(1-s)/s

n = P1/P2*100
