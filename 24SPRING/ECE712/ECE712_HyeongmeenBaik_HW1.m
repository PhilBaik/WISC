clc
clear

f = 60;
L = 600e-6;
VDC = 800;

Q = [-8000 0 8000]
P = [-10000 0 10000]

AC_grid_voltage_peak = 480/sqrt(3)*sqrt(2)
Vq = AC_grid_voltage_peak

iq = zeros(3,3)
id = zeros(3,3)
for l = 1:1:3
    for m = 1:1:3
        iq(l,m) = P(1,l)/3*2/Vq;
        id(l,m) = Q(1,m)/3*2/Vq;
    end
end
iqd = iq-1j*id
mag_iqd = abs(iqd)
ang_iqd = angle(iqd)/pi*180


max_rms_iqd = max(max(mag_iqd))/sqrt(2)
mqd = 2/VDC*(j*2*pi*f*L*(iqd)+Vq)

%% prob 5~7
P = 10e3;
Q = 8e3;
Iqd = iq(3,3)-j*id(3,3)
Mqd = mqd(3,3)
abs_Iqd = abs(Iqd)
a = exp(j*2*pi/3);

Ts = 0.000001
T_end = 1
t = 0:Ts:T_end-Ts;

Voa = Vq*cos(2*pi*f*t);
Vob = Vq*cos(2*pi*f*t-2*pi/3);
Voc = Vq*cos(2*pi*f*t-4*pi/3);

Voabc = [Voa;Vob;Voc];

iqd_extended = zeros(3,length(t));
iqd_extended(1,:) = real(Iqd);
iqd_extended(2,:) = -imag(Iqd);
iqd_extended(3,:) = 0; %%zero sequence

iabc = qd2abc(2*pi*f*t,iqd_extended);
Ipa = iabc(1,:);
Ipb = iabc(2,:);
Ipc = iabc(3,:);



% RMS verification
S_est = 3*rms(Ipa)*rms(Voa);
S = sqrt(P^2+Q^2);
P = sum(Voa.*Ipa)/T_end*3; % this should be 10e3

% Vpa with three winding neutral point
Vpqd = j*2*pi*f*L*Iqd + Vq
Vpqd_extended = zeros(3,length(t));
Vpqd_extended(1,:) = real(Vpqd);
Vpqd_extended(2,:) = -imag(Vpqd);
Vpqd_extended(3,:) = 0;
Vpabc = qd2abc(2*pi*f*t,Vpqd_extended);
Vpa = Vpabc(1,:);
Vpb = Vpabc(2,:);
Vpc = Vpabc(3,:);

% with ground reference
Vpa_ground = Vpa+400;
Vpb_ground = Vpb+400;
Vpc_ground = Vpc+400;

% Mqd
Mqd_extended = zeros(3,length(t));
Mqd_extended(1,:) = real(Mqd);
Mqd_extended(2,:) = -imag(Mqd);
Mqd_extended(3,:) = 0;
Mabc = qd2abc(2*pi*f*t,Mqd_extended);
Ma = Mabc(1,:);
Mb = Mabc(2,:);
Mc = Mabc(3,:);

da1 = zeros(size(Ma));
da2 = zeros(size(Ma));
da3 = zeros(size(Ma));
db1 = zeros(size(Ma));
db2 = zeros(size(Ma));
db3 = zeros(size(Ma));
dc1 = zeros(size(Ma));
dc2 = zeros(size(Ma));
dc3 = zeros(size(Ma));
% Sawtooth
% carrier1 = sawtooth(2*pi*f*t,1/2)/2+1/2;
% carrier2 = sawtooth(2*pi*f*t,1/2)/2-1/2;

% duty ratio
for i = 1:1:length(t)
    if Ma(1,i)>0
        da1(1,i) = Ma(1,i);
        da2(1,i) = 1-da1(1,i);
        da3(1,i) = 0;
    else
        da3(1,i) = -Ma(1,i);
        da2(1,i) = 1-da3(1,i);
        da1(1,i) = 0;
    end

    if Mb(1,i)>0
        db1(1,i) = Mb(1,i);
        db2(1,i) = 1-db1(1,i);
        db3(1,i) = 0;
    else
        db3(1,i) = -Mb(1,i);
        db2(1,i) = 1-db3(1,i);
        db1(1,i) = 0;
    end

    if Mc(1,i)>0
        dc1(1,i) = Mc(1,i);
        dc2(1,i) = 1-dc1(1,i);
        dc3(1,i) = 0;
    else
        dc3(1,i) = -Mc(1,i);
        dc2(1,i) = 1-dc3(1,i);
        dc1(1,i) = 0;
    end
end

%% Prob 6
ita1_average = sum(da1.*Ipa)*Ts/T_end
itb1_average = sum(db1.*Ipb)*Ts/T_end
itc1_average = sum(dc1.*Ipc)*Ts/T_end
it01 = ita1_average+itb1_average+itc1_average

ita2_average = sum(da2.*Ipa)*Ts/T_end
itb2_average = sum(db2.*Ipb)*Ts/T_end
itc2_average = sum(dc2.*Ipc)*Ts/T_end
it02 = ita2_average+itb2_average+itc2_average

ita3_average = sum(da3.*Ipa)*Ts/T_end;
itb3_average = sum(db3.*Ipb)*Ts/T_end;
itc3_average = sum(dc3.*Ipc)*Ts/T_end;

itopsum = da1.*Ipa+db1.*Ipb+dc1.*Ipc;
itop_average = sum(itopsum)*Ts/T_end
%% plot
close all;

fig_num=51;
fig=figure(fig_num);
fig.Position = [10 10 800 300];
plot(t,Voa,'DisplayName','Voa'); hold on;
title("V_{oa}")
grid on; legend; xlim([0 2/f])
fig_num= fig_num+1;

fig = figure(fig_num);
fig.Position = [10+(fig_num-51)*5 10+(fig_num-51)*5 800 300];
plot(t,Ipa,'DisplayName','Ipa'); hold on;
plot(t,Ipb,'DisplayName','Ipb');
plot(t,Ipc,'DisplayName','Ipc');
title("I_{pa}")
grid on; legend; xlim([0 2/f])
fig_num= fig_num+1;

fig = figure(fig_num);
fig.Position = [10+(fig_num-51)*5 10+(fig_num-51)*5 800 300];
plot(t,Vpa,'DisplayName','Vpa from neutral'); hold on;
plot(t,Voa,'DisplayName','Voa');;
title("V_{pa}")
grid on; legend; xlim([0 2/f])
fig_num= fig_num+1;

fig = figure(fig_num);
fig.Position = [10+(fig_num-51)*5 10+(fig_num-51)*5 800 300];
plot(t,Vpa_ground,'DisplayName','Vpa from ground'); hold on;
plot(t,Vpb_ground,'DisplayName','Vpb from ground');plot(t,Vpc_ground,'DisplayName','Vpc from ground');
title("V_{pa} ground reference")
grid on; legend; xlim([0 2/f])
fig_num= fig_num+1;

fig = figure(fig_num);
fig.Position = [10+(fig_num-51)*5 10+(fig_num-51)*5 800 300];
plot(t,Ma,'DisplayName','m_a'); hold on;
title("m_a, modulation function")
plot(t,Mb,'DisplayName','m_b');
plot(t,Mc,'DisplayName','m_c');
grid on; legend; xlim([0 2/f])
fig_num= fig_num+1;

fig = figure(fig_num);
fig.Position = [10+(fig_num-51)*5 10+(fig_num-51)*5 800 300];
plot(t,da1,'DisplayName','da1'); hold on;
title("d_{a123}")
plot(t,da2,'DisplayName','da2');
plot(t,da3,'DisplayName','da3');
grid on; legend; xlim([0 2/f])
fig_num= fig_num+1;

fig = figure(fig_num);
fig.Position = [10+(fig_num-51)*5 10+(fig_num-51)*5 800 300];
plot(t,db1,'DisplayName','db1'); hold on;
title("d_{b123}")
plot(t,db2,'DisplayName','db2');
plot(t,db3,'DisplayName','db3');
grid on; legend; xlim([0 2/f])
fig_num= fig_num+1;

fig = figure(fig_num);
fig.Position = [10+(fig_num-51)*5 10+(fig_num-51)*5 800 300];
plot(t,dc1,'DisplayName','dc1'); hold on;
title("d_{c123}")
plot(t,dc2,'DisplayName','dc2');
plot(t,dc3,'DisplayName','dc3');
grid on; legend; xlim([0 2/f])
fig_num= fig_num+1;

fig = figure(fig_num);
fig.Position = [10+(fig_num-51)*5 10+(fig_num-51)*5 800 300];
plot(t,da1,'DisplayName','da1'); hold on;
title("d_{abc1}")
plot(t,db1,'DisplayName','db1');
plot(t,dc1,'DisplayName','dc1');
grid on; legend; xlim([0 2/f])
fig_num= fig_num+1;

fig = figure(fig_num);
fig.Position = [10+(fig_num-51)*5 10+(fig_num-51)*5 800 300];
plot(t,da2,'DisplayName','da2'); hold on;
title("d_{abc2}")
plot(t,db2,'DisplayName','db2');
plot(t,dc2,'DisplayName','dc2');
grid on; legend; xlim([0 2/f])
fig_num= fig_num+1;

fig = figure(fig_num);
fig.Position = [10+(fig_num-51)*5 10+(fig_num-51)*5 800 300];
plot(t,da3,'DisplayName','da3'); hold on;
title("d_{abc3}")
plot(t,db3,'DisplayName','db3');
plot(t,dc3,'DisplayName','dc3');
grid on; legend; xlim([0 2/f])
fig_num= fig_num+1;


fig = figure(fig_num);
fig.Position = [10+(fig_num-51)*5 10+(fig_num-51)*5 800 300];
plot(t,itopsum,'DisplayName','i top throw'); hold on;
title("Top throw current")
grid on; legend; xlim([0 2/f])
fig_num= fig_num+1;
%% function
function output = qd2abc(theta,qd0)
    a = qd0(1,:).*cos(theta)+qd0(2,:).*sin(theta)+qd0(3,:) +qd0(3,:);
    b = cos(theta).*(-sqrt(3)/2*qd0(2,:)-1/2*qd0(1,:))+sin(theta).*(sqrt(3)/2*qd0(1,:)-1/2*qd0(2,:))+qd0(3,:);
    c = cos(theta).*(sqrt(3)/2*qd0(2,:)-1/2*qd0(1,:))+sin(theta).*(-sqrt(3)/2*qd0(1,:)-1/2*qd0(2,:))+qd0(3,:);
    output = [a;b;c];
end