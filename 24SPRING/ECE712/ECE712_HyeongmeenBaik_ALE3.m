clc
clear
close all
ts = 0.000001;
t = 0:ts:1-ts;
w = 2*pi*60;
Ts = 1/1000;

ma = 0.9*cos(w*t);
da1 = 1/2*(ma+1);
da2 = 1-da1;

carrier = abs(2/Ts*mod(t,Ts)-1);
ha1= (sign(da1-carrier)+1)/2;

%% a)
figure(11)
plot(t,carrier,'DisplayName','Carrier');hold on;grid on;legend
plot(t,da1,'DisplayName','da1');
plot(t,ha1,'LineWidth',2,'DisplayName','ha1');
xlim([0 1/60])

%% b)
vpa = ha1*800;
vpa_average = da1*800;

figure(2)
plot(t,vpa,'DisplayName','vpa(t)');hold on;grid on;legend
plot(t,vpa_average,'DisplayName','<vpa(t)>');
xlim([0 1/60])

%% c)
ipa_t = 30*cos(w*t-pi/3);
it_instant = ha1.*ipa_t;
it_average = da1.*ipa_t;
it_rms = sqrt(da1).*ipa_t;

figure(3)
plot(t,ipa_t,'DisplayName','i_{pa}(t)'); hold on; grid on; legend
plot(t,it_instant,'DisplayName','instant. throw current');
plot(t,it_average,'DisplayName','average throw current');
plot(t,it_rms,'DisplayName','rms throw current');
xlim([0 1/60])

%% d)
% igbt 
igbt_instant = it_instant;
igbt_average = it_average;
igbt_rms = it_rms;
igbt_instant(find(it_instant<=0)) = 0;
igbt_average(find(it_average<=0)) = 0;
igbt_rms (find(it_rms<=0)) = 0;

diode_instant = -it_instant;
diode_average = -it_average;
diode_rms = -it_rms;
diode_instant(find(-it_instant<=0)) = 0;
diode_average(find(-it_average<=0)) = 0;
diode_rms (find(-it_rms<=0)) = 0;

%%
figure(4)
plot(t,igbt_instant,'DisplayName','i_{igbt} instantaneous (t)'); hold on; grid on; legend
plot(t,igbt_average,'DisplayName','average igbt current');
plot(t,igbt_rms,'DisplayName','rms igbt throw current');
xlim([0 1/60])

figure(5)
plot(t,diode_instant,'DisplayName','i_{diode} instantaneous (t)'); hold on; grid on; legend
plot(t,diode_average,'DisplayName','average diode current');
plot(t,diode_rms,'DisplayName','rms diode throw current');
xlim([0 1/60])

%% e)

