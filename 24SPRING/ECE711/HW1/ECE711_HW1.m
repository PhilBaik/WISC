%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Author: Hyeongmeen Baik, Dheeraj
%%%%%%%%%%%% Date: 02/07/2024
%%%%%%%%%%%% Title: ECE 711 - HW1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;

%% Making Functions

%%% a) 
No = 30;
theta_step = 0.0001
theta = 0:theta_step:2*pi-theta_step;

%%% finding index for winding function
temp_index1 = find(theta>pi/2,1);
temp_index2 = find(theta>3*pi/2,1);

N1_a = zeros(size(theta));
N1_a(1,1:end) = No/2;
N1_a(1,temp_index1:temp_index2-1) = -No/2;

temp_index1 = find(theta>pi,1);
temp_index2 = find(theta>3*pi/2,1);
N2_b = zeros(size(theta));
for i = 1:1:temp_index1-1
    N2_b(1,i) = (temp_index1/2-i+1)/temp_index1;
end
for i = temp_index1:1:length(theta)
    N2_b(1,i) = i/temp_index1-3/2;
end

%%% c) 

temp_index1 = find(theta>pi,1);

N1_c = zeros(size(theta));
N1_c(1,1:end) = 30;
N1_c(1,1:temp_index1-1) = -30;

temp_index1 = find(theta>pi/6,1);
temp_index2 = find(theta>4*pi/6,1);
temp_index3 = find(theta>5*pi/6,1);
temp_index4 = find(theta>6*pi/6,1);
temp_index5 = find(theta>7*pi/6,1);
temp_index6 = find(theta>10*pi/6,1);
temp_index7 = find(theta>11*pi/6,1);

N2_c = N1_c;
N2_c(1,1:temp_index1-1) = -20;
N2_c(1,temp_index2:temp_index3-1) = -20;
N2_c(1,temp_index3:temp_index4-1) = 0;
N2_c(1,temp_index4:temp_index5-1) = 20;
N2_c(1,temp_index5:temp_index6-1) = 30;
N2_c(1,temp_index6:temp_index7-1) = 20;
N2_c(1,temp_index7:end) = 0;

%% For FFT, longer version of the signal
iter_length = 100;
theta_long = unwrap(repmat(theta,1,iter_length));
N1_c = repmat(N1_c,1,iter_length);
N2_c = repmat(N2_c,1,iter_length);
N2_bb = repmat(N2_b,1,iter_length);
N1_aa = repmat(N1_a,1,iter_length);


%% FFT test (For testing)
Fs = 1/theta_step;
ff = linspace(0,Fs,length(theta_long));
sin_test = cos(theta_long);
fft_sin_test=fft(sin_test);
[M max_index_sin]=max(abs(fft_sin_test));
M/length(theta_long)*2
fft_fundamental = zeros(size(sin_test));
fft_fundamental(1,max_index_sin)=fft_sin_test(1,max_index_sin);
fft_fundamental(1,end-max_index_sin+2)=fft_sin_test(1,max_index_sin);
ifft_sin = ifft(fft_fundamental);


%% Finding fundamental 
fft_N2b=fft(N2_bb);
[M max_index_fft_N2b]=max(abs(fft_N2b));
fft_fund_N2b= zeros(size(fft_N2b));
fft_fund_N2b(1,max_index_fft_N2b)=fft_N2b(1,max_index_fft_N2b);
fft_fund_N2b(1,end-max_index_fft_N2b+2)=fft_N2b(1,max_index_fft_N2b);
ifft_N2b = ifft(fft_fund_N2b);

figure(300)
plot(ff,abs(fft_N2b),'DisplayName','fft N2_b')
hold on;
plot(ff,fft_fund_N2b,'DisplayName','recovered fft N2_b')
grid on;
legend;
xlim([0 10])

figure(301)
plot(theta_long,N2_bb,'DisplayName','N2_b')
hold on;
plot(theta_long,ifft_N2b,'DisplayName','recovered N2_b')
grid on;
legend;

%% c) 
self_induc = 54/35

fft_N1_c=fft(N1_c);
[M max_index_fft_N1_c]=max(abs(fft_N1_c));
fft_fund_N1_c= zeros(size(fft_N1_c));
fft_fund_N1_c(1,max_index_fft_N1_c)=fft_N1_c(1,max_index_fft_N1_c);
fft_fund_N1_c(1,end-max_index_fft_N1_c+2)=fft_N1_c(1,max_index_fft_N1_c);
ifft_N1_c = ifft(fft_fund_N1_c);
peak_N1_c = M*2/length(N1_c)
30*4/pi

fft_N2_c=fft(N2_c);
[M max_index_fft_N2_c]=max(abs(fft_N2_c));
fft_fund_N2_c= zeros(size(fft_N2_c));
fft_fund_N2_c(1,max_index_fft_N2_c)=fft_N2_c(1,max_index_fft_N2_c);
fft_fund_N2_c(1,end-max_index_fft_N2_c+2)=fft_N2_c(1,max_index_fft_N2_c);
ifft_N2_c = ifft(fft_fund_N2_c);
peak_N2_c = M*2/length(N2_c)

ratio_c_fund = peak_N1_c^2/peak_N2_c^2

figure(400)
plot(ff,abs(fft_N1_c),'DisplayName','fft N1_c')
hold on;
plot(ff,abs(fft_fund_N1_c),'DisplayName','fundamental fft N1_c')
grid on;
legend;
xlim([0 10])

figure(401)
plot(theta_long,N1_c,'DisplayName','N1_c')
hold on;
plot(theta_long,imag(ifft_N1_c),'DisplayName','recovered N1_c')
grid on;
legend;

%% 
Ltot1 = pi/2
Ltot2 = pi/6
Ltot5 = 30^2*2*pi
Ltot6 = 3500*pi/3
Lh1 = Ltot1 - 4/pi
Lh2 = Ltot2-16/pi^3
Lh5 = Ltot5-14400/pi
Lh6 = Ltot6-1129*pi

Lper1 = Lh1/Ltot1
Lper2 = Lh2/Ltot2
Lper5 = Lh5/Ltot5
Lper6 = Lh6/Ltot6