%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Author: Hyeongmeen Baik, Dheerej
%%%%%%%%%%%% Date: 02/13/2024
%%%%%%%%%%%% Title: ECE 711 - HW2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all
syms theta No Lo


%% a)
Ps = sym(2);
Pr = sym(2);

Nsa = (No/Ps)*(cos(Ps/sym(2)*theta));
Nsb = (No/Ps)*(cos(Ps/sym(2)*theta-sym(2*pi/3)));
Nsc = (No/Ps)*(cos(Ps/sym(2)*theta-sym(4*pi/3)));

Nra = (No/Pr)*sign(cos(Pr/sym(2)*theta));
Nrb = (No/Pr)*sign(cos(Pr/sym(2)*theta-sym(2*pi/3)));
Nrc = (No/Pr)*sign(cos(Pr/sym(2)*theta-sym(4*pi/3)));

% No = 1;

% self inductance
Las = Lo*int(Nsa^2,0,2*pi)
Lar = Lo*int(Nra^2,0,2*pi)

% fundamental self inductance
as = fund_extract(Nsa,Ps,theta);
Lms = Lo*(int(as^2,0,2*pi))

ar = fund_extract(Nra,Ps,theta);
Lmr = Lo*(int(ar^2,0,2*pi))

Lls = Las - Lms;
Llr = Lar - Lmr;

% getting estimated values
Lasev = sym(300);
Lmsev = vpa(Lms*Lasev/Las,6)
Llsev = vpa(Lls*Lasev/Las,6)
Larev = vpa(Lar*Lasev/Las,6)
Lmrev = vpa(Lmr*Lasev/Las,6)
Llrev = vpa(Llr*Lasev/Las,6)

N_ratio = sqrt(Lmr/Lms)

%%
No = 1;
figure(1)
subplot(211);
fplot(subs(Nsa),[0 2*pi],'LineWidth',2,'DisplayName','Nsa')
hold on;grid on;ylabel('Phase A');ylim([-No No]);
fplot(subs(Nra),[0 2*pi],'LineWidth',2,'DisplayName','Nra')
set(gca,'XTick',0:pi/3:2*pi);set(gca,'XTickLabel',{'0','\pi/3','2\pi/3','\pi','4\pi/3','5\pi/3','2\pi'})
legend

%% b)
syms No theta_r t
assume(t,"real")
assume(theta_r,"real")

% for computational issue, assigning theta_r to stator winding function is
% much faster
Nsa = (No/Ps)*(cos(Ps/sym(2)*theta+theta_r));
Nsb = (No/Ps)*(cos(Ps/sym(2)*theta+theta_r-sym(2*pi/3)));
Nsc = (No/Ps)*(cos(Ps/sym(2)*theta+theta_r-sym(4*pi/3)));

Nra = (No/Pr)*sign(cos(Pr/sym(2)*theta));
Nrb = (No/Pr)*sign(cos(Pr/sym(2)*theta-sym(2*pi/3)));
Nrc = (No/Pr)*sign(cos(Pr/sym(2)*theta-sym(4*pi/3)));

% Mutual indutances
Lasbs = Lo*int(Nsa*Nsb,0,2*pi);
Larbr = Lo*int(Nra*Nrb,0,2*pi);
Lasar = Lo*int(Nsa*Nra,0,2*pi);
Lasbr = Lo*int(Nsa*Nrb,0,2*pi);
Lascr = Lo*int(Nsa*Nrb,0,2*pi);
Lar = Lo*int(Nra*Nra,0,2*pi);

% getting estimated values
Lasev = sym(300);
Larev = Lar*Lasev/Las;
Lasbsev = Lasbs*Lasev/Las;
Larbrev = Larbr*Lasev/Las;
Lasarev = Lasar*Lasev/Las;
Lasbrev = Lasbr*Lasev/Las;
Lascrev = Lascr*Lasev/Las;

% set up constraints
w = sym(377);
Im = sym(4);

% Ldi/dt = v
ias = Im*sin(w*t)
vas = vpa(Lasev/1000*diff(ias,t),6)
vbs = vpa(Lasbsev/1000*diff(ias,t),6)
vcs = vpa(Lasbsev/1000*diff(ias,t),6)
var = vpa(Lasarev/1000*diff(ias,t),6)
vbr = vpa(Lasbrev/1000*diff(ias,t),6)
vcr = vpa(Lasbrev/1000*diff(ias,t),6)

%% c)
Labcr_s = [Lasarev, Lasbrev, Lasbrev;Lasbrev, Lasarev, Lasbrev;Lasbrev, Lasbrev, Lasarev]/1000
Labcr_r = [Larev, Larbrev, Larbrev;Larbrev, Larev, Larbrev;Larbrev, Larbrev, Larev]/1000
Labcs_s = [Lasev, Lasbsev, Lasbsev;Lasbsev, Lasev, Lasbsev;Lasbsev, Lasbsev, Lasev]/1000
Labcs_r = [Lasarev, Lasbrev, Lasbrev;Lasbrev, Lasarev, Lasbrev;Lasbrev, Lasbrev, Lasarev]/1000

% lambda abcr should be zero
% iabcr = - inv(Labcr_r)*Labcr_s*iabcs

iabcs = [ias;0;0]
iabcr = -inv(Labcr_r)*Labcr_s*iabcs
iar = vpa(iabcr(1,1),6)
ibr = vpa(iabcr(2,1),6)
icr = vpa(iabcr(3,1),6)

Lambda_abcs_r = Labcs_r*iabcr
Lambda_abcs_s = Labcs_s*iabcs
Lambda_abcs = Lambda_abcs_s + Lambda_abcs_r

Lambdaabcs = diff(Lambda_abcs,t)
vas = vpa(Lambdaabcs(1,1),6)
vbs = vpa(Lambdaabcs(2,1),6)
vcs = vpa(Lambdaabcs(3,1),6)
% vas = Las*dias/dt + Nr/Ns*Lms*[d(iar*cos(theta_r))/dt 

%% d) and e)
ias = Im*sin(w*t)
ibs = Im*sin(w*t-2*sym(pi)/3)
ics = Im*sin(w*t-4*sym(pi)/3)
a = exp(j*2*pi/3)
Iabcs = simplify(2/3*(ias+a*ibs+a^2*ics))
Lambdaabcs = simplify((Lasev - Lasbsev)*Iabcs)/1000
Lambdabcr = expand((Lasarev - Lasbrev)*Iabcs/1000,'ArithmeticOnly',true)
vabcs = diff(Lambdaabcs,t)
vabcr = vpa(diff(Lambdabcr,t),4)
vas = vpa(real(vabcs),4)
vas_check = vpa(real(vabcs),4)
vbs = vpa(real(vabcs*a^2),4)
vcs = vpa(real(vabcs*a),4)

var = vpa(real(vabcr),4)
vbr = vpa(real(vabcr*a^2),4)
vcr = vpa(real(vabcr*a),4)


%% f)

iabcs = [ias;ibs;ics]
iabcr = -inv(Labcr_r)*Labcr_s*iabcs

iar = iabcr(1,1);
ibr = iabcr(2,1);
icr = iabcr(3,1);

i0r = simplify(iar+ibr+icr)

Lambda_abcs_r = Labcs_r*iabcr
Lambda_abcs_s = Labcs_s*iabcs
Lambda_abcs = Lambda_abcs_s + Lambda_abcs_r

Lambdaabcs = Lambda_abcs
vabcs = diff(Lambdaabcs,t)
vas = simplify(vabcs(1,1))
vbs = simplify(vabcs(2,1))
vcs = simplify(vabcs(3,1))

vabcs = simplify(2/3*(vas+vbs*a+vcs*a));
iabcr = simplify(2/3*(iar+ibr*a+icr*a));

vas_recovered = real(vabcs)
err = simplify(vas-vas_recovered)

%% Appendix
% verification between matrix method and complex vector method (e)
ias = Im*sin(w*t)
ibs = Im*sin(w*t-2*sym(pi)/3)
ics = Im*sin(w*t-4*sym(pi)/3)
iabcs = [ias;ibs;ics]
Lambda_abcs_s = Labcs_s*iabcs
Lambda_abcs = Lambda_abcs_s
Lambdaabcs = diff(Lambda_abcs,t)
vas_temp = vpa(Lambdaabcs(1,1),6)
vbs = vpa(Lambdaabcs(2,1),6)
vcs = vpa(Lambdaabcs(3,1),6)

theta_r = 0
No = 1;
figure(3)
subplot(211);
fplot(subs(vas_check),[0 pi/10],'LineWidth',2,'DisplayName','iar')
hold on;grid on;ylabel('Phase A');
fplot(subs(vas_temp),[0 pi/10],'--','LineWidth',2,'DisplayName','ias')
% set(gca,'XTick',0:pi/3:2*pi);set(gca,'XTickLabel',{'0','\pi/3','2\pi/3','\pi','4\pi/3','5\pi/3','2\pi'})
legend

% recovered signal from f)
subplot(212);
fplot(subs(vas),[0 pi/10],'LineWidth',2,'DisplayName','vas')
hold on;grid on;ylabel('Phase A');
fplot(subs(vas_recovered),[0 pi/10],'--','LineWidth',2,'DisplayName','vas _ recovered')
% set(gca,'XTick',0:pi/3:2*pi);set(gca,'XTickLabel',{'0','\pi/3','2\pi/3','\pi','4\pi/3','5\pi/3','2\pi'})
legend
%% functions
% fundamental
function out = fund_extract(Nsa,Ps,theta)
    % fundamental self inductance
    C_pos1 = 1/2/sym(pi)*int(Nsa*exp(-j*theta*2/Ps),0,2*pi);
    C_neg1 = 1/2/sym(pi)*int(Nsa*exp(j*theta*2/Ps),0,2*pi);
    k = sqrt((C_pos1+C_neg1)^2+(j*C_pos1-j*C_neg1)^2);
    if C_neg1+C_pos1 ==0
        if subs(j*(C_pos1-C_neg1),No,1)>0
            phi = pi/2;
        else 
            phi = -pi/2;
        end
    else
        phi = atan(j*((C_pos1-C_neg1)/(C_neg1+C_pos1))); %%% phi = -pi/2
    end
    out = k*sin(theta+phi);
end