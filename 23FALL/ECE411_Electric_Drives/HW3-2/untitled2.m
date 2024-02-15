%%%%%%%%%%%%%%%%%%%%%%
%%%%% ECE 411 HW 3-2
%%%%% Hyeongmeen Baik
%%%%%%%%%%%%%%%%%%%%%%

%% Q4 a
clc
clear
close all

Vs = 1;
we = 1;
r1 = 0.025;
r2 = 0.03;
x1 = 0.095;
x2 = 0.095;
xm = 3;
Zs = r1 + j*x1;
Zm = j*xm;

Zth = j*x2 + Zs*j*xm/(j*xm+Zs);
Zth_abs=abs(Zth)
r2_s = abs(Zth)

s_MT = r2/r2_s
Zr_MT = j*x2 + r2/s_MT;

wr = 1-s_MT

Is_ph = Vs/(Zs+Zr_MT*j*xm/(j*xm+Zr_MT))
Is_angle = angle(Is_ph)
Is_abs = abs(Is_ph)

Ir_ph = Is_ph*(j*xm/(j*xm+Zr_MT))
Ir_angle = angle(Ir_ph)
Ir_abs = abs(Ir_ph)

Te_MT = Ir_abs^2*r2/s_MT

s_MT_approx = r2/(x1+x2)
Te_MT_approx = Vs^2/(2*(x1+x2))


%% Q4 b

s = 0.1;
for i = 1:1:10
    iter(i) = fixed_freq(Vs,Zs,Zm,r2,x2,s);
    s = iter(i).s_new;
    
    if iter(i).error < 0.001
        break;
    end

end

final_error = iter(end).error;
s_final = iter(end-1).s_new;

wrR = 1-s_final;
IrR = iter(end).Ir;
IsR = iter(end).Is;
IsR_abs = abs(IsR);
IsR_angle = angle(IsR);
pfR = cos(angle(IsR));
PinR = abs(IsR)*pfR;
effR = iter(end).P/PinR*100;

%% Qb C
%% we = beta = 1
we = 1;
s = 0.02;
beta = 1
beta1 = variable_freq(1, xm,r1,x1,r2,x2,s,beta)
beta1.Te
beta1.Is
beta1.P
beta1.Te_abs = abs(beta1.Te)

we = 1;
s = 0.04;
beta = 0.5
beta_5 = variable_freq(1, xm,r1,x1,r2,x2,s,beta)
beta_5.Te
beta_5.Is
beta_5.P

we = 1;
s = -0.02;
beta = 1
beta1n = variable_freq(1, xm,r1,x1,r2,x2,s,beta)
beta1n.Te
beta1n.Is
beta1n.P
beta1n.Te_abs = abs(beta1n.Te)

we = 1;
s = -0.04;
beta = 0.5
beta_5n = variable_freq(1, xm,r1,x1,r2,x2,s,beta)
beta_5n.Te
beta_5n.Is
beta_5n.P


%% Functions
function out = fixed_freq(Vs, Zs, Zm,r2,x2, s)
    Zr = j*x2 +r2/s;
    Zs+Zm*Zr/(Zm+Zr)
    out.Is = Vs/(Zs+Zm*Zr/(Zm+Zr));
    out.Ir = out.Is*(Zm/(Zm+Zr));
    out.Te = abs(out.Ir)^2*r2/s;
    out.P = out.Te*(1-s);
    out.error = out.P-1;
    out.s_new = s/out.P;
end

function out = variable_freq(Vs0,xm, r1, x1 ,r2,x2, s, beta)
    Zs = r1 + j*beta*x1;
    Zm = j*beta*xm;
    Zr = j*x2*beta +r2/s;
    Vs = Vs0*beta;
    out.Zr = Zr
    out.Zs = Zs
    out.Zm = Zm
    out.Is = Vs/(Zs+Zm*Zr/(Zm+Zr));
    out.Ir = out.Is*(Zm/(Zm+Zr));
    out.Ir_abs = abs(out.Ir)
    out.Te = abs(out.Ir)^2*r2/s/beta;
    out.P = out.Te*(1-s)*beta;
    out.Te_abs = abs(out.Te);
    out.Is_abs = abs(out.Is);
    out.wr = (1-s)*beta

%     out.error = out.P-1;
%     out.s_new = s/out.P;
end