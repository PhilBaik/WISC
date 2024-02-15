%%%%%%%%%%%%%%%%
%% Hyeongmeen Baik
%% Project #1
RPM_rads = 0.10472
J = 0.2
wNL_RPM = 4200 %% wNL = Va/K
Va = 72;

P_R = 5000
N_R = 3700
w_R = 3700*RPM_rads
I_R = 85

Tau_R = P_R/w_R
K = Tau_R/I_R

w_plot = linspace(0,3750*RPM_rads,10000);
Pm_plot = zeros(size(w_plot));
for k = 1:1:length(w_plot)
    if w_plot(1,k) < 3000*RPM_rads
        Pm_plot(1,k) = ((w_plot(1,k)/3000/RPM_rads)^3)*5e3;
    else
        Pm_plot(1,k) = 5e3;
    end
end

plot(w_plot,Pm_plot)

