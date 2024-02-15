%% Parameters for ECE 411 HW 2-2
clear all
clc

%% PMSM Machine Parameters
P   = 8;            % Number of poles
Rs   = 0.0124;       % Stator resistance per phase [Ohm]
Lq   = 150e-6;      % Stator q-axis inductance [H]
Ld   = 150e-6;      % Stator d-axis inductance [H]
lambda_f = 0.0185;     % PMSM Rotor flux   [N-m/A]
RPM_max = 4259;     % Maximum Speed     [RPM]
Vstep = 1.0;        % Step voltage amplitude [V]



