%% Parameters for ECE 411 Final Project - Solution
clear all
%% PMSM Machine Parameters
% P   =;            % Number of poles
% Rs   =;       % Stator resistance per phase [Ohm]
% Lq   = ;      % Stator q-axis inductance [H]
% Ld   = ;      % Stator d-axis inductance [H]
% lambda_f =;     % PMSM Rotor flux   [N-m/A]
Jm = 0.2;       % Rotational inertia [Nms]

RPM_max = 3700;     % Maximum Speed     [RPM]

lambda_f_est = 1.0*lambda_f;  % Set this to slightly higher or lower than 1.0 to model 
% the effects of mis-estimated parameters

% tau_kpI;       % Time const of proportional regulator
                        % I found the results are a little better with
                        % slightly faster proportional loop (smaller tau)
                        % than recommended in project statement.
% w_kpI = ;      % BW of the proprotional regulator
% w_kiI ;      % Freq where Kp and Ki/s cross over each other

Kp_ireg = 0.1416;  % qd current regulator proportional gain
Ki_ireg = 141.6    % qd current regulator intergator gain

Vdc = 72;             % Nominal Vdc level
Kp_wreg = Jm*30;
Ki_wreg = Kp_wreg*(2*pi*1);
Ls = Lq;
t_end = 10;    % Use something like 2 seconds for Part 4a, longer for Part 4b, 
% maybe 20 sec for extra credit section
Te_cmd_max_pos = 10;  
Te_cmd_max_neg = -20;  % Limit of torque command Nm

Te_test_part4a = -15;
w_test_Part4a_RPM = 3000;  % 1000 RPM also

v_wind_ramp_m_per_s_per_s = 0.5;    % 16 sec to increase up to 8 m/s


