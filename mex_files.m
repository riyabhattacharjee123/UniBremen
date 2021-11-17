close all
clear
% compile the s-functions and the c functions for the Simulink Model
%mex Simple_gravity_s_function.c Simple_grav.c;
%mex SRP_s_function.c SRP_acc.c;
%mex J2_effect_s_function.c J2_effect.c;
mex sun_grav_s_function.c sun_grav.c;
%mex air_drag_force_s_function.c air_drag.c;
mex moon_grav_s_function.c moon_grav.c;


% Run the hcw_cgo.m file
%run('D:\MSc. Project Zarm\Sample_s_function\solution_20_May\gco_hcw_triangle.m');
run('D:\MSc. Project Zarm\Sample_s_function\solution_20_May\gco_hcw_cartwheel.m');
%run('D:\MSc. Project Zarm\Sample_s_function\solution_20_May\gco_hcw_leader_follower.m');
%run('D:\MSc. Project Zarm\Sample_s_function\solution_20_May\gco_hcw_pendulum_two.m');
%run('D:\MSc. Project Zarm\Sample_s_function\solution_20_May\gco_hcw_pentagon.m');

% Run the generic init file
% Runs the 2017b version
%run('D:\MSc. Project Zarm\Sample_s_function\solution_20_May\init_generic.m');
% Runs the 2021b version
run('D:\MSc. Project Zarm\Sample_s_function\solution_20_May\init_general_formation.m');

%run the acceleration of each perturbation
run('D:\MSc. Project Zarm\Sample_s_function\solution_20_May\perturbed_acc_plots.m');
