clc
clear
close all
% initialize orbit calculation
r_start = [8000000 0 0 10000000 0 0 8500000 0 0];
v_start = [0 7000 0 0 5000 0 0 7500 0];

% start model
sim('model.slx')

 
% Plot orbit
plot(position(:,1),position(:,2))
hold on
plot(position(:,4),position(:,5))
hold on
plot(position(:,7),position(:,8))
axis equal
