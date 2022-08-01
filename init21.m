clear all
close all
% initialize orbit calculation
r_start = [8000000 0 0 10000000 0 0 12000000 0 0 14000000 0 0];
v_start = [0 7000 0 0 5000 0 0 2000 0 0 1800 0];

% start model
sim('model1.slx')

 
% Plot orbit
plot(position(:,1),position(:,2))
hold on
plot(position(:,4),position(:,5))
plot(position(:,7),position(:,8))
plot(position(:,10),position(:,11))
axis equal
