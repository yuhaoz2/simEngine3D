% This is the main script to be excuted for the double Pendulum model
% of Assignment 8 Problem 2

clear

% build the model
sys = simEngine3D("double_pen"); % double_pen.m is the model info file

% Assign start time, ending time and time step
t_start = 0;
t_end = 10;
t_step = 1e-3;
BDF_order = 1;

% Perform the dynamic analysis
start = tic;
results = sys.dynamic_analysis(t_start,t_step,t_end,BDF_order);
run_time = toc(start);

disp(['The simulation time is ' num2str(run_time) ' seconds, and step size is ' num2str(t_step)]);

t = t_start:t_step:t_end;


%% Plots

% Plot coordiantes of Body 1's point O'

r = zeros(6,length(results));
for k = 1:length(results)
    r(:,k) = results{k}.r; % position of Body 1's and Body 2's O' 
end

figure 
subplot(3,1,1)
plot(t,r(1,:));
title('x coordinate function plot of Body 1`s point O`');
xlabel('Time (sec)');
ylabel('x coordinate (m)');

subplot(3,1,2)
plot(t,r(2,:));
title('y coordinate function plot of Body 1`s point O`');
xlabel('Time (sec)');
ylabel('y coordinate (m)');

subplot(3,1,3)
plot(t,r(3,:));
title('z coordinate function plot of Body 1`s point O`');
xlabel('Time (sec)');
ylabel('z coordinate (m)');

% Plot coordiantes of Body 2's point O'

figure 
subplot(3,1,1)
plot(t,r(4,:));
title('x coordinate function plot of Body 2`s point O`');
xlabel('Time (sec)');
ylabel('x coordinate (m)');

subplot(3,1,2)
plot(t,r(5,:));
title('y coordinate function plot of Body 2`s point O`');
xlabel('Time (sec)');
ylabel('y coordinate (m)');

subplot(3,1,3)
plot(t,r(6,:));
title('z coordinate function plot of Body 2`s point O`');
xlabel('Time (sec)');
ylabel('z coordinate (m)');

% Plot the angular velocity of Body 1 and Body 2

p = zeros(8,length(results));
dp = zeros(8,length(results));
omega1 = zeros(3,length(results));
omega2 = zeros(3,length(results));

for k = 1:length(results)
    % compute orientation matrix
    p     = results{k}.p; 
    dp    = results{k}.dp; 

    omega1(:,k) = 2*p2E(p(1:4))*dp(1:4);
    omega2(:,k) = 2*p2E(p(5:8))*dp(5:8);
end

figure 
subplot(3,1,1)
plot(t,omega1(1,:));
title('Plot of Body 1`s angular velocity in G-RF x axis');
xlabel('Time (sec)');
ylabel('omega_x (rad/s)');

subplot(3,1,2)
plot(t,omega1(2,:));
title('Plot of Body 1`s angular velocity in G-RF y axis');
xlabel('Time (sec)');
ylabel('omega_y (rad/s)');

subplot(3,1,3)
plot(t,omega1(3,:));
title('Plot of Body 1`s angular velocity in G-RF z axis');
xlabel('Time (sec)');
ylabel('omega_z (rad/s)');

figure 
subplot(3,1,1)
plot(t,omega2(1,:));
title('Plot of Body 2`s angular velocity in G-RF x axis');
xlabel('Time (sec)');
ylabel('omega_x (rad/s)');

subplot(3,1,2)
plot(t,omega2(2,:));
title('Plot of Body 2`s angular velocity in G-RF y axis');
xlabel('Time (sec)');
ylabel('omega_y (rad/s)');

subplot(3,1,3)
plot(t,omega2(3,:));
title('Plot of Body 2`s angular velocity in G-RF z axis');
xlabel('Time (sec)');
ylabel('omega_z (rad/s)');


% Plot the 2-norm of the violation of the velocity constraint equations for the revolute joint between Body 1 and Body 2

vio = zeros(1,length(results));
for k = 1:length(results)
    vio(k) = norm(results{k}.violation_vel(6:10));
end

figure
plot(t,vio);
title('2-norm of the violation of the velocity constraint equations');
xlabel('Time (sec)');
ylabel('2-norm of violation');