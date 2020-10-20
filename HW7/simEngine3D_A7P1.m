% This is the main script to be excuted for the first mechanism model
% (Pendulum with revolute joint) of Assignment 7 Problem 1

clear

% build the model
sys = simEngine3D("revJoint"); % revJoint.m is the model info file

t_start = 0;
t_end = 10;
t_step = 1e-2;

torques = sys.inverse_dynamics(t_start,t_step,t_end);

t = t_start:t_step:t_end;

for i = 1: length(t)
    torque(i) = torques{i}{6,1}(3);
end

plot(t,torque)
title('Torque applied to pendulum')
xlabel('Time (sec)')
ylabel('Torque (Nm)')