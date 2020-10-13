% This is the main script to be excuted for the first mechanism model
% (Pendulum with revolute joint) of Assignment 6 Problem 2

clear

% build the model
sys = simEngine3D("revJoint"); % revJoint.m is the model info file
sys.t=0;
% compute the 4 quantities
sys.compute_cons();
% print out
disp("Phi(q,t) at t=0 is");
display(sys.Phi);
disp("Phi_q at t=0 is");
display(sys.Phi_q);
disp("nu at t=0 is");
display(sys.nu);
disp("gamma at t=0 is");
display(sys.gamma);