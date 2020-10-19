% Define number of bodies
nb = 2;

% body 1
body(1).id = 1;
r_0 = [0, sqrt(2), -sqrt(2)]';
p_0 = [0.653281482438188, 0.270598050073099, 0.653281482438188, 0.270598050073099]';
body(1).q_0 = [r_0; p_0];
dr_0 = [0, 0, 0]';
dp_0 = [0, 0, 0, 0]';
body(1).dq_0 = [dr_0; dp_0];

body(1).mass = 7800*4*0.05*0.05;
body(1).inertia = [body(1).mass/12*2*0.05^2 0 0;...
                   0 body(1).mass/12*(4^2+0.05^2) 0;...
                   0 0 body(1).mass/12*(4^2+0.05^2)];

% ground
body(2).id = 0;
body(2).q_0 = [0 0 0 1 0 0 0]';
body(2).dq_0 = [0 0 0 0 0 0 0]';

% Define number of constaints
nc = 6;
% 3 CD constraints
% constraint 1
constraints(1).type = "CD";
constraints(1).c = [1, 0, 0]';
constraints(1).s_i_P_bar = [-2, 0, 0]';
constraints(1).s_j_Q_bar =[0, 0, 0]';
constraints(1).f = @(t)(0*t);
constraints(1).df = @(t)(0*t);
constraints(1).ddf = @(t)(0*t);

% constraint 2
constraints(2).type = "CD";
constraints(2).c = [0, 1, 0]';
constraints(2).s_i_P_bar = [-2, 0, 0]';
constraints(2).s_j_Q_bar =[0, 0, 0]';
constraints(2).f = @(t)(0*t);
constraints(2).df = @(t)(0*t);
constraints(2).ddf = @(t)(0*t);

% constraint 3
constraints(3).type = "CD";
constraints(3).c = [0, 0, 1]';
constraints(3).s_i_P_bar = [-2, 0, 0]';
constraints(3).s_j_Q_bar =[0, 0, 0]';
constraints(3).f = @(t)(0*t);
constraints(3).df = @(t)(0*t);
constraints(3).ddf = @(t)(0*t);

% 2 DP1 constraints
% constraint 4
constraints(4).type = "DP1";
constraints(4).a_i_bar = [1, 0, 0]';
constraints(4).a_j_bar =[1, 0, 0]';
constraints(4).f = @(t)(0*t);
constraints(4).df = @(t)(0*t);
constraints(4).ddf = @(t)(0*t);

% constraint 5
constraints(5).type = "DP1";
constraints(5).a_i_bar = [0, 1, 0]';
constraints(5).a_j_bar =[1, 0, 0]';
constraints(5).f = @(t)(0*t);
constraints(5).df = @(t)(0*t);
constraints(5).ddf = @(t)(0*t);

% driving constraint
constraints(6).type = "DP1";
constraints(6).a_i_bar = [0, 1, 0]';
constraints(6).a_j_bar =[0, 0, -1]';
constraints(6).f = @(t)cos((pi/2)+(pi/4)*cos(2*t));
constraints(6).df = @(t)((pi*sin(2*t)*sin(pi/2 + (pi*cos(2*t))/4))/2);
constraints(6).ddf = @(t)(pi*cos(2*t)*sin(pi/2 + (pi*cos(2*t))/4) - (pi^2*sin(2*t)^2*cos(pi/2 + (pi*cos(2*t))/4))/4);