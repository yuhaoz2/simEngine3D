% This is a script to test two functions

%%data

r_i = [8, 6, -3];
p_i = [4, 3, -5, 1];
p_i = p_i / norm(p_i);
p_i_dot = [-0.2, 1.3, 3.4, 0];
p_i_dot(4) = -dot(p_i_dot, p_i)/p_i(4);
p_i_dot = p_i_dot / norm(p_i_dot);

a_i_bar = [-1.2, 1 ,0.3];
s_i_P_bar =[0.1, -0.3, 6.0];

r_j = [-0.5, 1.6, -6.3];
p_j = [3.3, -4, 5.1, 6];
p_j = p_j/norm(p_j);

p_j_dot = [0.6, -3.7, 5.1, 0];
p_j_dot(4) = -dot(p_j_dot, p_j) / p_j(4);
p_j_dot = p_j_dot / norm(p_j_dot);

a_j_bar = [1.2, 4.5, 3.1];
s_j_Q_bar = [0.2, -1.0, 1.5];

c = [0.3, 0.4, -6];

f = 1.2;
df = 2.5;
ddf = 0.2;
    
    
i = 1;
j = 2;
qi = [r_i p_i];
dqi = [0 0 0 p_i_dot];
qj = [r_j p_j];
dqj = [0 0 0 p_j_dot];


output = 0;


[Phi_DP2,nu_DP2,gamma_DP2,Phi_r_DP2,Phi_p_DP2] = GCon_DP2(i,qi,dqi,a_i_bar,s_i_P_bar,j,qj,dqj,s_j_Q_bar,f,df,ddf,output)

[Phi_D,nu_D,gamma_D,Phi_r_D,Phi_p_D] = GCon_D(i,qi,dqi,s_i_P_bar,j,qj,dqj,s_j_Q_bar,f,df,ddf,output)