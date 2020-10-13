% This is the main script to be excuted for the first mechanism model
% (Pendulum with revolute joint) of Assignment 6 Problem 2

clear

% read data
revJoint;

t = 0;
output = 0;

% compute kinematic constraints
[Phi_G_1,nu_G_1,gamma_G_1,Phi_r_G_1,Phi_p_G_1] = GCon_CD(c_1,i,q_i_0,dq_i_0,s_i_P_bar_1,j,q_j_0,dq_j_0,s_j_Q_bar_1,f_1,df_1,ddf_1,output);
[Phi_G_2,nu_G_2,gamma_G_2,Phi_r_G_2,Phi_p_G_2] = GCon_CD(c_2,i,q_i_0,dq_i_0,s_i_P_bar_2,j,q_j_0,dq_j_0,s_j_Q_bar_2,f_2,df_2,ddf_2,output);
[Phi_G_3,nu_G_3,gamma_G_3,Phi_r_G_3,Phi_p_G_3] = GCon_CD(c_3,i,q_i_0,dq_i_0,s_i_P_bar_3,j,q_j_0,dq_j_0,s_j_Q_bar_3,f_3,df_3,ddf_3,output);

[Phi_G_4,nu_G_4,gamma_G_4,Phi_r_G_4,Phi_p_G_4] = GCon_DP1(i,q_i_0,dq_i_0,a_i_bar_4,j,q_j_0,dq_j_0,a_j_bar_4,f_4,df_4,ddf_4,output);
[Phi_G_5,nu_G_5,gamma_G_5,Phi_r_G_5,Phi_p_G_5] = GCon_DP1(i,q_i_0,dq_i_0,a_i_bar_5,j,q_j_0,dq_j_0,a_j_bar_5,f_5,df_5,ddf_5,output);

% compute driving constraints
f_6 = double(f(t));
df_6 = double(df(t));
ddf_6 = double(ddf(t));
[Phi_D_6,nu_D_6,gamma_D_6,Phi_r_D_6,Phi_p_D_6] = GCon_DP1(i,q_i_0,dq_i_0,a_i_bar_6,j,q_j_0,dq_j_0,a_j_bar_6,f_6,df_6,ddf_6,output);

% compute euler parameter constraints
Phi_P_7 = p_i_0'*p_i_0/2 -1/2;
nu_P_7 = 0;
gamma_P_7 = -dp_i_0'*dp_i_0;
Phi_r_P_7 = zeros(1,3);
Phi_p_P_7 = p_i_0';

Phi = [Phi_G_1;Phi_G_2;Phi_G_3;Phi_G_4;Phi_G_5;Phi_D_6;Phi_P_7];
nu = [nu_G_1;nu_G_2;nu_G_3;nu_G_4;nu_G_5;nu_D_6;nu_P_7];
gamma = [gamma_G_1;gamma_G_2;gamma_G_3;gamma_G_4;gamma_G_5;gamma_D_6;gamma_P_7];
Phi_q = [Phi_r_G_1,Phi_p_G_1; Phi_r_G_2,Phi_p_G_2; Phi_r_G_3,Phi_p_G_3; Phi_r_G_4,Phi_p_G_4; Phi_r_G_5,Phi_p_G_5; Phi_r_D_6,Phi_p_D_6; Phi_r_P_7,Phi_p_P_7];

display(Phi);
display(Phi_q);
display(nu);
display(gamma);