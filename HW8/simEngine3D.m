classdef simEngine3D < handle
% A collection of all the functions  for kinematics and dynamics analysis
% on 3D model.

% Author: Yuhao Zhang
% Date: Oct, 2020
    
    properties
        nb; % number of bodies
        body; % collection of body info
        r; % position vector
        dr; % derivative of r
        ddr; % second order derivative of r
        p; % euler parameter vector
        dp; % derivative of p
        ddp; % second order derivative of p
        q; % generalized coordinates
        dq; % derivative of q
        ddq; % second order derivative of q
        nc; % number of kinematic constraints (not including Euler parameter normalization constraint)
        constraints; % collection of kinematic constraints
        t; % current time in system
        lambda_p; % lagrange multipliers from Euler parameter normalization constarints
        lambda; % lagrange multipliers from kinematic constarints
    end
    
    methods (Access = public)
        
      function this = simEngine3D(model_name) % build the model
          run(model_name); % this step run the driver script
          this.nb =nb;
          this.body = body;
          
          this.r = zeros(3*this.nb,1);
          this.dr = zeros(3*this.nb,1);
          this.ddr = zeros(3*this.nb,1);
          this.p = zeros(4*this.nb,1);
          this.dp = zeros(4*this.nb,1);
          this.ddp = zeros(4*this.nb,1);
          this.q = zeros(7*this.nb,1);
          this.dq = zeros(7*this.nb,1);
          this.ddq = zeros(7*this.nb,1);
          
          for i = 1:this.nb
              this.r((i-1)*3+1:3*i) = this.body(i).r_0;
              this.dr((i-1)*3+1:3*i) = this.body(i).dr_0;
              this.ddr((i-1)*3+1:3*i) = this.body(i).ddr_0;
              this.p((i-1)*4+1:4*i) = this.body(i).p_0;
              this.dp((i-1)*4+1:4*i) = this.body(i).dp_0;
              this.ddp((i-1)*4+1:4*i) = this.body(i).ddp_0;
              
              this.q((i-1)*7+1:7*i) = [this.body(i).r_0;this.body(i).p_0];
              this.dq((i-1)*7+1:7*i) = [this.body(i).dr_0;this.body(i).dp_0];
              this.ddq((i-1)*7+1:7*i) = [this.body(i).ddr_0;this.body(i).ddp_0];
          end
          
          this.nc=nc;
          this.constraints = constraints;
          this.t = 0;

      end
      
      
      
      
      
      function PhiK = get_PhiK(this)
          % compute kinematic constraints set
          PhiK = zeros(this.nc,1);
           for k = 1:this.nc
               i = this.constraints(k).body_i;
               j = this.constraints(k).body_j;
              switch  this.constraints(k).type
                  case 'DP1'
                      if j ~=0
                          [Phi_GCon,~,~,~,~]=GCon_DP1(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,...
                              j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).a_j_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),1);
                      else
                          [Phi_GCon,~,~,~,~]=GCon_DP1(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,...
                              j,[],[],this.constraints(k).a_j_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),1);
                      end
                  case 'DP2'
                      if j ~=0
                          [Phi_GCon,~,~,~,~]=GCon_DP2(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,this.constraints(k).s_i_P_bar,...
                              j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),1);
                      else
                          [Phi_GCon,~,~,~,~]=GCon_DP2(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,this.constraints(k).s_i_P_bar,...
                              j,[],[],this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),1);
                      end
                  case 'D'
                      if j ~=0
                          [Phi_GCon,~,~,~,~]=GCon_D(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                              j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),1);
                      else
                          [Phi_GCon,~,~,~,~]=GCon_D(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                              j,[],[],this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),1);
                      end
                  case 'CD'
                      if j ~=0
                          [Phi_GCon,~,~,~,~]=GCon_CD(this.constraints(k).c,i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                              j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),1);
                      else
                          [Phi_GCon,~,~,~,~]=GCon_CD(this.constraints(k).c,i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                              j,[],[],this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),1);
                      end
                  otherwise
                      error('Constraint type incorrect.');
              end
              PhiK(k) = Phi_GCon;
           end
      end
      
      function PhiP = get_PhiP(this)
          % compute Euler parameter normalization constraints
          PhiP = zeros(this.nb,1);
          for i = 1:this.nb
              PhiP(i) = this.p((i-1)*4+1:4*i)'*this.p((i-1)*4+1:4*i)/2 -1/2;
          end
      end
      
      function PhiF = get_PhiF(this)
          % compute overall constarints set
          PhiK = this.get_PhiK();
          PhiP = this.get_PhiP();
          PhiF = [PhiK; PhiP];
      end
      
      function nuK = get_nuK(this)
          % compute right hand side of velocity eqaution for kinematic constarints
           nuK = zeros(this.nc,1);
           for k = 1:this.nc
               i = this.constraints(k).body_i;
               j = this.constraints(k).body_j;
               if j ~=0
                   switch  this.constraints(k).type
                       case 'DP1'
                           [~,nu_GCon,~,~,~]=GCon_DP1(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,...
                               j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).a_j_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),2);
                       case 'DP2'
                           [~,nu_GCon,~,~,~]=GCon_DP2(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,this.constraints(k).s_i_P_bar,...
                               j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),2);
                       case 'D'
                           [~,nu_GCon,~,~,~]=GCon_D(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                               j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),2);
                       case 'CD'
                           [~,nu_GCon,~,~,~]=GCon_CD(this.constraints(k).c,i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                               j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),2);
                       otherwise
                           error('Constraint type incorrect.');
                   end
               else
                   switch  this.constraints(k).type
                       case 'DP1'
                           [~,nu_GCon,~,~,~]=GCon_DP1(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,...
                               j,[],[],this.constraints(k).a_j_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),2);
                       case 'DP2'
                           [~,nu_GCon,~,~,~]=GCon_DP2(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,this.constraints(k).s_i_P_bar,...
                               j,[],[],this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),2);
                       case 'D'
                           [~,nu_GCon,~,~,~]=GCon_D(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                               j,[],[],this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),2);
                       case 'CD'
                           [~,nu_GCon,~,~,~]=GCon_CD(this.constraints(k).c,i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                               j,[],[],this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),2);
                       otherwise
                           error('Constraint type incorrect.');
                   end
               end
               nuK(k) = nu_GCon;
           end
      end
      
      function nuP = get_nuP(this)
          % compute right hand side of velocity eqaution for Euler parameter normalization constarints
          nuP = zeros(this.nb,1);
      end
      
      function nuF = get_nuF(this)
          % compute right hand side of velocity eqaution for overall constarints
          nuK = this.get_nuK();
          nuP = this.get_nuP();
          nuF = [nuK; nuP];
      end
      
      function gammaK = get_gammaK(this)
          % compute right hand side of acceleration eqaution for kinematic constarints
           gammaK = zeros(this.nc,1);
           for k = 1:this.nc
               i = this.constraints(k).body_i;
               j = this.constraints(k).body_j;
               if j ~=0
                   switch  this.constraints(k).type
                       case 'DP1'
                           [~,~,gamma_GCon,~,~]=GCon_DP1(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,...
                               j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).a_j_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),3);
                       case 'DP2'
                           [~,~,gamma_GCon,~,~]=GCon_DP2(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,this.constraints(k).s_i_P_bar,...
                               j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),3);
                       case 'D'
                           [~,~,gamma_GCon,~,~]=GCon_D(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                               j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),3);
                       case 'CD'
                           [~,~,gamma_GCon,~,~]=GCon_CD(this.constraints(k).c,i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                               j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),3);
                       otherwise
                           error('Constraint type incorrect.');
                   end
               else
                   switch  this.constraints(k).type
                       case 'DP1'
                           [~,~,gamma_GCon,~,~]=GCon_DP1(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,...
                               j,[],[],this.constraints(k).a_j_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),3);
                       case 'DP2'
                           [~,~,gamma_GCon,~,~]=GCon_DP2(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,this.constraints(k).s_i_P_bar,...
                               j,[],[],this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),3);
                       case 'D'
                           [~,~,gamma_GCon,~,~]=GCon_D(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                               j,[],[],this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),3);
                       case 'CD'
                           [~,~,gamma_GCon,~,~]=GCon_CD(this.constraints(k).c,i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                               j,[],[],this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),3);
                       otherwise
                           error('Constraint type incorrect.');
                   end
               end
               gammaK(k) = gamma_GCon;
           end
      end
      
      function gammaP = get_gammaP(this)
          % compute right hand side of acceleration eqaution for Euler parameter normalization constarints
          gammaP = zeros(this.nb,1);
          for i = 1:this.nb
              gammaP(i) = -this.dp((i-1)*4+1:4*i)'*this.dp((i-1)*4+1:4*i);
          end
      end
      
      function gammaF = get_gammaF(this)
          % compute right hand side of acceleration eqaution for overall constarints
          gammaK = this.get_gammaK();
          gammaP = this.get_gammaP();
          gammaF = [gammaK; gammaP];
      end
      
      function PhiK_r = get_PhiK_r(this)
          % compute jacobian matrix of PhiK
          PhiK_r = zeros(this.nc,3*this.nb);
          
          for k = 1:this.nc
               i = this.constraints(k).body_i;
               j = this.constraints(k).body_j;
               if j ~=0
                   switch  this.constraints(k).type
                       case 'DP1'
                           [~,~,~,Phi_r_GCon,~]=GCon_DP1(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,...
                               j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).a_j_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),4);
                       case 'DP2'
                           [~,~,~,Phi_r_GCon,~]=GCon_DP2(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,this.constraints(k).s_i_P_bar,...
                               j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),4);
                       case 'D'
                           [~,~,~,Phi_r_GCon,~]=GCon_D(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                               j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),4);
                       case 'CD'
                           [~,~,~,Phi_r_GCon,~]=GCon_CD(this.constraints(k).c,i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                               j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),4);
                       otherwise
                           error('Constraint type incorrect.');
                   end
                   PhiK_r(k,3*(i-1)+1:3*i) = Phi_r_GCon(1:3);
                   PhiK_r(k,3*(j-1)+1:3*j) = Phi_r_GCon(4:6);
               else
                   switch  this.constraints(k).type
                       case 'DP1'
                           [~,~,~,Phi_r_GCon,~]=GCon_DP1(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,...
                               j,[],[],this.constraints(k).a_j_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),4);
                       case 'DP2'
                           [~,~,~,Phi_r_GCon,~]=GCon_DP2(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,this.constraints(k).s_i_P_bar,...
                               j,[],[],this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),4);
                       case 'D'
                           [~,~,~,Phi_r_GCon,~]=GCon_D(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                               j,[],[],this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),4);
                       case 'CD'
                           [~,~,~,Phi_r_GCon,~]=GCon_CD(this.constraints(k).c,i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                               j,[],[],this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),4);
                       otherwise
                           error('Constraint type incorrect.');
                   end
                   PhiK_r(k,3*(i-1)+1:3*i) = Phi_r_GCon(1:3);
               end

           end
      end
      
      function PhiK_p = get_PhiK_p(this)
          % compute jacobian matrix of PhiK
          PhiK_p = zeros(this.nc,4*this.nb);
          
          for k = 1:this.nc
               i = this.constraints(k).body_i;
               j = this.constraints(k).body_j;
               if j ~=0
                   switch  this.constraints(k).type
                       case 'DP1'
                           [~,~,~,~,Phi_p_GCon]=GCon_DP1(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,...
                               j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).a_j_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),4);
                       case 'DP2'
                           [~,~,~,~,Phi_p_GCon]=GCon_DP2(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,this.constraints(k).s_i_P_bar,...
                               j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),4);
                       case 'D'
                           [~,~,~,~,Phi_p_GCon]=GCon_D(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                               j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),4);
                       case 'CD'
                           [~,~,~,~,Phi_p_GCon]=GCon_CD(this.constraints(k).c,i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                               j,this.q((j-1)*7+1:7*j),this.dq((j-1)*7+1:7*j),this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),4);
                       otherwise
                           error('Constraint type incorrect.');
                   end
                   PhiK_p(k,4*(i-1)+1:4*i) = Phi_p_GCon(1:4);
                   PhiK_p(k,4*(j-1)+1:4*j) = Phi_p_GCon(5:8);
               else
                   switch  this.constraints(k).type
                       case 'DP1'
                           [~,~,~,~,Phi_p_GCon]=GCon_DP1(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,...
                               j,[],[],this.constraints(k).a_j_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),4);
                       case 'DP2'
                           [~,~,~,~,Phi_p_GCon]=GCon_DP2(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).a_i_bar,this.constraints(k).s_i_P_bar,...
                               j,[],[],this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),4);
                       case 'D'
                           [~,~,~,~,Phi_p_GCon]=GCon_D(i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                               j,[],[],this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),4);
                       case 'CD'
                           [~,~,~,~,Phi_p_GCon]=GCon_CD(this.constraints(k).c,i,this.q((i-1)*7+1:7*i),this.dq((i-1)*7+1:7*i),this.constraints(k).s_i_P_bar,...
                               j,[],[],this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),4);
                       otherwise
                           error('Constraint type incorrect.');
                   end
                   PhiK_p(k,4*(i-1)+1:4*i) = Phi_p_GCon(1:4);
               end

           end
      end
      
      function PhiF_q = get_PhiF_q(this)
          % jacobian matrix of PhiF
          PhiF_q = zeros(this.nc + this.nb,7*this.nb);
          PhiK_p = get_PhiK_p();
          PhiK_r = get_PhiK_r();
          
          for i = 1:this.nb
              PhiF_q(1:this.nc,7*(i-1)+1:7*(i-1)+3) = PhiK_r(:,3*(i-1)+1:3*i);
              PhiF_q(1:this.nc,7*(i-1)+4:7*i) = PhiK_p(:,4*(i-1)+1:4*i);
          end
              
          for i = 1:this.nb
              PhiF_q(this.nc+i,7*(i-1)+4:7*i) = this.p((i-1)*4+1:4*i)';
          end
      end
      
      
      
      
      
      function results = kinematic_analysis(this,t_start,t_step,t_end)
          % perform kinematic analysis

          counts = length(t_start:t_step:t_end);
          results = cell(counts,1);% store states at each step
          
          for k = 1:counts
              time = t_start + (k-1)*t_step;
              this.t = time; % update system time
              
              if time > t_start
                  this.position_analysis();
                  this.velocity_analysis();
                  this.acceleration_analysis();
              end
              results{k}.time = this.t;
              results{k}.q = this.q;
              results{k}.dq = this.dq;
              results{k}.ddq = this.ddq;
          end
      end
      
      
      function position_analysis(this)
          % Kinematic Analysis Stage 1
          
          eps = 1e-6; % tolerance
          itCountMax =20; % max iteration
          
          % solve for position
          new_q = this.q;
          for h=1:itCountMax
              PhiF_q = get_PhiF_q();
              PhiF = get_PhiF();
              delta_q = PhiF_q\PhiF;
              new_q =new_q - delta_q;
              
              this.q = new_q;
              for k = 1:this.nb
                  this.r(3*(k-1)+1:3*k)=this.q(7*(k-1)+1:7*(k-1)+3);
                  this.p(4*(k-1)+1:4*k)=this.q(7*(k-1)+4:7*k);
              end
              
              if norm(delta_q) < eps %check for convergence
                  break;
              end
          end
          
          if h >= itCountMax
              error('Maximum number of iterations reached');
          end
          
      end
      
      function velocity_analysis(this)
          % Kinematic Analysis Stage 2
          
          PhiF_q = get_PhiF_q();
          nuF = get_nuF();
          new_dq = PhiF_q\nuF;
          
          this.dq = new_dq;
          for k = 1:this.nb
              this.dr(3*(k-1)+1:3*k)=this.dq(7*(k-1)+1:7*(k-1)+3);
              this.dp(4*(k-1)+1:4*k)=this.dq(7*(k-1)+4:7*k);
          end
          
      end
      
      function acceleration_analysis(this)
          % Kinematic Analysis Stage 3
          
          PhiF_q = get_PhiF_q();
          gammaF = get_gammaF();
          new_ddq = PhiF_q\gammaF;
          
          this.ddq = new_ddq;
          for k = 1:this.nb
              this.ddr(3*(k-1)+1:3*k)=this.ddq(7*(k-1)+1:7*(k-1)+3);
              this.ddp(4*(k-1)+1:4*k)=this.ddq(7*(k-1)+4:7*k);
          end

      end
      
      function torques = inverse_dynamics(this,t_start,t_step,t_end)
          % compute the torques you need to get the prescribed motion via
          % inverse dynamics analysis
          
          counts = length(t_start:t_step:t_end);
          torques = cell(counts,1);% store torques at each step
          
          for k = 1:counts
              time = t_start + (k-1)*t_step;
              this.t = time; % update system time
              
              % STEP 1: solve for ddr and ddp
              
              %if time > t_start
                  this.position_analysis();
                  this.velocity_analysis();
                  this.acceleration_analysis();
              %end
              
              % STEP 2: solve for lagrange multiplier
              
              % Compute the needed matrices
              M = this.computeM();
              J_p = this.computeJ_p();
              P = this.computeP();
              F = this.computeF();
              tau_caret = this.compute_tau_caret();
              PhiK_r = this.get_PhiK_r();
              PhiK_p = this.get_PhiK_p();
              
              % Build the matrix for inverse dynamics
              % Left
              L = [zeros(3*this.nb,this.nb) PhiK_r';...
                   P' PhiK_p']; 
              % Right
              R = -[M*this.ddr - F;...
                    J_p*this.ddp - tau_caret]; 
                
              lambda_comb = L\R;
              this.lambda_p = lambda_comb(1:this.nb);
              this.lambda = lambda_comb(this.nb+1:end);
                                        
              % STEP 3: recover the forces/torques
              
              torques{k} = cell(this.nc,this.nb);
              for h = 1:this.nb
                  PhiK_p_h = PhiK_p(:,(4*h-3):4*h);
                  G = p2G(this.p(4*(h-1)+1:4*h));
                  for i = 1:this.nc
                      torques{k}{i,h} = -1/2*G*PhiK_p_h(i,:)'*this.lambda(i);
                  end
              end
          end
      end
      

      function results = dynamic_analysis(this,t_start,t_step,t_end,BDF_order)
          % This function perform a dynamic analysis using a Quasi-Newton approach
          
          counts = length(t_start:t_step:t_end);
          results = cell(counts,1);% store unknowns' values at each step
          
          %%% First, need to check initial position and velocity %%%
          k = 1;
          time = t_start + (k-1)*t_step;
          this.t = time; % update system time to t_0
          this.check_pos_vel();
          
          
          %%% Second, need to compute initial acceleration %%%
          this.compute_acceleration();
          
          results{k}.time = this.t;
          results{k}.r = this.r;
          results{k}.dr = this.dr;
          results{k}.ddr = this.ddr;
          results{k}.p = this.p;
          results{k}.dp = this.dp;
          results{k}.ddp = this.ddp;
          results{k}.lambda_p = this.lambda_p;
          results{k}.lambda = this.lambda;
          results{k}.torques = cell(this.nc,this.nb);
          
          PhiK_r = this.get_PhiK_r();
          PhiK_p = this.get_PhiK_p();
          nuK = this.get_nuK();
          for h = 1:this.nb
              PhiK_p_h = PhiK_p(:,(4*h-3):4*h);
              G = p2G(this.p(4*(h-1)+1:4*h));
              for i = 1:this.nc
                  results{k}.torques{i,h} = -1/2*G*PhiK_p_h(i,:)'*this.lambda(i);
              end
          end
          results{k}.violation_vel = PhiK_r*this.dr + PhiK_p*this.dp - nuK;
          
          
          %%% Third, apply BDF to solve the dynamic problem %%%
          k = 2;
          time = t_start + (k-1)*t_step;
          this.t = time; % update system time to t_0
          this.BDF_solve_dynamic(1,t_step,results,k); % BDF method of order 1
          
          results{k}.time = this.t;
          results{k}.r = this.r;
          results{k}.dr = this.dr;
          results{k}.ddr = this.ddr;
          results{k}.p = this.p;
          results{k}.dp = this.dp;
          results{k}.ddp = this.ddp;
          results{k}.lambda_p = this.lambda_p;
          results{k}.lambda = this.lambda;
          results{k}.torques = cell(this.nc,this.nb);
          
          PhiK_r = this.get_PhiK_r();
          PhiK_p = this.get_PhiK_p();
          nuK = this.get_nuK();
          for h = 1:this.nb
              PhiK_p_h = PhiK_p(:,(4*h-3):4*h);
              G = p2G(this.p(4*(h-1)+1:4*h));
              for i = 1:this.nc
                  results{k}.torques{i,h} = -1/2*G*PhiK_p_h(i,:)'*this.lambda(i);
              end
          end
          results{k}.violation_vel = PhiK_r*this.dr + PhiK_p*this.dp - nuK;
          
          for k=3:counts
              time = t_start + (k-1)*t_step;
              this.t = time; % update system time to t_0
              this.BDF_solve_dynamic(BDF_order,t_step,results,k); % BDF method of order 2
              
              results{k}.time = this.t;
              results{k}.r = this.r;
              results{k}.dr = this.dr;
              results{k}.ddr = this.ddr;
              results{k}.p = this.p;
              results{k}.dp = this.dp;
              results{k}.ddp = this.ddp;
              results{k}.lambda_p = this.lambda_p;
              results{k}.lambda = this.lambda;
              results{k}.torques = cell(this.nc,this.nb);
              
              PhiK_r = this.get_PhiK_r();
              PhiK_p = this.get_PhiK_p();
              nuK = this.get_nuK();
              for h = 1:this.nb
                  PhiK_p_h = PhiK_p(:,(4*h-3):4*h);
                  G = p2G(this.p(4*(h-1)+1:4*h));
                  for i = 1:this.nc
                      results{k}.torques{i,h} = -1/2*G*PhiK_p_h(i,:)'*this.lambda(i);
                  end
              end
              results{k}.violation_vel = PhiK_r*this.dr + PhiK_p*this.dp - nuK;
          end
          
          
      end
      
      function check_pos_vel(this)
          % Check whether the initial positions and velocities satisfy
          % level one and level two constraint equations.
          
          eps = 1e-3;
          
          PhiK = this.get_PhiK();
          if max(abs(PhiK))>eps
              error("level 0 constraint is not satisfied");
          end
          
          for k=1:this.nb
              if abs(this.p(4*(k-1)+1:4*k)'*this.p(4*(k-1)+1:4*k)-1)>eps 
                  error("level 0 constraint is not satisfied");
              end
          end
          
          PhiK_r = this.get_PhiK_r();
          PhiK_p = this.get_PhiK_p();
          nuK = this.get_nuK();
          if max(abs(PhiK_r*this.dr+PhiK_p*this.dp - nuK))>eps
              error("level 1 constraint is not satisfied");
          end
          
          P = this.computeP();
          if max(abs(P*this.dp))>eps 
              error("level 1 constraint is not satisfied");
          end
              
      end
      
      function compute_acceleration(this)
          % This function solve the initial acceleration using EOM and two
          % constraint equations

          M = this.computeM();
          J_p = this.computeJ_p();
          P = this.computeP();
          F = this.computeF();
          tau_caret = this.compute_tau_caret();
          PhiK_r = this.get_PhiK_r();
          PhiK_p = this.get_PhiK_p();
          gammaP = this.get_gammaP();
          gammaK = this.get_gammaK();
          % Left matrix
          L = [M  zeros(3*this.nb,4*this.nb)  zeros(3*this.nb,this.nb)  PhiK_r';
               zeros(4*this.nb,3*this.nb)  J_p  P'  PhiK_p';
               zeros(this.nb,3*this.nb)  P  zeros(this.nb)  zeros((this.nb),this.nc);
               PhiK_r  PhiK_p  zeros(this.nc,this.nb)  zeros(this.nc) ];
          
          % Right matrix
          R = [F; tau_caret; gammaP; gammaK];
          
          unknown = L\R;
          
          % need to update when there are more bodies
          this.ddr = unknown(1:3*this.nb);
          this.ddp = unknown(3*this.nb+1:7*this.nb);
          for k = 1:this.nb
              this.ddq(7*(k-1)+1:7*(k-1)+3) = this.ddr(3*(k-1)+1:3*k);
              this.ddq(7*(k-1)+4:7*k) = this.ddp(4*(k-1)+1:4*k);
          end
          this.lambda_p = unknown(7*this.nb+1:8*this.nb);
          this.lambda = unknown(8*this.nb+1:end);
      end
      
      function BDF_solve_dynamic(this,BDF_order,step_size,recent_solution,n)
          % This function use BDF method to solve dynamic problem
          
          eps = 1e-2; % tolerance
          itCountMax =20; % max iteration
          h = step_size;
                    
          switch BDF_order
              case 1
                  s1 = recent_solution{n-1}; % most recent solution
                  beta = 1;
                  a1 = -1;
                  C_r_dot = -a1*s1.dr;
                  C_r = -a1*s1.r + beta*h*C_r_dot;
                  C_p_dot = -a1*s1.dp;
                  C_p = -a1*s1.p + beta*h*C_p_dot;
                  
              case 2
                  s1 = recent_solution{n-1}; % most recent solution
                  s2 = recent_solution{n-2}; % second most recent solution
                  beta = 2/3;
                  a1 = -4/3;
                  a2 = 1/3;
                  C_r_dot = -a1*s1.dr -a2*s2.dr;
                  C_r = -a1*s1.r -a2*s2.r + beta*h*C_r_dot;
                  C_p_dot = -a1*s1.dp -a2*s2.dp;
                  C_p = -a1*s1.p -a2*s2.p + beta*h*C_p_dot;
                  
              otherwise
                  error("only support BDF method of order 1 or 2");
          end
          
          %%% Stage 0: prime new time step %%%
          
          this.ddr = s1.ddr;
          this.ddp = s1.ddp;
          for i = 1:this.nb
              this.ddq(7*(i-1)+1:7*(i-1)+3) = this.ddr(3*(i-1)+1:3*i);
              this.ddq(7*(i-1)+4:7*i) = this.ddp(4*(i-1)+1:4*i);
          end
          this.lambda_p = s1.lambda_p;
          this.lambda  = s1.lambda;
          
          for k = 1:itCountMax
              
              %%% Stage 1: compute position and velocity using BDF and most recent accelerations %%%
              
              this.r = C_r + beta^2*h^2*this.ddr;
              this.dr = C_r_dot + beta*h*this.ddr;
              this.p = C_p + beta^2*h^2*this.ddp;
              this.dp = C_p_dot + beta*h*this.ddp;
              for i = 1:this.nb
                  this.q(7*(i-1)+1:7*(i-1)+3) = this.r(3*(i-1)+1:3*i);
                  this.q(7*(i-1)+4:7*i) = this.p(4*(i-1)+1:4*i);
              end
              for i = 1:this.nb
                  this.dq(7*(i-1)+1:7*(i-1)+3) = this.dr(3*(i-1)+1:3*i);
                  this.dq(7*(i-1)+4:7*i) = this.dp(4*(i-1)+1:4*i);
              end
              
              %%% Stage 2: compute the residual %%%
              
              M = this.computeM();
              J_p = this.computeJ_p();
              P = this.computeP();
              F = this.computeF();
              tau_caret = this.compute_tau_caret();
              PhiK_r = this.get_PhiK_r();
              PhiK_p = this.get_PhiK_p();
              PhiP = this.get_PhiP();
              PhiK = this.get_PhiK();
              
              g = [M*this.ddr + PhiK_r'*this.lambda - F;
                  J_p*this.ddp + PhiK_p'*this.lambda + P'*this.lambda_p - tau_caret;
                  PhiP/(beta^2*h^2);
                  PhiK/(beta^2*h^2)];
              
              %%% Stage 3: solve linear system to get correction %%%
              
              % use Quasi Newton
              if k == 1
                  Psi = [M  zeros(3*this.nb,4*this.nb)  zeros(3*this.nb,this.nb)  PhiK_r';
                        zeros(4*this.nb,3*this.nb)  J_p  P'  PhiK_p';
                        zeros(this.nb,3*this.nb)  P  zeros(this.nb)  zeros(this.nb,this.nc);
                        PhiK_r  PhiK_p  zeros(this.nc,this.nb)  zeros(this.nc) ]; % approximate the Jacobian
              end
              
              delta = Psi\-g;
              
              %%% Stage 4: improve the quality of approximate solution %%%
              
              this.ddr   = this.ddr + delta(1:3*this.nb);
              this.ddp   = this.ddp + delta(3*this.nb+1:7*this.nb);
              for i = 1:this.nb
                  this.ddq(7*(i-1)+1:7*(i-1)+3) = this.ddr(3*(i-1)+1:3*i);
                  this.ddq(7*(i-1)+4:7*i) = this.ddp(4*(i-1)+1:4*i);
              end
              this.lambda_p = this.lambda_p + delta(7*this.nb+1:8*this.nb);
              this.lambda  = this.lambda + delta(8*this.nb+1:end);
              
              %%% Stage 5: Go to Stage 6 if norm(delta) is small otherwise go back to stage 1 %%%
              
              if norm(delta) < eps
                  break;
              end
              
          end
          
          if k >= itCountMax
              error("Maximum number of iterations reached");
          end
          
          %%% Stage 6: get level 0 and level 1 variables %%%
          
          this.r     = C_r + beta^2*h^2*this.ddr;
          this.p     = C_p + beta^2*h^2*this.ddp;
          this.dr  = C_r_dot + beta*h*this.ddr;
          this.dp  = C_p_dot + beta*h*this.ddp;
          for i = 1:this.nb
              this.q(7*(i-1)+1:7*(i-1)+3) = this.r(3*(i-1)+1:3*i);
              this.q(7*(i-1)+4:7*i) = this.p(4*(i-1)+1:4*i);
          end
          for i = 1:this.nb
              this.dq(7*(i-1)+1:7*(i-1)+3) = this.dr(3*(i-1)+1:3*i);
              this.dq(7*(i-1)+4:7*i) = this.dp(4*(i-1)+1:4*i);
          end
          
      end
      
      
      % The following are some supplement functions
      
      function M = computeM(this)
         % compute the mass matrix
         M = zeros(3*this.nb,3*this.nb);
         for k = 1:this.nb
             M(3*k-2:3*k, 3*k-2:3*k) = this.body(k).mass*eye(3);
         end
      end
      
      function J_p = computeJ_p(this)
          % compute the inertia matrix
          J_p = zeros(4*this.nb,4*this.nb);
          for k = 1:this.nb
             G = p2G(this.p(4*(k-1)+1:4*k));
             J_p(4*k-3:4*k, 4*k-3:4*k) = 4*G'*this.body(k).inertia*G;
         end
      end
          
      function P = computeP(this)
          % compute the P matrix
          P = zeros(this.nb,4*this.nb);
          for k = 1:this.nb
              P(k, 4*k-3:4*k) = this.p(4*(k-1)+1:4*k)';
          end
      end
         
      function F = computeF(this)
          % compute the F vector
          F = zeros(3*this.nb,1);
          for k = 1:this.nb
              F(3*k-2:3*k) = this.body(k).force; % total force acting on each body
          end
      end
      
      function tau_caret = compute_tau_caret(this)
          % compute the tau_caret vector
          tau_caret = zeros(4*this.nb,1);
          for k = 1:this.nb
              G = p2G(this.p(4*(k-1)+1:4*k));
              dG = p2G(this.dp(4*(k-1)+1:4*k));
              tau_caret(4*k-3:4*k) = 2*G'*this.body(k).torque + 8*dG'*this.body(k).inertia*dG*this.p(4*(k-1)+1:4*k); % total torque acting on each body
          end
      end
     
      
    end
    
end