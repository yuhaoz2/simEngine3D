classdef simEngine3D_backup < handle
% A collection of all the functions  for kinematics and dynamics analysis
% on 3D model.

% Author: Yuhao Zhang
% Date: Oct, 2020
    
    properties
        nb; % number of bodies
        body; % collection of body info
        i; % Body 1 index, typically is 1
        j; % Body 2 index, typically is 0 or 2
        q_i; % generalized coordinates of body i, q = [r;p]
        dq_i; % derivative of generalized coordinates of body i, dq = [dr;dp]
        ddq_i; % second derivative of generalized coordinates of body i, ddq = [ddr;ddp]
        q_j; % generalized coordinates of body j, q = [r;p]
        dq_j; % derivative of generalized coordinates of body j, dq = [dr;dp]
        ddq_j; % second derivative of generalized coordinates of body j, ddq = [ddr;ddp]
        r_i;
        dr_i;
        ddr_i;
        p_i;
        dp_i;
        ddp_i;
        r_j;
        dr_j;
        ddr_j;
        p_j;
        dp_j;
        ddp_j;
        nc; % number of constraints (not including Euler parameter normalization constraint)
        constraints; % collection of kinematic and driving constraints
        t; % current time in system
        PhiF;
        PhiF_r;
        PhiF_p;
        PhiK;
        PhiK_r;
        PhiK_p;
        nu;
        gamma;
        lambda_p;
        lambda;
    end
    
    methods (Access = public)
      function this = simEngine3D(model_name) % build the model
          run(model_name); % this step run the driver script
          this.nb =nb;
          this.body = body;
          this.i=body(1).id;
          this.j=body(2).id;
          
          this.q_i=body(1).q_0;
          this.dq_i=body(1).dq_0;
          this.ddq_i=zeros(7,1);
          
          this.r_i=this.q_i(1:3);
          this.dr_i=this.dq_i(1:3);
          this.ddr_i=this.ddq_i(1:3);
          this.p_i=this.q_i(4:7);
          this.dp_i=this.dq_i(4:7);
          this.ddp_i=this.ddq_i(4:7);
          
          this.q_j=body(2).q_0;
          this.dq_j=body(2).dq_0;
          this.ddq_j=zeros(7,1);
          
          this.r_j=this.q_j(1:3);
          this.dr_j=this.dq_j(1:3);
          this.ddr_j=this.ddq_j(1:3);
          this.p_j=this.q_j(4:7);
          this.dp_j=this.dq_j(4:7);
          this.ddp_j=this.ddq_j(4:7);
          
          this.nc=nc;
          this.constraints = constraints;
          this.t=0;
          this.Phi = [];
          this.Phi_q=[];
          this.nu=[];
          this.gamma=[];
      end
      
      function compute_cons(this) % this function computer kinematic constraint values
          this.Phi = [];
          this.Phi_q=[];
          this.nu=[];
          this.gamma=[];
          for k = 1:this.nc
              switch  this.constraints(k).type
                  case 'DP1'
                      [Phi,nu,gamma,Phi_r,Phi_p]=GCon_DP1(this.i,this.q_i,this.dq_i,this.constraints(k).a_i_bar,this.j,this.q_j,this.dq_j,this.constraints(k).a_j_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),0);
                  case 'DP2'
                      [Phi,nu,gamma,Phi_r,Phi_p]=GCon_DP2(this.i,this.q_i,this.dq_i,this.constraints(k).a_i_bar,this.constraints(k).s_i_P_bar,this.j,this.q_j,this.dq_j,this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),0);
                  case 'D'
                      [Phi,nu,gamma,Phi_r,Phi_p]=GCon_D(this.i,this.q_i,this.dq_i,this.constraints(k).s_i_P_bar,this.j,this.q_j,this.dq_j,this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),0);
                  case 'CD'
                      [Phi,nu,gamma,Phi_r,Phi_p]=GCon_CD(this.constraints(k).c,this.i,this.q_i,this.dq_i,this.constraints(k).s_i_P_bar,this.j,this.q_j,this.dq_j,this.constraints(k).s_j_Q_bar,this.constraints(k).f(this.t),this.constraints(k).df(this.t),this.constraints(k).ddf(this.t),0);
                  otherwise
                      error('Constraint type incorrect.');
              end
              this.Phi = [this.Phi;Phi];
              this.Phi_q = [this.Phi_q;Phi_r, Phi_p];
              this.nu = [this.nu;nu];
              this.gamma = [this.gamma;gamma];
          end
          % add euler parameter normalization constraint
          Phi_euler = this.p_i'*this.p_i/2 -1/2;
          Phi_q_euler = [zeros(1,3) this.p_i'];
          nu_euler = 0;
          gamma_euler =  -this.dp_i'*this.dp_i;
          this.Phi = [this.Phi;Phi_euler];
          this.Phi_q = [this.Phi_q;Phi_q_euler];
          this.nu = [this.nu;nu_euler];
          this.gamma = [this.gamma;gamma_euler];
      end
      
      
      function states = kinematic_analysis(this,t_start,t_step,t_end)
          % perform kinematic analysis

          counts = length(t_start:t_step:t_end);
          states = cell(counts,1);% store states at each step
          
          for k = 1:counts
              time = t_start + (k-1)*t_step;
              this.t = time; % update system time
              
              if time > t_start
                  this.position_analysis();
                  this.velocity_analysis();
                  this.acceleration_analysis();
              end
              states{k}.time = this.t;
              states{k}.q = [this.q_i;this.q_j];
              states{k}.dq = [this.dq_i;this.dq_j];
              states{k}.ddq = [this.ddq_i;this.ddq_j];
          end
      end
      
      
      function position_analysis(this)
          % Kinematic Analysis Stage 1
          
          eps = 1e-6; % tolerance
          itCountMax =20; % max iteration
          
          if this.j~=0 % bodi j is not ground
              % solve for position
              new_q = [this.q_i;this.q_j];
              for h=1:itCountMax
                  this.compute_cons();
                  delta_q = this.Phi_q\this.Phi;
                  new_q =new_q - delta_q;
                  
                  this.q_i = new_q(1:7);
                  this.q_j = new_q(8:14);
                  this.r_i=this.q_i(1:3);
                  this.p_i=this.q_i(4:7);
                  this.r_j=this.q_j(1:3);
                  this.p_j=this.q_j(4:7);
                  
                  if norm(delta_q) < eps %check for convergence
                      break;
                  end
              end
              
              if h >= itCountMax
                  error('Maximum number of iterations reached');
              end
          else %  bodi j is ground
              new_q = this.q_i;
              for h=1:itCountMax
                  this.compute_cons();
                  delta_q = this.Phi_q\this.Phi;
                  new_q =new_q - delta_q;
                  
                  this.q_i = new_q(1:7);
                  this.r_i=this.q_i(1:3);
                  this.p_i=this.q_i(4:7);
                  
                  if norm(delta_q) < eps %check for convergence
                      break;
                  end
              end
              
              if h >= itCountMax
                  error('Maximum number of iterations reached');
              end
          end
      end
      
      function velocity_analysis(this)
          % Kinematic Analysis Stage 2
          
          if this.j~=0 % bodi j is not ground
              this.compute_cons();
              new_dq = this.Phi_q\this.nu;
              
              this.dq_i = new_dq(1:7);
              this.dq_j = new_dq(8:14);
              this.dr_i=this.dq_i(1:3);
              this.dp_i=this.dq_i(4:7);
              this.dr_j=this.dq_j(1:3);
              this.dp_j=this.dq_j(4:7);
          else %  bodi j is ground
              this.compute_cons();
              new_dq = this.Phi_q\this.nu;
              
              this.dq_i = new_dq(1:7);
              this.dr_i=this.dq_i(1:3);
              this.dp_i=this.dq_i(4:7);
              
          end
      end
      
      function acceleration_analysis(this)
          % Kinematic Analysis Stage 3
          
          if this.j~=0 % bodi j is not ground
              this.compute_cons();
              new_ddq = this.Phi_q\this.gamma;
              
              this.ddq_i = new_ddq(1:7);
              this.ddq_j = new_ddq(8:14);
              this.ddr_i=this.ddq_i(1:3);
              this.ddp_i=this.ddq_i(4:7);
              this.ddr_j=this.ddq_j(1:3);
              this.ddp_j=this.ddq_j(4:7);
          else %  bodi j is ground
              this.compute_cons();
              new_ddq = this.Phi_q\this.gamma;
              
              this.ddq_i = new_ddq(1:7);
              this.ddr_i=this.ddq_i(1:3);
              this.ddp_i=this.ddq_i(4:7);
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
              this.compute_cons();
              M = this.computeM();
              J_p = this.computeJ_p();
              P = this.computeP();
              F = this.computeF();
              tau_caret = this.compute_tau_caret();
              Phi_r = this.Phi_q(1:this.nc,1:3*(this.nb-1));
              Phi_p = this.Phi_q(1:this.nc,3*(this.nb-1)+1:end);
              
              % Build the matrix for inverse dynamics
              % Left
              L = [zeros(3*(this.nb-1),this.nb-1) Phi_r';...
                   P' Phi_p']; 
              %right
              R = -[M*this.ddr_i - F;...
                    J_p*this.ddp_i - tau_caret];  % need to update when there are more bodies
                
              lambda_comb = L\R;
              lambda_p = lambda_comb(1:this.nb-1);
              lambda = lambda_comb(this.nb:end);
                                        
              % STEP 3: recover the forces/torques
              
              torques{k} = cell(this.nc,this.nb-1);
              for h = 1:this.nb-1
                  Phi_p_h = Phi_p(:,(4*h-3):4*h);
                  G = p2G(this.p_i); % need to update when there are more bodies
                  for hh = 1:this.nc
                      torques{k}{hh,h} = -1/2*G*Phi_p_h(hh,:)'*lambda(hh);
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
          
          this.compute_cons();
          results{k}.time = this.t;
          results{k}.r = [this.r_i];
          results{k}.dr = [this.dr_i];
          results{k}.ddr = [this.ddr_i];
          results{k}.p = [this.p_i];
          results{k}.dp = [this.dp_i];
          results{k}.ddp = [this.ddp_i];
          results{k}.lambda_p = this.lambda_p;
          results{k}.lambda = this.lambda;
          results{k}.torques = cell(this.nc,this.nb-1);
          Phi_r = this.Phi_q(1:this.nc,1:3*(this.nb-1));
          Phi_p = this.Phi_q(1:this.nc,3*(this.nb-1)+1:end);
          for h = 1:this.nb-1
              Phi_p_h = Phi_p(:,(4*h-3):4*h);
              G = p2G(this.p_i); % need to update when there are more bodies
              for hh = 1:this.nc
                  results{k}.torques{hh,h} = -1/2*G*Phi_p_h(hh,:)'*this.lambda(hh);
              end
          end
          results{k}.violation_vel = Phi_r*this.dr_i + Phi_p*this.dp_i - this.nu(1:this.nc); % need to update when there are more bodies
          
          
          %%% Third, apply BDF to solve the dynamic problem %%%
          k = 2;
          time = t_start + (k-1)*t_step;
          this.t = time; % update system time to t_0
          this.BDF_solve_dynamic(1,t_step,results,k); % BDF method of order 1
          
          this.compute_cons();
          results{k}.time = this.t;
          results{k}.r = [this.r_i];
          results{k}.dr = [this.dr_i];
          results{k}.ddr = [this.ddr_i];
          results{k}.p = [this.p_i];
          results{k}.dp = [this.dp_i];
          results{k}.ddp = [this.ddp_i];
          results{k}.lambda_p = this.lambda_p;
          results{k}.lambda = this.lambda;
          results{k}.torques = cell(this.nc,this.nb-1);
          Phi_r = this.Phi_q(1:this.nc,1:3*(this.nb-1));
          Phi_p = this.Phi_q(1:this.nc,3*(this.nb-1)+1:end);
          for h = 1:this.nb-1
              Phi_p_h = Phi_p(:,(4*h-3):4*h);
              G = p2G(this.p_i); % need to update when there are more bodies
              for hh = 1:this.nc
                  results{k}.torques{hh,h} = -1/2*G*Phi_p_h(hh,:)'*this.lambda(hh);
              end
          end
          results{k}.violation_vel = Phi_r*this.dr_i + Phi_p*this.dp_i - this.nu(1:this.nc); % need to update when there are more bodies
          
          for k=3:counts
              time = t_start + (k-1)*t_step;
              this.t = time; % update system time to t_0
              this.BDF_solve_dynamic(BDF_order,t_step,results,k); % BDF method of order 2
              
              this.compute_cons();
              results{k}.time = this.t;
              results{k}.r = [this.r_i];
              results{k}.dr = [this.dr_i];
              results{k}.ddr = [this.ddr_i];
              results{k}.p = [this.p_i];
              results{k}.dp = [this.dp_i];
              results{k}.ddp = [this.ddp_i];
              results{k}.lambda_p = this.lambda_p;
              results{k}.lambda = this.lambda;
              results{k}.torques = cell(this.nc,this.nb-1);
              Phi_r = this.Phi_q(1:this.nc,1:3*(this.nb-1));
              Phi_p = this.Phi_q(1:this.nc,3*(this.nb-1)+1:end);
              for h = 1:this.nb-1
                  Phi_p_h = Phi_p(:,(4*h-3):4*h);
                  G = p2G(this.p_i); % need to update when there are more bodies
                  for hh = 1:this.nc
                      results{k}.torques{hh,h} = -1/2*G*Phi_p_h(hh,:)'*this.lambda(hh);
                  end
              end
              results{k}.violation_vel = Phi_r*this.dr_i + Phi_p*this.dp_i - this.nu(1:this.nc); % need to update when there are more bodies
          end
          
          
      end
      
      function check_pos_vel(this)
          % Check whether the initial positions and velocities satisfy
          % level one and level two constraint equations.
          
          eps = 1e-3;
          
          this.compute_cons();
          Phi_r = this.Phi_q(1:this.nc,1:3*(this.nb-1));
          Phi_p = this.Phi_q(1:this.nc,3*(this.nb-1)+1:end);
          P = this.computeP();
          
          if max(abs(this.Phi))>eps
              error("level 0 constraint is not satisfied");
          end
          
          for k=1:this.nb-1
              if abs(this.p_i'*this.p_i-1)>eps   % need to update when there are more bodies
                  error("level 0 constraint is not satisfied");
              end
          end
          
          if max(abs(Phi_r*this.dr_i+Phi_p*this.dp_i - this.nu(1:this.nc)))>eps % need to update when there are more bodies
              error("level 1 constraint is not satisfied");
          end
          
          if max(abs(P*this.dp_i))>eps 
              error("level 1 constraint is not satisfied");
          end
              
      end
      
      function compute_acceleration(this)
          % This function solve the initial acceleration using EOM and two
          % constraint equations
          
          this.compute_cons();
          M = this.computeM();
          J_p = this.computeJ_p();
          P = this.computeP();
          F = this.computeF();
          tau_caret = this.compute_tau_caret();
          Phi_r = this.Phi_q(1:this.nc,1:3*(this.nb-1));
          Phi_p = this.Phi_q(1:this.nc,3*(this.nb-1)+1:end);
          gamma_p = this.gamma(this.nc+1:end);
          gamma_caret = this.gamma(1:this.nc);
          % Left matrix
          L = [M  zeros(3*(this.nb-1),4*(this.nb-1))  zeros(3*(this.nb-1),(this.nb-1))  Phi_r';
               zeros(4*(this.nb-1),3*(this.nb-1))  J_p  P'  Phi_p';
               zeros((this.nb-1),3*(this.nb-1))  P  zeros(this.nb-1)  zeros((this.nb-1),this.nc);
               Phi_r  Phi_p  zeros(this.nc,(this.nb-1))  zeros(this.nc) ];
          
          % Right matrix
          R = [F; tau_caret; gamma_p; gamma_caret];
          
          unknown = L\R;
          
          % need to update when there are more bodies
          this.ddr_i = unknown(1:3*(this.nb-1));
          this.ddp_i = unknown(3*(this.nb-1)+1:7*(this.nb-1));
          this.ddq_i = [this.ddr_i; this.ddp_i];
          this.lambda_p = unknown(7*(this.nb-1)+1:8*(this.nb-1));
          this.lambda = unknown(8*(this.nb-1)+1:end);
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
                  s2 = recent_solution{n-2};
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
          
          this.ddr_i = s1.ddr;
          this.ddp_i = s1.ddp;
          this.ddq_i = [this.ddr_i;this.ddp_i];
          this.lambda_p = s1.lambda_p;
          this.lambda  = s1.lambda;
          
          for k = 1:itCountMax
              
              %%% Stage 1: compute position and velocity using BDF and most recent accelerations %%%
              
              this.r_i = C_r + beta^2*h^2*this.ddr_i;
              this.dr_i = C_r_dot + beta*h*this.ddr_i;
              this.p_i = C_p + beta^2*h^2*this.ddp_i;
              this.dp_i = C_p_dot + beta*h*this.ddp_i;
              this.q_i = [this.r_i;this.p_i];
              this.dq_i = [this.dr_i;this.dp_i];
              
              %%% Stage 2: compute the residual %%%
              
              this.compute_cons();
              M = this.computeM();
              J_p = this.computeJ_p();
              P = this.computeP();
              F = this.computeF();
              tau_caret = this.compute_tau_caret();
              Phi_r = this.Phi_q(1:this.nc,1:3*(this.nb-1));
              Phi_p = this.Phi_q(1:this.nc,3*(this.nb-1)+1:end);
              
              g = [M*this.ddr_i + Phi_r'*this.lambda - F;
                  J_p*this.ddp_i + Phi_p'*this.lambda + P'*this.lambda_p - tau_caret;
                  this.Phi(this.nc+1:end)/(beta^2*h^2);
                  this.Phi(1:this.nc)/(beta^2*h^2)];
              
              %%% Stage 3: solve linear system to get correction %%%
              
              % use Quasi Newton
              if k == 1
                  Psi = [M  zeros(3*(this.nb-1),4*(this.nb-1))  zeros(3*(this.nb-1),(this.nb-1))  Phi_r';
                        zeros(4*(this.nb-1),3*(this.nb-1))  J_p  P'  Phi_p';
                        zeros((this.nb-1),3*(this.nb-1))  P  zeros(this.nb-1)  zeros((this.nb-1),this.nc);
                        Phi_r  Phi_p  zeros(this.nc,(this.nb-1))  zeros(this.nc) ]; % approximate the Jacobian
              end
              
              delta = Psi\-g;
              
              %%% Stage 4: improve the quality of approximate solution %%%
              
              this.ddr_i   = this.ddr_i + delta(1:3*(this.nb-1));
              this.ddp_i   = this.ddp_i + delta(3*(this.nb-1)+1:7*(this.nb-1));
              this.ddq_i = [this.ddr_i;this.ddp_i];
              this.lambda_p = this.lambda_p + delta(7*(this.nb-1)+1:8*(this.nb-1));
              this.lambda  = this.lambda + delta(8*(this.nb-1)+1:end);
              
              %%% Stage 5: Go to Stage 6 if norm(delta) is small otherwise go back to stage 1 %%%
              
              if norm(delta) < eps
                  break;
              end
              
          end
          
          if k >= itCountMax
              error("Maximum number of iterations reached");
          end
          
          %%% Stage 6: get level 0 and level 1 variables %%%
          
          this.r_i     = C_r + beta^2*h^2*this.ddr_i;
          this.p_i     = C_p + beta^2*h^2*this.ddp_i;
          this.dr_i  = C_r_dot + beta*h*this.ddr_i;
          this.dp_i  = C_p_dot + beta*h*this.ddp_i;
          this.q_i = [this.r_i;this.p_i];
          this.dq_i = [this.dr_i;this.dp_i];
          
      end
      
      
      % The following are some supplement functions
      
      function M = computeM(this)
         % compute the mass matrix
         M = zeros(3*(this.nb-1),3*(this.nb-1)); % -1 because body 2 is ground
         for k = 1:this.nb-1
             M(3*k-2:3*k, 3*k-2:3*k) = this.body(k).mass*eye(3);
         end
      end
      
      function J_p = computeJ_p(this)
          % compute the inertia matrix
          J_p = zeros(4*(this.nb-1),4*(this.nb-1)); % -1 because body 2 is ground
          for k = 1:this.nb-1
             G = p2G(this.p_i); % need to update when there are more bodies
             J_p(4*k-3:4*k, 4*k-3:4*k) = 4*G'*this.body(k).inertia*G;
         end
      end
          
      function P = computeP(this)
          % compute the P matrix
          P = zeros(this.nb-1,4*(this.nb-1));
          for k = 1:this.nb-1
              P(k, 4*k-3:4*k) = this.p_i'; % need to update when there are more bodies
          end
      end
         
      function F = computeF(this)
          % compute the F vector
          F = zeros(3*(this.nb-1),1);
          for k = 1:this.nb-1
              F(3*k-2:3*k) = this.body(k).force; % total force acting on each body
          end
      end
      
      function tau_caret = compute_tau_caret(this)
          % compute the tau_caret vector
          tau_caret = zeros(4*(this.nb-1),1);
          for k = 1:this.nb-1
              G = p2G(this.p_i); % need to update when there are more bodies
              dG = p2G(this.dp_i); % need to update when there are more bodies
              tau_caret(4*k-3:4*k) = 2*G'*this.body(k).torque + 8*dG'*this.body(k).inertia*dG*this.p_i; % total torque acting on each body
          end
      end
     
      
    end
    
end