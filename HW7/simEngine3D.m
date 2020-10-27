classdef simEngine3D < handle
% A collection of all the functions  for kinematics and dynamics analysis
% on 3D model. For now it only deal with 2-body system.

% Author: Yuhao Zhang
% Date: Oct, 2020
    
    properties
        nb; % number of bodies, 2 for now
        body; % collection of body info
        i; % Body 1 index, typically is 1
        j; % Body 2 index, typically is 0 or 2
        q_i; % generalized coordinates of body i, q = [r;p]
        dq_i; % derivative of generalized coordinates of body i, dq = [dr;dp]
        ddq_i; % second derivative of generalized coordinates of body i, dq = [dr;dp]
        q_j; % generalized coordinates of body j, q = [r;p]
        dq_j; % derivative of generalized coordinates of body j, dq = [dr;dp]
        ddq_j; % second derivative of generalized coordinates of body j, dq = [dr;dp]
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
        constraints; % collection of kinematic constraints
        t; % current time in system
        Phi;
        Phi_q;
        nu;
        gamma;
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