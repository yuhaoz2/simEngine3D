classdef simEngine3D < handle
% A collection of all the functions  for kinematics and dynamics analysis
% on 3D model. For now it only deal with 2 bodies system.

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
                      
                      % solve for velocity
                      this.compute_cons();
                      new_dq = this.Phi_q\this.nu;
                      
                      this.dq_i = new_dq(1:7);
                      this.dq_j = new_dq(8:14);
                      this.dr_i=this.dq_i(1:3);
                      this.dp_i=this.dq_i(4:7);
                      this.dr_j=this.dq_j(1:3);
                      this.dp_j=this.dq_j(4:7);
                      
                      % solve for acceleration
                      this.compute_cons();
                      new_ddq = this.Phi_q\this.gamma;
                      
                      this.ddq_i = new_ddq(1:7);
                      this.ddq_j = new_ddq(8:14);
                      this.ddr_i=this.ddq_i(1:3);
                      this.ddp_i=this.ddq_i(4:7);
                      this.ddr_j=this.ddq_j(1:3);
                      this.ddp_j=this.ddq_j(4:7);
                      
                  else
                      
                      % solve for position
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
                      
                      % solve for velocity
                      this.compute_cons();
                      new_dq = this.Phi_q\this.nu;
                      
                      this.dq_i = new_dq(1:7);
                      this.dr_i=this.dq_i(1:3);
                      this.dp_i=this.dq_i(4:7);
                      
                      % solve for acceleration
                      this.compute_cons();
                      new_ddq = this.Phi_q\this.gamma;
                      
                      this.ddq_i = new_ddq(1:7);
                      this.ddr_i=this.ddq_i(1:3);
                      this.ddp_i=this.ddq_i(4:7);
                  end
                  %disp(['t = ' num2str(time) ' sec.']);
              end
              states{k}.time = this.t;
              states{k}.q = [this.q_i;this.q_j];
              states{k}.dq = [this.dq_i;this.dq_j];
              states{k}.ddq = [this.ddq_i;this.ddq_j];
          end
      end
      

    end
end