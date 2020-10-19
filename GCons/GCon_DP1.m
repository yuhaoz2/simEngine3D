function [Phi,nu,gamma,Phi_r,Phi_p] = GCon_DP1(i,qi,dqi,ai_bar,j,qj,dqj,aj_bar,ft,dft,ddft,output)
% Input: Body number: i, j;
%        Generalized coordinates: qi,qj (qi=[ri;pi]), and their
%        derivatives: dqi, dqj;
%        Two algebric vectors in L-RF: ai_bar, aj_bar;
%        The prescribed value at current time t: ft, and its first order
%        derivative: dft, and second order derivative: ddft;
%        Output flag: output (1 for Phi, 2 for nu, 3 for gamma, 4 for
%        Phi_r/Phi_p, 0 for all).

% Output: Phi: the value of the constraint;
%         nu: the right-hand side of the velocity equation;
%         gamma: The right-hand side of the acceleration equation;
%         Phi_r/Phi_p: The expression of the partial derivatives.

%% initialize

Phi = 0;
nu = 0;
gamma = 0;
Phi_r = [];
Phi_p = [];

% Return all quantities if not specified
if nargin < 12
    output = 0;
elseif ~(output == 0 || 1 || 2 || 3 || 4)
    warning('output should be 1 for Phi, 2 for nu, 3 for gamma, 4 for Phi_r/Phi_p, 0 for all');
end

if j==0
    qj = [0 0 0 1 0 0 0]';
    dqj = zeros(7,1);
end
% separate ri and pi from qi
if size(qi,2)>1
    qi = qi';
end
if size(qj,2)>1
    qj = qj';
end
if size(dqi,2)>1
    dqi = dqi';
end
if size(dqj,2)>1
    dqj = dqj';
end
%ri=qi(1:3);
%rj=qj(1:3);
pi=qi(4:7);
pj=qj(4:7);
%dri=dqi(1:3);
%drj=dqj(1:3);
dpi=dqi(4:7);
dpj=dqj(4:7);

% make sure a_bar are column vectors
if size(ai_bar,2)>1
    ai_bar = ai_bar';
end
if size(aj_bar,2)>1
    aj_bar = aj_bar';
end

% compute the orientation matrices
if j ==0
    Aj = eye(3);
else
    Aj = p2A(pj);
end
Ai = p2A(pi);

%% compute the quantities

% compute Phi
if output == 0 || 1
   Phi = ai_bar'*Ai'*Aj*aj_bar - ft; 
end

% compute nu
if output == 0 || 2
   nu = dft; 
end
    
% compute gamma
if output == 0 || 3
    if j~=0
        gamma = -ai_bar'*Ai'*pa2B(dpj,aj_bar)*dpj - aj_bar'*Aj'*pa2B(dpi,ai_bar)*dpi - 2*(pa2B(pi,ai_bar)*dpi)'*(pa2B(pj,aj_bar)*dpj) + ddft;
    else
        gamma = - aj_bar'*Aj'*pa2B(dpi,ai_bar)*dpi + ddft;
    end
end

% compute Phi_r/Phi_p
if output == 0 || 4
    if j ~=0
        Phi_r = zeros(1,6);
        Phi_p = [aj_bar'*Aj'*pa2B(pi,ai_bar), ai_bar'*Ai'*pa2B(pj,aj_bar)];
    else
        Phi_r = zeros(1,3);
        Phi_p = [aj_bar'*Aj'*pa2B(pi,ai_bar)];
    end
end

end