function B = pa2B(p,a_bar)
% This function compute the B matrix corresponding to the Euler parameter p
% and a_bar in L-RF

e = p(2:4);
if size(e,2)>1
    e=e';
end
if size(a_bar,2)>1
    a_bar=a_bar';
end

B = 2*[(p(1)*eye(3)+cross_matrix(e))*a_bar, e*a_bar'-(p(1)*eye(3)+cross_matrix(e))*cross_matrix(a_bar)];
end