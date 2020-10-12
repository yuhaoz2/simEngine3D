function a_tilde=cross_matrix(a)
% This function convert a vector "a" into its cross-product matrix form
% a_tilde.

a_tilde = [0 -a(3) a(2);...
           a(3) 0 -a(1);...
           -a(2) a(1) 0];
end