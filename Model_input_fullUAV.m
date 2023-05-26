function [u_real] = Model_input_fullUAV(X,u_ref, q, q_p, q_pp)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

 S = S_fuction(X);
 Q = Q_fuction(X);
 E = E_fuction(X);
 T = T_fuction(X);
% A = A_fuction(X);
% B = B_fuction(X);
u_real = S*u_ref-Q*q-E*q_p-T*q_pp;
%u_real = A*q_p+B*u_ref;
end