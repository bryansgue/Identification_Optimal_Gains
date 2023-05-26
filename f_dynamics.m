function [xp] = f_dynamics(h, u, chi)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
q = h(1:6);
q_p = h(7:12);

Mbar = M_matrix_bar(chi,q);
Cbar = C_matrix_bar(chi,q,q_p);
Gbar = G_matrix_bar(chi,q);

R = Rot_zyx(q(4:6));
% 
S = S_fuction(chi);
Q = Q_fuction(chi);
E = E_fuction(chi);
T = T_fuction(chi);

R_T = [R*T(1:3,1:3) T(1:3,4:6);T(4:6,1:3) T(4:6,4:6)];
Aux = (S*u-Q*q-E*q_p);
Aux1 = R*Aux(1:3,1);
Aux2 = Aux(4:6,1);
Input_model = [Aux1;Aux2];
q_pp = inv(Mbar+R_T)*(Input_model-Cbar*q_p-Gbar);


%  A = A_fuction(chi);
%  B = B_fuction(chi);
% 
% Aux1 = R*u(1:3,1);
% Aux2 = u(4:6,1);
% Input_model = A*q_p+B*[Aux1;Aux2];
% q_pp = inv(Mbar)*(Input_model-Cbar*q_p-Gbar);


xp = [q_p; q_pp];

end

