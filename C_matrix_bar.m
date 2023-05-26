function [Cbar] = C_matrix_bar(X,q,q_p)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
phi = q(4);
theta = q(5);
psi = q(6);

phi_p = q_p(4);
theta_p = q_p(5);
psi_p = q_p(6);

m1 = X(1);
Ixx = X(2);
Iyy = X(3);
Izz = X(4);

C11 = 0;
C12 = (Iyy - Izz)*(theta_p*cos(phi)*sin(phi)+psi_p*sin(phi)^2*cos(theta))+(Izz-Iyy)*(psi_p*cos(phi)^2*cos(theta)-Ixx*psi_p*cos(theta));
C13 = (Izz-Iyy)*psi_p*cos(phi)*sin(phi)*cos(theta)^2;
C21 = (Izz - Iyy)*(theta_p*cos(phi)*sin(phi)+phi_p*sin(phi)^2*cos(theta)) + (Iyy - Izz)*(psi_p*cos(phi)^2*cos(theta) + Ixx*psi_p*cos(theta));
C22 = (Izz-Iyy)*(phi_p*cos(phi)*sin(psi));
C23 = -Ixx*psi_p*sin(theta)*cos(theta)+Iyy*psi_p*sin(psi)^2*cos(theta)*sin(theta)+Izz*psi_p*cos(psi)^2*sin(theta)*cos(theta);
C31 = (Iyy-Izz)*(theta_p*cos(phi)*sin(phi)*cos(theta)^2)-Ixx*theta_p*cos(theta);
C32 = (Izz-Iyy)*(theta_p*cos(phi)*sin(phi)*sin(theta)+phi_p*sin(phi)^2*cos(theta))+(Iyy-Izz)*(phi_p*cos(phi)^2*cos(theta)) + Ixx*psi_p*cos(theta)*sin(theta)-Iyy*psi_p*sin(phi)^2*cos(theta)*sin(theta)-Izz*psi_p*cos(phi)^2*cos(theta)*sin(theta);
C33 = (Iyy-Izz)*(phi_p*cos(phi)*sin(phi)*cos(theta)^2) - Iyy*theta_p*sin(phi)^2*cos(theta)*sin(theta)-Izz*theta_p*cos(phi)^2*cos(theta)*sin(theta)+Ixx*theta_p*cos(theta)*sin(theta);

C = [C11 C12 C13;
    C21 C22 C23;
    C31 C32 C33];

Cbar = [zeros(3,3) zeros(3,3);
        zeros(3,3) C];
end

