function [VMref] = Control_UAV_s_for_ident(X,hd_p,hd,h)
                                         
psi = h(4);

J11 = cos(psi);
J12 = -sin(psi);
J13 = 0;
J14 = 0;

J21 = sin(psi);
J22 = cos(psi);
J23 = 0;
J24 = 0;

J31 = 0;
J32 = 0;
J33 = 1;
J34 = 0;

J41 = 0;
J42 = 0;
J43 = 0;
J44 = 1;

J = [[J11 J12 J13 J14];[J21 J22 J23 J24];[J31 J32 J33 J34];[J41 J42 J43 J44]];
% J = [[J11 J12 J13 J14];[J21 J22 J23 J24];[J31 J32 J33 J34]];

%3) Calculos del Error

% xd = hd(1);
% yd = hd(2);
% zd = hd(3);
% psid = hd(4);
% 
% hx = h(1);
% hy = h(2);
% hz = h(3);
% psi = h(4);
% 
%   hxe= xd - hx;
%   hye= yd - hy;
%   hze= zd - hz;
%   psie= Angulo(psid-psi);   
%   he= [hxe ;hye ;hze ;psie];
  
 he = hd - h;
 he(4) = Angulo(he(4));
  
  %he = hd-h;
  %he(4) = Angulo(h(4));
% Constantes de ganancia ( ROS DJI_SDK)
K1 = zeros(4,4);
K1(1,1) = X(1);
K1(2,2) = X(2);
K1(3,3) = X(3);
K1(4,4) = X(4);

K2 = zeros(4,4);
K2(1,1) = X(5);
K2(2,2) = X(6);
K2(3,3) = X(7);
K2(4,4) = X(8);
   


% 7) Ley de control completa,  solucion = [u omega qpunto1 qpunto2]    
 VMref = pinv(J)*(hd_p +K2*tanh(inv(K2)*K1*he));
%    VMref = pinv(J)*(K1*tanh(K2*he'));
  
end