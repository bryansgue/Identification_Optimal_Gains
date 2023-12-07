function [cost] = funcion_costo_Ident_KinGain_Dinamica(X,hd_p,hd,N,ts,chi_real)

u(:,1) = [0;0;0;0];
uc_p(:,1) = [0;0;0;0];
h(:,1) = [0;0;1;0];
h_p(:,1) = [0;0;0;0];
he = [];

a=0;
b=0;
L=[a;b];

for k=1:N
    
    uc(:,k) = Control_UAV_s_for_ident(X,hd_p(:,k),hd(:,k),h(:,k));
    
    %% 2) ACELERATIONS VCP
   
    
    %% 2) DINAMICA DEL UAV (VELOCIDAD Y POSICION)
    u(:, k+1) = system_dynamic_s(chi_real, u(:,k), uc(:,k), h(:,k), L,ts);
    
    h(:,k+1) = RK4_UAV_simple(h(:,k),u(:,k),ts);
    
    he = [he; hd(:,k) - h(:,k)];
end
cost = norm(he,2);
%cost = norm(he'*he,2); %;
%cost = he'*he;
end

