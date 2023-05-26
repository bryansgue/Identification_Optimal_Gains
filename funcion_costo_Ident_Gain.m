function [cost] = funcion_costo_Ident_Gain(X,hd_p,hd,N,ts)

h(:,1) = [0;0;1;0];
h_p(:,1) = [0;0;0;0];
he = [];
for k=1:N
    
    vc(:,k) = Control_UAV_s_for_ident(X,hd_p(:,k),hd(:,k),h(:,k));     
    [h(:,k+1),h_p(:,k+1)] = RK4_UAV_simple(h(:,k),vc(:,k),ts);

    he = [he; hd(:,k) - h(:,k)];
end
cost = norm(he,2);% + 0.05*norm(x,1);
%cost = norm(he'*he,2); %;
%cost = he'*he;
end

