function [q_next,k1]=RK4_UAV_simple(q,con,ts)

    k1 = f_uav(q, con);   % new 
    k2 = f_uav(q + ts/2*k1, con); % new
    k3 = f_uav(q + ts/2*k2, con); % new
    k4 = f_uav(q + ts*k3, con); % new
    q_next= q +ts/6*(k1 +2*k2 +2*k3 +k4);
    q_next(4) = Angulo(q_next(4));
    
end