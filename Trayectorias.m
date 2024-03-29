function [xd, yd, zd, psid, xd_p, yd_p, zd_p, psid_p] = Trayectorias(n,t,mul)


% 2) Seleccion de trayectoria deseada (POSICIÓN y VELOCIDADES)
     switch n
   % a)Trayectoria Churo     
        case 1   
            xd= 0.06 *t.* cos(0.2*t) ;   
            xd_p= 0.06 * cos(0.2*t) + 0.02 *0.2 *t.* (-sin(0.2*t));
            xd_pp= -0.06 * sin(0.2*t) + 0.02 *0.2 * (-sin(0.2*t)) - 0.02 * 0.2 * 0.2 * t.* cos(0.2*t);
            yd= 0.06 *t.* sin (0.2 * t);  
            yd_p= 0.06 * sin(0.2*t) + 0.02 *0.2*t.* cos(0.2*t);
            yd_pp= 0.06 * 0.2 * cos(0.2*t) + 0.02 *0.2* cos(0.2*t) - 0.02 * 0.2 * 0.2 *t.* sin(0.2*t);
            zd= 1 * sin (0.3 * t) + 10 ;
            zd_p= 1 * 0.3* cos(0.3*t);
            zd_pp= -1 * 0.3 * 0.3 * sin(0.3*t); 
   % b) Trayectoria seno   
        case 2 
            yd= 5 * sin(0.1*t) + 2;           yd_p= 5*0.1 * cos(0.1*t);         yd_pp= -5* 0.1* 0.1* sin(0.1*t);
            xd= 0.1*t;                        xd_p= 0.1* ones(1,length(t));     xd_pp= 0* ones(1,length(t));
            zd= 2 * ones(1,length(t)) +4;     zd_p= 0 * ones(1,length(t)); 
   % c) Trayectoria de un 8 
        case 3
            
            xd = 4 * sin(5*0.04*t) + 3;         xd_p = 4 * 5 * 0.04 * cos(5*0.04*t);     xd_pp = -4 * 5 * 0.04 * 5*0.04 * sin(5*0.04*t);
            yd = 4 * sin(5*0.08*t);              yd_p = 4 * 5 * 0.08 * cos(5*0.08*t);     yd_pp = -4 * 5 * 0.08 *5*0.08* sin(5*0.08*t);               
            zd = 2 * sin(5*0.08*t) + 5;             zd_p =2 * 5 * 0.08 * cos(5*0.08*t);
   % d) Trayectoria Silla de Montar
        case 4  
            xd= 5 * cos(0.05*t) + 5;                xd_p=-0.25*sin(0.05*t);           xd_pp=-0.0125*cos(0.05*t);
            yd= 5 * sin (0.05 * t) ;                yd_p=0.25*cos(0.05*t);            yd_pp=-0.0125*sin(0.05*t);
            zd= 1 * sin (0.3 * t) +18 ;             zd_p=0.3*cos(0.3*t);              %zd_pp=-0.09*sin(0.3*t);
   % e) Otra opcion
        otherwise
            disp("Ninuna de las anteriores");
     end

% 3) Cálculo de orientación 
     psid= (atan2(yd_p,xd_p));
     psid_p = (1./((yd_p./xd_p).^2+1)).*((yd_pp.*xd_p-yd_p.*xd_pp)./xd_p.^2);
     psid_p(1)=0;
end