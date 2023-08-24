function [cost] = funcion_costo_simple_MPC(chi,t,ts,N,chi_uav)

%% CONSTANTS VALUES OF THE ROBOT
a = 0.0; 
b = 0.0;
c = 0.0;
L = [a, b, c];

%% Definicion de los estados iniciales del sistema
h = [0;0;1;0];
v = [0; 0;0;0];

x = [h;v];

%% GENERAL VECTOR DEFINITION
H = [h;v];

%% Variables definidas por la TRAYECTORIA y VELOCIDADES deseadas
[hxd, hyd, hzd, hpsid, ul_d, um_d, un_d, r_d] = Trayectorias(3,t,5);

%% GENERALIZED DESIRED SIGNALS
%hd = [hxd; hyd; hzd; hpsid];
hd = [hxd;hyd;hzd;0*hpsid; ul_d; um_d; un_d; 0*r_d];
%hdp = [hxdp;hydp;hzdp;hpsidp];

%% Definicion de los limites de las acciondes de control
bounded = 3*[1.2; -1.2; 1.2; -1.2; 1.2; -1.2; 5.5; -5.5];

%% Definicion del vectro de control inicial del sistema
v_N = zeros(N,4);
H0 = repmat(H,1,N+1)'; 
x_N = H0;

% Definicion del optimizador
[f, solver, args] = mpc_drone_S_simple_gains(chi_uav,bounded, N, L, ts, chi);

tic
for k=1:length(t)-N

    %% Generacion del; vector de error del sistema
    he(:,k)=hd(1:4,k)-h(:,k);
     
    [u_opt,x_opt] = SolverUAV_MPC_din_gains(h,v,hd,N,x_N,v_N,args,solver,k);
      
    vref(:,k)= u_opt(1,:)';
    h_N(:,1:4,k) = x_opt(:,1:4);
 
    %% Dinamica del sistema   
    x(:,k+1) = UAV_Dinamica_RK4(chi_uav,x(:,k),vref(:,k),L,ts);
    
    h(:,k+1) = x(1:4,k+1);
    v(:,k+1) = x(5:8,k+1);
       
    %% Actualizacion de los resultados del optimizador para tener una soluciona aproximada a la optima
    v_N = [u_opt(2:end,:);u_opt(end,:)];
    x_N = [x_opt(2:end,:);x_opt(end,:)];
    
end
toc
cost = norm(he,2);

end

