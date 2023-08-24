function [cost] = funcion_costo_Ident_Gain_MPC(gains,t,ts,N,chi_real)

%% CONSTANTS VALUES OF THE ROBOT
a = 0.0; 
b = 0.0;
c = 0.0;
L = [a, b, c];

% Definicion de los estados iniciales del sistema

h = [0;0;1;0];

%% INITIAL GENERALIZE VELOCITIES
v = [0; 0;0;0];
x = [h;v];

%% GENERAL VECTOR DEFINITION
H = [h;v];

%% Variables definidas por la TRAYECTORIA y VELOCIDADES deseadas
[hxd, hyd, hzd, hpsid, hxdp, hydp, hzdp, hpsidp] = Trayectorias(3,t,5);

%% GENERALIZED DESIRED SIGNALS
%hd = [hxd; hyd; hzd; hpsid];
hd = [hxd;hyd;hzd;0*hpsid;hxdp; hydp; hzdp; 0*hpsidp];

%hdp = [hxdp;hydp;hzdp;hpsidp];


%% Definicion de los limites de las acciondes de control
bounded = 3*[1.2; -1.2; 1.2; -1.2; 1.2; -1.2; 5.5; -5.5];

%% Definicion del vectro de control inicial del sistema

v_N = zeros(N,4);
H0 = repmat(H,1,N+1)'; 
x_N = H0;

A = [   -1.5105    0.4017    1.7112   -0.1643;
    0.9038   -0.7523   -1.4176   -0.2979;
    -0.7307    0.3331   -0.7422   -0.2671;
    -0.0582    0.1306    0.2789   -2.4667];

B =[    1.9717   -0.5415   -1.5929   -0.0495;
    -0.2361    1.3144    0.8441    0.3255;
    0.7182   -0.2532    0.9526    0.1146;
    -0.5441   -0.4366    0.1335    2.6097];

% Definicion del optimizador
[f, solver, args] = mpc_drone_DMD_gains(A,B,bounded, N, L, ts, gains);

tic
for k=1:length(t)-N

    %% Generacion del; vector de error del sistema
    he(:,k)=hd(1:4,k)-h(:,k);
    
    
    [u_opt,x_opt] = SolverUAV_MPC_din_gains(h,v,hd,N,x_N,v_N,args,solver,k);
    
    
    vref(:,k)= u_opt(1,:)';
    h_N(:,1:4,k) = x_opt(:,1:4);

    
    %% Dinamica del sistema 
    
    x(:,k+1) = UAV_DMD_RK4(A,B,x(:,k),vref(:,k),L,ts);
    
    h(:,k+1) = x(1:4,k+1);
    v(:,k+1) = x(5:8,k+1);
    
    %% Simulacion del sistema
%     h=h+system(h,[ul(k+1);um(k+1);un(k+1);w(k+1)],f,ts);
    

        
    %% Actualizacion de los resultados del optimizador para tener una soluciona aproximada a la optima
    
    v_N = [u_opt(2:end,:);u_opt(end,:)];
    x_N = [x_opt(2:end,:);x_opt(end,:)];
    
end
toc
cost = norm(he,2);% + 0.05*norm(x,1);
%cost = norm(he'*he,2); %;
%cost = he'*he;
end

