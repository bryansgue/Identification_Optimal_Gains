%%******************************************************************************************************************
%************************************ SEGUIMIENTO DE TRAYECTORIA **************************************************
%************************************* ROBOT MANIPULADOR AÉREO *****************************************************
%******************************************************************************************************************
clc; clear all; close all; warning off % Inicializacion

ts = 1/30;       % Tiempo de muestreo
tfin = 45;      % Tiempo de simulación
t = 0:ts:tfin;
load("chi_values.mat");
a=0;
b=0;
L=[a;b];

chi_real(:,1) = chi';
%% Variables definidas por la TRAYECTORIA y VELOCIDADES deseadas
[xd, yd, zd, psid, xdp, ydp, zdp, psidp] = Trayectorias(3,t);
%% GENERALIZED DESIRED SIGNALS
hd = [xd; yd; zd; psid];
hd_p = [xdp;ydp;zdp; psidp];

%% a) Posiciones iniciales del UAV
xu(1) = 0; 
yu(1) = 0; 
zu(1) = 1; 
psi(1)= 0;
h=[xu(1);yu(1);zu(1);psi(1)]
%% Velocidad inicial real del UAV
u = [0;0;0;0];
u_realchi = [0;0;0;0];
%% Ganancia Compensacion Dinamica
K = 1;
G = 1;
%******************************************************************************************************************
%***************************************** CONTROLADOR ***********************************************************
%*****************************************************************************************************************
disp('Empieza el programa')
uc_p(:,1) = [0;0;0;0];

X = [5.0000
    1.0032
    5.0000
    4.6478
    1.1463
    1.3356
    1.1794
    1.0002
    4.9999
    2.1339
    5.0000
    4.9988
    4.9986
    4.9973
    5.0000
    4.9937];


for k=1:length(t)
tic

%% 1) LEY DE CONTROL
he(:,k) = hd(:,k)-h(:,k);
uc(:,k) = Control_UAV_s(X,hd_p(:,k),hd(:,k),h(:,k));
%% 2) ACELERATIONS VCP 
if k>1
uc_p(:,k)=(uc(:,k)- uc(:,k-1))/ts;    
else
uc_p(:,k)=0;    
end
%% DYNAMIC COMPENSATION
uref(:,k) = dynamic_compensation_UAV_s(X, uc_p(:,k), uc(:,k), u(:,k), chi_real, L, ts);
%% 2) DINAMICA DEL UAV (VELOCIDAD Y POSICION)
u(:, k+1) = system_dynamic_s(chi_real, u(:,k), uref(:,k), h(:,k), L,ts);
h(:,k+1) = RK4_UAV_simple(h(:,k),u(:,k),ts);

%% 3) Tiempo de máquina   
dt(k) = toc;
end
disp('Fin de los calculos')

%*************************************************************************%
%**************ANIMACION SEGUIMIENTO DE TRAYECTORIA **********************%
%% ***********************************************************************%
disp('Animacion RUN')

% 1) Parámetros del cuadro de animacion
figure(1)
axis equal
view(0,0) % Angulo de vista
cameratoolbar
title ("Simulacion")

% 2) Configura escala y color del UAV
Drone_Parameters(0.02);
H1 = Drone_Plot_3D(xu(1),yu(1),zu(1),0,0,psi(1));hold on


% c) Gráfica de la trayectoria deseada
plot3(xd,yd,zd,'--')
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');

% 5) Simulación de movimiento del manipulador aéreo
for k=1:50:length(t)  
% a) Eliminas los dibujos anteriores del manipulador aéreo
delete(H1);
H1 = Drone_Plot_3D(h(1,k),h(2,k),h(3,k),0,0,h(4,k)); hold on
% b) Gráfica la posición deseada vs actual en cada frame
plot3(h(1,1:k),h(2,1:k),h(3,1:k),'r')
hold on
plot3(hd(1,1:k),hd(2,1:k),hd(3,1:k),'r')

pause(0.1)
end

disp('FIN Simulación RUN')  

%%
%******************************************************************************************************************
%********************************************* GR�?FICAS ***********************************************************
%% ****************************************************************************************************************


% 2) Cálculos del Error
figure(2)
plot(he(1,:)); hold on, grid on
plot(he(2,:)); hold on
plot(he(3,:)); hold on
plot(he(4,:)); hold on
legend("hxe","hye","hze","psie")
title ("Errores de posición")

% 3) Posiciones deseadas vs posiciones reales del extremo operativo del manipulador aéreo
figure(3)

subplot(4,1,1)
plot(hd(1,:))
hold on
plot(h(1,:))
legend("xd","hx")
ylabel('x [m]'); xlabel('s [ms]');
title ("Posiciones deseadas y reales del extremo operativo del manipulador aéreo")

subplot(4,1,2)
plot(hd(2,:))
hold on
plot(h(2,:))
legend("yd","hy")
ylabel('y [m]'); xlabel('s [ms]');

subplot(4,1,3)
plot(hd(3,:))
hold on
plot(h(3,:))
grid on
legend("zd","hz")
ylabel('z [m]'); xlabel('s [ms]');

subplot(4,1,4)
plot(hd(4,:))
hold on
plot(h(4,:))
legend("psid","psi")
ylabel('psi [rad]'); xlabel('s [ms]');

% 3) Posiciones deseadas vs posiciones reales del extremo operativo del manipulador aéreo
figure(4)


plot(uc(1,1:end))
hold on
plot(u(1,1:end))
hold on
plot(uref(1,1:end))
legend("ulc","ul","ul_{ref}")
ylabel('x [m/s]'); xlabel('s [ms]');
title ("Posiciones deseadas y reales del extremo operativo del manipulador aéreo")

figure(5)
plot(uc(2,1:end))
hold on
plot(u(2,1:end))
hold on
plot(uref(2,1:end))
legend("umc","um","um_{ref}")
ylabel('y [m/s]'); xlabel('s [ms]');

figure(6)
plot(uc(3,1:end))
hold on
plot(u(3,1:end))
hold on
plot(uref(3,1:end))
legend("unc","un","un_{ref}")
ylabel('z [m/ms]'); xlabel('s [ms]');

figure(7)
plot(uc(4,1:end))
hold on
plot(u(4,1:end))
hold on
plot(uref(4,1:end))
legend("wc","w","w_{ref}")
ylabel('psi [rad/s]'); xlabel('s [ms]');
  