%******************************************************************************************************************
%************************************ SEGUIMIENTO DE TRAYECTORIA **************************************************
%************************************* UAV *****************************************************
%******************************************************************************************************************
clc; clear all; close all; warning off % Inicializacion
 
ts = 1/30;       % Tiempo de muestreo
tfin = 60;      % Tiempo de simulación
t = 0:ts:tfin;

%% Condiciones Iniciales
h(:,1) = [0;0;1;0];
h_p(:,1) = [0;0;0;0];

%% Variables definidas por la TRAYECTORIA y VELOCIDADES deseadas
[xd, yd, zd, psid, xdp, ydp, zdp, psidp] = Trayectorias(3,t);
                                                      
hd = [xd;yd;zd;psid]; 
hd_p = [xdp;ydp;zdp;psidp];       
                                                      
disp('Empieza el programa')
X = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
%******************************************************************************************************************
%***************************************** CONTROLADOR ***********************************************************
%*****************************************************************************************************************

for k=1:length(t)
    tic  
%% 1) Controlador    
he(:,k) = hd(:,k)-h(:,k);
uc(:,k) = Control_UAV_s(X,hd_p(:,k),hd(:,k),h(:,k)); 
           
h(:,k+1) = RK4_UAV_simple(h(:,k),uc(:,k),ts);
    
%% 4) Tiempo de muestreo   
     dt(k) = toc;
     
end
disp('Fin de los cálculos')
%******************************************************************************************************************
%********************************* ANIMACIÓN SEGUIMIENTO DE TRAYECTORIA ******************************************
%% ******************************************************************************************************************
disp('Simulación RUN')

% 1) Parámetros del cuadro de animacion
  figure(1)
  axis equal
  view(-15,15) % Angulo de vista
  cameratoolbar
  title ("Simulación")
  
% 2) Configura escala y color del UAV
 Drone_Parameters(0.02);
 H1 = Drone_Plot_3D(h(1,1),h(2,1),h(3,1),0,0,h(4,1));hold on
 
 

% 3) Dimensiones del brazo robótico
  

% 4) Dibujo inicial del manipulador aéreo

  % a) Gráfica inicial del brazo
  
  
  % b) Gráfica inicial del UAV  

  
  % c) Gráfica de la trayectoria deseada
  p1 = plot3(xd,yd,zd,'--');
  p2 = plot3(0,0,0,'--');
  
  xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
  
% 5) Simulación de movimiento del manipulador aéreo
  for k=1:10:length(t)  
  % a) Eliminas los dibujos anteriores del manipulador aéreo
    delete(H1);
    delete(p1);
    delete(p2);
    H1 = Drone_Plot_3D(h(1,k),h(2,k),h(3,k),0,0,h(4,k)); hold on
  % b) Gráfica la posición deseada vs actual en cada frame
    p1 = plot3(h(1,1:k),h(2,1:k),h(3,1:k),'r');
    hold on
    p2 = plot3(hd(1,1:k),hd(2,1:k),hd(3,1:k),'b');
    
  % c) Grafica de brazo robótico
    
    
  % d) Gráfica del UAV
    
  
  pause(0.1)
  end
%%
%******************************************************************************************************************
%********************************************* GR�?FICAS ***********************************************************
%% ****************************************************************************************************************

% 1) Igualar columnas de los vectores creados
%   xu(:,end)=[];
%   yu(:,end)=[];
%   zu(:,end)=[];
%   psi(:,end)=[];
  
% 2) Cálculos del Error
  figure(2)
  hxe= he(1,1:end);
  hye= he(2,1:end);
  hze= he(3,1:end);
  psie= he(4,1:end);
  plot(hxe), hold on, grid on
  plot(hye)
  plot(hze)
  plot(psie)
  legend("hxe","hye","hze","psie")
  title ("Errores de posición")
  
% 3) Posiciones deseadas vs posiciones reales del extremo operativo del manipulador aéreo
  figure(3)
  
  subplot(4,1,1)
  plot(hd(1,1:end))
  hold on
  plot(h(1,1:end))
  legend("xd","hx")
  ylabel('x [m]'); xlabel('s [ms]');
  title ("Posiciones deseadas y reales del extremo operativo del manipulador aéreo")
  
  subplot(4,1,2)
  plot(hd(2,1:end))
  hold on
  plot(h(2,1:end))
  legend("yd","hy")
  ylabel('y [m]'); xlabel('s [ms]');

  subplot(4,1,3)
  plot(hd(3,1:end))
  hold on
  plot(h(3,1:end))
  grid on
  legend("zd","hz")
  ylabel('z [m]'); xlabel('s [ms]');

  subplot(4,1,4)
  plot(hd(4,1:end))
  hold on
  plot(h(4,1:end))
  legend("psid","psi")
  ylabel('psi [rad]'); xlabel('s [ms]');
  
    
% 3) Posiciones deseadas vs posiciones reales del extremo operativo del manipulador aéreo
  figure(4)
  
  subplot(4,1,1)
  plot(uc(1,1:end))
  hold on
  
  legend("xd","hx")
  ylabel('x [m]'); xlabel('s [ms]');
  title ("Posiciones deseadas y reales del extremo operativo del manipulador aéreo")
  
  subplot(4,1,2)
  plot(uc(2,1:end))
  legend("yd","hy")
  ylabel('y [m]'); xlabel('s [ms]');

  subplot(4,1,3)
  plot(uc(3,1:end))
  grid on
  legend("zd","hz")
  ylabel('z [m]'); xlabel('s [ms]');

  subplot(4,1,4)
  plot(uc(4,1:end))
  legend("psid","psi")
  ylabel('psi [rad]'); xlabel('s [ms]');