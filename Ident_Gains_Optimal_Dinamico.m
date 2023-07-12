%% IDENTIFICATION OF PARAMETERS DORNE DYNAMIC %%

clc; clear all; close all; warning off % Inicializacion
 
ts = 1/10;       % Tiempo de muestreo
tfin = 20;      % Tiempo de simulación
t = 0:ts:tfin;
N = length(t);

load("chi_values.mat");
chi_real(:,1) = chi';

%% Condiciones Iniciales
h(:,1) = [0;0;1;0];
h_p(:,1) = [0;0;0;0];

%% Variables definidas por la TRAYECTORIA y VELOCIDADES deseadas
[xd, yd, zd, psid, xdp, ydp, zdp, psidp] = Trayectorias(3,t);
                                                      
hd = [xd;yd;zd;psid];    
hd_p = [xdp;ydp;zdp;psidp]; 
                                                      
disp('Empieza el programa')

%% Parametros del optimizador
% options = optimset('Display','iter',...
%                 'TolFun', 1e-8,...
%                 'MaxIter', 10000,...
%                 'Algorithm', 'active-set',...
%                 'FinDiffType', 'forward',...
%                 'RelLineSrchBnd', [],...
%                 'RelLineSrchBndDuration', 1,...
%                 'TolConSQP', 1e-6); 
% x0=zeros(8,1);           
% f_obj1 = @(x)  funcion_costo_Ident_Gain(x, hd_p, hd, N, ts);
% x = fmincon(f_obj1,x0,[],[],[],[],[],[],[],options);
% values_final = x;

%%
% OPTIMIZATION PARAMETERS IDENTIFICATION
options = optimoptions(@fmincon, 'Algorithm','interior-point'); 
options.MaxFunctionEvaluations = 10000;   
rng default;
ms = MultiStart('FunctionTolerance',2e-4,'UseParallel',true,'Display','iter', 'MaxTime', 600);

% INITIAL VALUES
chi = ones(16,1);  
f_obj1 = @(x)  funcion_costo_Ident_Gain_Dinamica(x, hd_p, hd, N, ts, chi_real); 
vc_min = -1.5;
vc_max = 1.5;
Delta_hd_p_min = -0.001;
Delta_hd_p_max = 0.001;
problem = createOptimProblem('fmincon','objective',f_obj1,'x0',chi,...
                             'lb',1*[ones(16,1)],'ub',[],...
                             'nonlcon',{@(u)vc_constraint(u,vc_min,vc_max)},...
                             'options',options);

gs = GlobalSearch(ms);
[x, f] = run(gs, problem);
values_final = x;
X = values_final'


%%
% %save("/home/bryansgue/Doctorado/Matlab/UAV/DynamicControllers/chi_values.mat","chi");
% 
% %save("/home/bryan/Proyectos_Matlab/Doctorado/Matlab/IdentificacionM100/IdentificacionUAV/chi_values.mat","chi");
% save("chi_values_fullUAV.mat","chi");
%load("parameters.mat");


%% ACCIONES DE CONTROL TEST 
% F_ref =     0.0*ones(1, length(t));
% phi_ref =   0*ones(1, length(t));
% theta_ref = 0*ones(1, length(t));
% psi_ref =   0*ones(1, length(t));
% u_ref = [0*F_ref;0*F_ref;F_ref; phi_ref; theta_ref; psi_ref];

% SIMULATION DYNAMICS

figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 3]);
luz = light;
luz.Color=[0.65,0.65,0.65];
luz.Style = 'infinite';

%b) Dimenciones del Robot
Drone_Parameters(0.02);
%c) Dibujo del Robot
%G2=Drone_Plot_3D(q_estimate(1,1),q_estimate(2,1),q_estimate(3,1),q_estimate(4,1), q_estimate(5,1), q_estimate(6,1));hold on

%plot3(H(1,1),H(2,1),H(3,11),'--','Color',[56,171,217]/255,'linewidth',1.5);hold on,grid on

close all

axis equal
for k = 1:100:length(t)-1
    drawnow
    delete(G2);
   
    G2=Drone_Plot_3D(q_estimate(1,k),q_estimate(2,k),q_estimate(3,k),q_estimate(4,k), q_estimate(5,k), q_estimate(6,k));hold on

    plot3(q_estimate(1,1:k),q_estimate(2,1:k),q_estimate(3,1:k),'--','Color',[56,171,217]/255,'linewidth',1.5);
    %plot3(obs(1,:),obs(2,:),obs(3,:),'x','Color',[0,171,217]/255,'linewidth',2);
    legend({'$\mathbf{h}$','$\mathbf{h}_{des}$'},'Interpreter','latex','FontSize',11,'Location','northwest','Orientation','horizontal');
    legend('boxoff')
    title('$\textrm{Movement Executed by the Aerial Robot}$','Interpreter','latex','FontSize',11);
    xlabel('$\textrm{X}[m]$','Interpreter','latex','FontSize',9); ylabel('$\textrm{Y}[m]$','Interpreter','latex','FontSize',9);zlabel('$\textrm{Z}[m]$','Interpreter','latex','FontSize',9);
    
end





% 
% Parameters fancy plots
% define plot properties
lw = 2; % linewidth 1
lwV = 2; % linewidth 2
fontsizeLabel = 9; %11
fontsizeLegend = 9;
fontsizeTicks = 9;
fontsizeTitel = 9;
sizeX = 900; % size figure
sizeY = 300; % size figure

% color propreties
C1 = [246 170 141]/255;
C2 = [51 187 238]/255;
C3 = [0 153 136]/255;
C4 = [238 119 51]/255;
C5 = [204 51 17]/255;
C6 = [238 51 119]/255;
C7 = [187 187 187]/255;
C8 = [80 80 80]/255;
C9 = [140 140 140]/255;
C10 = [0 128 255]/255;
C11 = [234 52 89]/255;
C12 = [39 124 252]/255;
C13 = [40 122 125]/255;
%C14 = [86 215 219]/255;
C14 = [252 94 158]/255;
C15 = [244 171 39]/255;
C16 = [100 121 162]/255;
C17 = [255 0 0]/255;
% 
% 
% 
% 
figure('Position', [10 10 sizeX sizeY])
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8.5 11]);
set(gcf, 'PaperPositionMode', 'manual');

subplot(2,2,1)
%plot(t(1:length(F_ref)),F_ref,'-.','Color',C9,'LineWidth',lw*1.2); hold on
%plot(uv(1,:),uv(2,:),'-','Color',C11,'LineWidth',lw);
plot(t,hx,'-','Color',C11,'LineWidth',lw); hold on;
plot(t,q_estimate(1,1:length(t)),'--','Color',C12,'LineWidth',lw);
grid minor;
grid minor;
set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
%xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(a)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$hx_{ref}$','$hx_{m}$'},'interpreter','latex','fontsize',fontsizeLegend)
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
grid minor;


subplot(2,2,2)
%plot(t(1:length(phi_ref)),phi_ref,'-.','Color',C9,'LineWidth',lw*1.2); hold on
plot(t,hy,'-','Color',C13,'LineWidth',lw); hold on;
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
plot(t,q_estimate(2,1:length(t)),'--','Color',C14,'LineWidth',lw);
grid minor;

set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
%xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(b)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$hy_{ref}$','$hy_{m}$'},'interpreter','latex','fontsize',fontsizeLegend)

subplot(2,2,3)
%plot(t(1:length(theta_ref)),theta_ref,'-.','Color',C9,'LineWidth',lw*1.2); hold on
plot(t,hz,'-','Color',C2,'LineWidth',lw); hold on;
plot(t,q_estimate(3,1:length(t)),'--','Color',C15,'LineWidth',lw);
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
grid minor;

set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(c)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$hz_{ref}$','$hz_{m}$'},'interpreter','latex','fontsize',fontsizeLegend)



figure('Position', [10 10 sizeX sizeY])
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8.5 11]);
set(gcf, 'PaperPositionMode', 'manual');

subplot(2,2,1)
%plot(t(1:length(F_ref)),F_ref,'-.','Color',C9,'LineWidth',lw*1.2); hold on
%plot(uv(1,:),uv(2,:),'-','Color',C11,'LineWidth',lw);
plot(t,q_p(1,:),'-','Color',C11,'LineWidth',lw); hold on;
plot(t,q_estimate(7,1:length(t)),'--','Color',C12,'LineWidth',lw);
grid minor;
grid minor;
set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
%xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(a)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$xp_{ref}$','$xp_{m}$'},'interpreter','latex','fontsize',fontsizeLegend)
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
grid minor;


subplot(2,2,2)
%plot(t(1:length(phi_ref)),phi_ref,'-.','Color',C9,'LineWidth',lw*1.2); hold on
plot(t,q_p(2,:),'-','Color',C13,'LineWidth',lw); hold on;
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
plot(t,q_estimate(8,1:length(t)),'--','Color',C14,'LineWidth',lw);
grid minor;

set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
%xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(b)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$yp_{ref}$','$yp_{m}$'},'interpreter','latex','fontsize',fontsizeLegend)

subplot(2,2,3)
%plot(t(1:length(theta_ref)),theta_ref,'-.','Color',C9,'LineWidth',lw*1.2); hold on
plot(t,q_p(3,:),'-','Color',C2,'LineWidth',lw); hold on;
plot(t,q_estimate(9,1:length(t)),'--','Color',C15,'LineWidth',lw);
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
grid minor;

set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(c)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$zp_{ref}$','$zp_{m}$'},'interpreter','latex','fontsize',fontsizeLegend)


figure('Position', [10 10 sizeX sizeY])
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8.5 11]);
set(gcf, 'PaperPositionMode', 'manual');

subplot(2,2,1)
%plot(t(1:length(F_ref)),F_ref,'-.','Color',C9,'LineWidth',lw*1.2); hold on
%plot(uv(1,:),uv(2,:),'-','Color',C11,'LineWidth',lw);
plot(t,q_p(4,:),'-','Color',C11,'LineWidth',lw); hold on;
plot(t,q_estimate(10,1:length(t)),'--','Color',C12,'LineWidth',lw);
grid minor;
grid minor;
set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
%xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(a)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$\phi_{pref}$','$\phi_{pm}$'},'interpreter','latex','fontsize',fontsizeLegend)
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
grid minor;


subplot(2,2,2)
%plot(t(1:length(phi_ref)),phi_ref,'-.','Color',C9,'LineWidth',lw*1.2); hold on
plot(t,q_p(5,:),'-','Color',C13,'LineWidth',lw); hold on;
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
plot(t,q_estimate(11,1:length(t)),'--','Color',C14,'LineWidth',lw);
grid minor;

set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
%xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(b)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$\theta_{pref}$','$\theta_{pm}$'},'interpreter','latex','fontsize',fontsizeLegend)

subplot(2,2,3)
%plot(t(1:length(theta_ref)),theta_ref,'-.','Color',C9,'LineWidth',lw*1.2); hold on
plot(t,q_p(6,:),'-','Color',C2,'LineWidth',lw); hold on;
plot(t,q_estimate(12,1:length(t)),'--','Color',C15,'LineWidth',lw);
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
grid minor;

set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(c)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$\psi_{pref}$','$\psi_{pm}$'},'interpreter','latex','fontsize',fontsizeLegend)

