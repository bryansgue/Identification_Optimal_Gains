 %% ESTIMATION FO PARAMETERS DORNE DYNAMIC %%

% clear variables
clc, clear all, close all;

% LOAD VALUES FROM MATRICES
load('M100_ident.mat')
clear tf;
desface = 200;
t = t(1,1:end-desface);
N = length(t);

% REFERENCE SIGNALS FOR WEBOTS
ul_ref = v_ref(1,1:end-desface);
um_ref = v_ref(2,1:end-desface);
un_ref = v_ref(3,1:end-desface);
w_ref = v_ref(4,1:end-desface);

% REAL SYSTEM VELICITIES FOR WEBOTS DATA
ul = double(v(1,1:length(ul_ref)));
um = double(v(2,1:length(um_ref)));
un = double(v(3,1:length(un_ref)));

roll = double(euler(1,1:length(w_ref)));
pitch = double(euler(2,1:length(w_ref)));
yaw = double(euler(3,1:length(w_ref)));
euler_p = euler_p(:,1:end-1);
% REFERENCE SIGNALS
v_ref = [ul_ref; um_ref; un_ref; w_ref];

% REAL SYSTEM VELICITIES
v_real = [ul; um; un; roll;pitch; yaw];

% REAL SYSTEM ACCELERATIONS
ulp_real = [0, diff(ul)/ts];
ump_real = [0 , diff(um)/ts];
unp_real = [0 , diff(un)/ts];

roll_p = euler_p(1,1:length(ul_ref));
pitch_p = euler_p(2,1:length(ul_ref));
yaw_p = euler_p(3,1:length(ul_ref));

vp_real = [ulp_real; ump_real; unp_real; roll_p ; pitch_p; yaw_p];
% calculo propi

tic 
Gamma = v_ref;
Omega = [v_real(:,1:end-1); 
         v_ref(:,1:end-1)];
[U,S,V] = svd(Omega,'econ') ;


G = vp_real(:,2:end)*V*S^(-1)*U';

U1=U(1:6,1:10);
U2=U(7:10,1:10);


%
A = G(:,1:6);
B = G(:,6+1:end);


A_ = vp_real(:,2:end)*V*S^(-1)*U1';
B_ = vp_real(:,2:end)*V*S^(-1)*U2';

%A = U'* vp_real(:,2:end)*V*S^(-1)*U1'*U
%B = U* vp_real(:,2:end)'*V'*S^(-1)*U2'*U
%
%[A,B] = Ident(v_real,v_ref,vp_real,ts);
toc
% 
%% SIMULATION DYNAMICS
v_estimate = v_real(:,1);
for k=1:length(t)
    a = (A*v_estimate(:,k)+B*v_ref(:,k));
    v_estimate(:, k+1) = v_estimate(:, k) + a*ts;
%     v_estimate(:, k+1) = sysmodel_DMDc.A*v_estimate(:,k)+sysmodel_DMDc.B*vref(:,k);
end
% 
% save("IdentDMD_test.mat","v_estimate","v_ref","v_real","t");
% save("A_B_values.mat","A","B");
%% Parameters fancy plots
% define plot properties
lw = 2; % linewidth 1
lwV = 2; % linewidth 2
fontsizeLabel = 9 ; %11
fontsizeLegend = 20;
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
C14 = [86 215 219]/255;
C14 = [252 94 158]/255;
C15 = [244 171 39]/255;
C16 = [100 121 162]/255;
C17 = [255 0 0]/255;


figure('Position', [10 10 sizeX sizeY])
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8.5 11]);
set(gcf, 'PaperPositionMode', 'manual');

subplot(3,1,1)
plot(t(1:length(ul_ref)),ul_ref,'-.','Color',C9,'LineWidth',lw*1.2); hold on
%plot(uv(1,:),uv(2,:),'-','Color',C11,'LineWidth',lw);
plot(t,ul,'-','Color',C11,'LineWidth',lw);
plot(t,v_estimate(1,1:length(t)),'--','Color',C12,'LineWidth',lw);
grid minor;
grid minor;
set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
%xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(a)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$\mu_{lref}$','$\mu_{l}$','$\mu_{lm}$'},'interpreter','latex','fontsize',fontsizeLegend)
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
grid minor;


subplot(3,1,2)
plot(t(1:length(um_ref)),um_ref,'-.','Color',C9,'LineWidth',lw*1.2); hold on
plot(t,um,'-','Color',C13,'LineWidth',lw);
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
plot(t,v_estimate(2,1:length(t)),'--','Color',C14,'LineWidth',lw);
grid minor;

set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
%xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(b)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$\mu_{mref}$','$\mu_{m}$','$\mu_{mm}$'},'interpreter','latex','fontsize',fontsizeLegend)

subplot(3,1,3)
plot(t(1:length(un_ref)),un_ref,'-.','Color',C9,'LineWidth',lw*1.2); hold on
plot(t,un,'-','Color',C2,'LineWidth',lw);
plot(t,v_estimate(3,1:length(t)),'--','Color',C15,'LineWidth',lw);
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
grid minor;

set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(c)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$\mu_{nref}$','$\mu_{n}$','$\mu_{nm}$'},'interpreter','latex','fontsize',fontsizeLegend)

%%%%
figure(2)

subplot(3,1,1)
%plot(t(1:length(ul_ref)),ul_ref,'-.','Color',C9,'LineWidth',lw*1.2); hold on
%plot(uv(1,:),uv(2,:),'-','Color',C11,'LineWidth',lw);
plot(t,roll,'-','Color',C11,'LineWidth',lw);hold on
plot(t,v_estimate(4,1:length(t)),'--','Color',C12,'LineWidth',lw);
grid minor;
grid minor;
set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
%xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(a)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$\mu_{lref}$','$\mu_{l}$','$\mu_{lm}$'},'interpreter','latex','fontsize',fontsizeLegend)
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
grid minor;


subplot(3,1,2)

plot(t,pitch,'-','Color',C13,'LineWidth',lw);hold on
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
plot(t,v_estimate(5,1:length(t)),'--','Color',C14,'LineWidth',lw);
grid minor;

set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
%xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(b)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$\mu_{mref}$','$\mu_{m}$','$\mu_{mm}$'},'interpreter','latex','fontsize',fontsizeLegend)

subplot(3,1,3)
plot(t,yaw,'-','Color',C2,'LineWidth',lw);hold on
plot(t,v_estimate(6,1:length(t)),'--','Color',C15,'LineWidth',lw);
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
grid minor;

set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(c)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$\mu_{nref}$','$\mu_{n}$','$\mu_{nm}$'},'interpreter','latex','fontsize',fontsizeLegend)



print -dpng Model_dmd_identification
print -depsc Model_dmd_identification

