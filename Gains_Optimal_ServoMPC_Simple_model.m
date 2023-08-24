%% IDENTIFICATION OF PARAMETERS DORNE DYNAMIC %%

clc; clear all; close all; warning off % Inicializacion
 
load("chi_simple.mat");
chi_uav = chi';

%% DEFINITION OF TIME VARIABLES
f = 30 % Hz 
ts = 1/f;
to = 0;
tf = 15;
t = (to:ts:tf);

%% Definicion del horizonte de prediccion
N = 10; %20 
                                                      
disp('Empieza el programa')

%%
% OPTIMIZATION PARAMETERS IDENTIFICATION
% options = optimoptions(@fmincon, 'Algorithm','interior-point'); 
% options.MaxFunctionEvaluations = 10000;   
% rng default;
% ms = MultiStart('FunctionTolerance',2e-4,'UseParallel',true,'Display','iter', 'MaxTime', 600);
% 
% % INITIAL VALUES
% chi = ones(8,1);  
% 
% f_obj1 = @(gains)  funcion_costo_S_simple_MPC(gains, t, ts, N, chi_uav); 
% vc_min = -5.0;
% vc_max = 5;
% Delta_hd_p_min = -0.001;
% Delta_hd_p_max = 0.001;
% % problem = createOptimProblem('fmincon','objective',f_obj1,'x0',chi,...
% %                              'lb',1*[ones(8,1)],'ub',[],...
% %                              'nonlcon',{@(u)vc_constraint(u,vc_min,vc_max)},...
% %                              'options',options);
% %                          
% problem = createOptimProblem('fmincon','objective',f_obj1,'x0',chi,...
%                              'lb',[],'ub',[],...
%                              'nonlcon',{},...
%                              'options',options);
% 
% gs = GlobalSearch(ms);
% [x, f] = run(gs, problem);
% values_final = x;
% X = values_final'

options = optimset('Display','iter',...
    'TolFun', 1e-8,...
    'MaxIter', 60000,...
    'MaxFunEvals', 10000,... % Agregado nuevo_valor aquí
    'Algorithm', 'active-set',...
    'FinDiffType', 'forward',...
    'RelLineSrchBnd', [],...
    'RelLineSrchBndDuration', 1,...
    'TolConSQP', 2e-8);

x0 = ones(1,8);

f_obj1 = @(gains)  funcion_costo_S_simple_MPC(gains, t, ts, N, chi_uav); 
                            
chi = fmincon(f_obj1,x0,[],[],[],[],[],[],[],options);



