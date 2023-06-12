%% test convergency for a given variable
% The variable needs to be a struc variable abtained from the program
% "calcul_solution_alfven_balayage_rayon"
clear all
clc
%% toolbox flottant haute précision
%addpath('C:\Users\lalloz\Documents\Multiprecision Computing Toolbox')
precision= 8*32; % Setup default precision to 40 decimal digits (quadruple).
digits(precision);  

set(0,'defaultTextInterpreter','latex')
opengl software

precis_b= 32;

save_video= 0;
show_fig = 1;

%% ajout des chemins utiles au fonctionnement du code
% addpath(['C:\Users\lalloz\Documents\these\onde alfven\etude theorique numerique\'...
%     'programme_calcul_ratio_ampl_phase_shift']);
addpath(['C:\Users\lalloz\Documents\these\onde alfven\etude theorique numerique\'...
    '1_electrode_model_program']);
addpath(['C:\Users\lalloz\Documents\these\onde alfven\etude theorique numerique\'...
    'common_functions']);
% fit files downloading : allows the user to set values of delta kt and
% nk_kt regarding the sought error (in %)

%For r_investigation <0.2
%addpath(['C:\Users\lalloz\Documents\these\onde alfven\etude theorique numerique\'...
%    'modélisation_b_induit\fichier_fit_approx_b\radius_elec_0p5mm']);

%For r_investigation > 0.2
addpath(['C:\Users\lalloz\Documents\these\onde alfven\etude theorique numerique\'...
    'modélisation_b_induit\fichier_fit_approx_b\radius_elec_0p5mm\r_inv_max_0p8\bis'])
%% fitting files loading
load('fit_r_investigation__delta_kt__to__err_relative_pc.mat')
load('fit_r_investigation__delta_kt__to__nb_kt')
load('fit_r_investigation__err_relative_pc__to__delta_kt')


%% Import data

folder= ['C:\Users\lalloz\Documents\these\onde alfven\etude theorique numerique\'...
    'analyse de convergence\test_convergence_sur_gradV\Ha_38e3\r_0p001to0p6\'...
    'Reta_100__368Hz\'];
subfolder= 'nb_kt_2500_delta_kt_4p8\';
file_name= 'donnees.mat';
file_path = fullfile(folder,subfolder,file_name);
load(file_path)

%% Fin import
nb_mode_max= nb_kt-1;
%nb_mode_test= floor(logspace(log10(1),log10(nb_mode_max),20));
nb_mode_test= floor(linspace(1,nb_mode_max,10));
V_env_mode_struct=struct;
V_phase_mode_struct= struct;
r_loc= 1:20:length(r_investigation);
i=1
%%
for nb_mode_i= nb_mode_test
    V_env_mode= [];
    V_phase_mode=[];
    for r_loc_i= r_loc 
    [V_env_i,V_phase_i]=amplitudes_somme_onde(V_struc(r_loc_i).value(1:nb_mode_i,:),s1_struc(r_loc_i).value(1:nb_mode_i),...
    k1_struc(r_loc_i).value(1:nb_mode_i),s2_struc(r_loc_i).value(1:nb_mode_i),k2_struc(r_loc_i).value(1:nb_mode_i)...
    ,hadm,nb_point_z,precis_b,epsilon);
    V_env_mode= [V_env_mode; V_env_i];
    V_phase_mode= [V_phase_mode; V_phase_i];
    r_loc_i
    end
    V_env_mode_struct(i).value= V_env_mode;
    V_phase_mode_struct(i).value= V_phase_mode;
    fprintf('%4.2f %%\n',100*i/length(nb_mode_test))
    i=i+1;
end
%%
err_V_ampl_L2= ones(length(V_env_mode_struct),1);
for i=1:length(V_env_mode_struct)
    err_V_ampl_L2(i)=  vpa(norm(V_env_mode_struct(i).value- vpa(V_env_mode_struct(end).value,32)),32)/...
       vpa(norm(V_env_mode_struct(end).value));
end

err_V_ampl_Linf= ones(length(V_env_mode_struct),1);
for i=1:length(V_env_mode_struct)
    err_V_ampl_Linf(i)=  norm(V_env_mode_struct(i).value- V_env_mode_struct(end).value,'Inf')/...
       norm(V_env_mode_struct(end).value,'Inf');
end

%% data save
%selpath = uigetdir(['C:\Users\lalloz\Documents\these\onde alfven\etude theorique numerique'...
%    '\1_electrode_model_Results\one plate forcing case']);
save('donnees_conv.mat','V_env_mode_struct','V_phase_mode_struct','nb_mode_test',...
    "err_V_ampl_Linf","err_V_ampl_L2")

%%
figure;
plot(nb_mode_test(1:end),vpa(err_V_ampl_L2),'-+')
hold on
plot(nb_mode_test(1:end),vpa(err_V_ampl_Linf),'-+')
legend('$L^{2} norm$','$L^{\infty} norm$','interpreter','latex')
set(gca, 'YScale', 'log')
xlabel('$N_{\kappa_\perp}$')
ylabel('$\varepsilon$')
f= gca
f.XLabel.FontSize= 14;
f.YLabel.FontSize= 14;
f.TickLabelInterpreter= 'latex';

