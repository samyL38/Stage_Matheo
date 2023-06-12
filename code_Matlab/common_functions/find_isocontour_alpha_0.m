%%%
% Study of the evolution of the isocontour attenuation coeffficent alpha =0
% regarding the location r, for different Reta (i.e. Rnu) and Ha
%
% Remarks :
% - all data are obtained from the "calcul_solution_alfven_balayage_rayon"
% program

%%
clear all
clc
%% toolbox flottant haute précision
%addpath('C:\Users\lalloz\Documents\Multiprecision Computing Toolbox')
precision= 8*32; % Setup default precision to 40 decimal digits (quadruple).
digits(precision);  

set(0,'defaultTextInterpreter','latex')
opengl software

precis_b= 16;

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

%% Aprameter initialisation
alpha_z1_struct=struct;
D_r_alpha_z1_struct=struct;
Reta_vect=[];
Rnu_vect=[];
Ha= [];
r_investigation_struct= struct;
    %% import
load(['donnees_nb_z_2.mat'],'epsilon','gradV_env','gradV_phase',...
    'r_investigation','Ha_num','Rnu','Reta')

%% Parameter                

%Paramètres physiques
mu_0= vpa(1.25*sym(1e-6));
sigma_gal= vpa(3.46*sym(1e6));
visco_gal= vpa(0.0024); %Viscosité dynamique du Galstn en Pa.s
eta= 1/(sigma_gal * mu_0); %diffusivité mhd
rho_gal= vpa(6440);% masse volumique Galinstan kg/m3
nu_gal= visco_gal/rho_gal;

% Grandeurs caractéristiques
h= vpa(0.1); % distance entre la plaque inf et sup en mètre
L_carac=h;%(eta/omega_0)^(1/2);%longueur caractéristique ==> hauteur de la cuve
tau_eta= L_carac^2/eta;
tau_nu= L_carac^2/nu_gal;
I0=3.18; % Injection courant en Ampère
j0= I0/(L_carac)^2; %%% forme de jz = j0/(pi*sigma_r^2)*exp(-(r/sigma_r)^2)
b0= mu_0*I0/(L_carac); %%% forme de b_theta= b0*(1-exp(-((r)/sigma_r).^2))./(r);
E0= eta/L_carac * b0;


%Reta_i= epsilon;
%Rnu_i= tau_nu*2*vpa(pi)*frequence_forcage;
%frequence_forcage= Reta_i/(tau_eta*2*vpa(pi));

%epsilon= Reta_i;

%% attenuation coef at z=1

%initialisation
ind= 11;

%
Ha(ind)= Ha_num;
Reta_vect(ind)= Reta;
Rnu_vect(ind)= Rnu;
alpha_z1_struct(ind).value= log(transpose(gradV_env(:,end)./gradV_env(:,1)));
D_r_alpha_z1_struct(ind).value= gradient(alpha_z1_struct(ind).value,double(r_investigation));
r_investigation_struct(ind).value= double(r_investigation);


%% figure
lgd_label= {};
figure('Name','attenuation coef vs r z=1');
for i=10
plot(r_investigation_struct(i).value,alpha_z1_struct(i).value,'-*')
hold on
%lgd_label= {lgd_label, sprintf('$R_\\eta = %4.2f$',Reta(i))};
end  
%legend(lgd_label,'interpreter','latex')
xlim([0,0.2])

f= gca;
f.TickLabelInterpreter= 'latex';
xlabel('$r$')
ylabel('$\alpha(r,z=1)$')
f.FontSize= 12;      
f.XLabel.FontSize= 14;
f.YLabel.FontSize= 14;
title(sprintf(['Evolution of $\\alpha \\left(r, z=1\\right)$ \n $\\left(Ha= %4.2e\\right)$'],Ha(1)))
%%
figure('Name','Grad_r attenuation coef vs r z=1');
for i=1:4
plot(r_investigation,D_r_alpha_z1(i,:))
hold on
%lgd_label= {lgd_label, sprintf('$R_\\eta = %4.2f$',Reta(i))};

end  
%legend(lgd_label,'interpreter','latex')
f= gca;
f.TickLabelInterpreter= 'latex';
xlabel('$r$')
ylabel('$\alpha(r,z=1)$')
f.FontSize= 12;      
f.XLabel.FontSize= 14;
f.YLabel.FontSize= 14;
title(sprintf(['Evolution of $\\frac{\\partial}{\\partial\\,r} \\alpha \\left( r, z=1\\right)$ \n $Ha= %4.2e,\\, R_{\\eta}/2\\pi = %4.2f,',...
           ' \\,R_\\nu/2\\pi = %4.2e$'],Ha_num,Reta_vect(i)/(2*vpa(pi)),Rnu_vect(i)/(2*vpa(pi))))


%% calcul  root of alpha, maximum of Grad_r alpha
r_investigation_root= zeros(length(Reta_vect),1);
r_investigation_max= zeros(length(Reta_vect),1);
for i=1:length(Rnu_vect)
    ind_root_i= find(alpha_z1_struct(i).value>=0,1);
    r_investigation_root(i)= r_investigation_struct(i).value(ind_root_i);
end


for i=1:length(Rnu_vect)
    [~,ind_max]=max(D_r_alpha_z1_struct(i).value);
    r_investigation_max(i)= r_investigation_struct(i).value(ind_max);
end
%%
scaled_r_root= r_investigation_root./sqrt(Rnu_vect');


%% Saving

save("donnees_scaling_with_Ha.mat",'Ha','Reta_vect','Rnu_vect','alpha_z1','D_r_alpha_z1',...
'r_investigation_mat','r_investigation_max','r_investigation_root')
%%
% Figure vs Reta

%linear
figure('Name','root of attenuation coef vs Reta');
plot(Reta_vect(1:end),r_investigation_root(1:end),'+')
hold on
%lgd_label= {lgd_label, sprintf('$R_\\eta = %4.2f$',Reta(i))};
%legend(lgd_label,'interpreter','latex')
f= gca;
f.TickLabelInterpreter= 'latex';
xlabel('$R_\eta$')
ylabel('$r_{root}$')
f.FontSize= 12;      
f.XLabel.FontSize= 14;
f.YLabel.FontSize= 14;
Reta_i= epsilon;
frequence_forcage= Reta_i/(tau_eta*2*vpa(pi));title(sprintf(['Evolution of the radial location of $\\alpha \\left( r, z=1\\right) = 0$ \n ($Ha= %4.2e$)'],Ha_num))


figure('Name','root of attenuation coef vs Rnu');
plot(Rnu_vect(1:end),r_investigation_root(1:end),'+')
hold on
%lgd_label= {lgd_label, sprintf('$R_\\eta = %4.2f$',Reta(i))};
%legend(lgd_label,'interpreter','latex')
f= gca;
f.TickLabelInterpreter= 'latex';
xlabel('$R_\nu$')
ylabel('$r_{root}$')
f.FontSize= 12;      
f.XLabel.FontSize= 14;
f.YLabel.FontSize= 14;
title(sprintf(['Evolution of the radial location of $\\alpha \\left( r, z=1\\right) = 0$ \n ($Ha= %4.2e$)'],Ha_num))

%%
% log
figure('Name','root of attenuation coef vs Reta');
plot(log10(Reta_vect(1:end)),log10(r_investigation_root(1:end)),'+')
hold on
%lgd_label= {lgd_label, sprintf('$R_\\eta = %4.2f$',Reta(i))};

%legend(lgd_label,'interpreter','latex')
f= gca;
f.TickLabelInterpreter= 'latex';
xlabel('$\log\left(R_\eta\right)$')
ylabel('$\log\left(r_{root}\right)$')
f.FontSize= 12;      
f.XLabel.FontSize= 14;
f.YLabel.FontSize= 14;
title(sprintf(['Evolution of the radial location of $\\alpha \\left( r, z=1\\right) = 0$ \n ($Ha= %4.2e$)'],Ha_num))

%%
figure('Name','root of attenuation coef vs Rnu');
plot((Rnu_vect(1:end)),(scaled_r_root(1:end)),'+')
hold on
%lgd_label= {lgd_label, sprintf('$R_\\eta = %4.2f$',Reta(i))};

%legend(lgd_label,'interpreter','latex')
f= gca;
f.TickLabelInterpreter= 'latex';
xlabel('$\log\left(R_\nu\right)$')
ylabel('$\log\left(r_{root}\right)$')
f.FontSize= 12;      
f.XLabel.FontSize= 14;
f.YLabel.FontSize= 14;
title(sprintf(['Evolution of the radial location of $\\alpha \\left( r, z=1\\right) = 0$ \n ($Ha= %4.2e$)'],Ha_num))


%% Figure vs Ha

%linear
figure('Name','root of attenuation coef vs Ha');
plot(Ha(1:end),r_investigation_root(1:end),'+')
hold on
%lgd_label= {lgd_label, sprintf('$R_\\eta = %4.2f$',Reta(i))};
%legend(lgd_label,'interpreter','latex')
f= gca;
f.TickLabelInterpreter= 'latex';
xlabel('$Ha$')
ylabel('$r_{root}$')
f.FontSize= 12;      
f.XLabel.FontSize= 14;
f.YLabel.FontSize= 14;

title(sprintf(['Evolution of the radial location of $\\alpha \\left( r, z=1\\right) = 0$ \n ($Ha= %4.2e$)'],Ha_num))



%% log
figure('Name','root of attenuation coef vs Ha');
plot(log10(Ha(1:end)),log10(r_investigation_root(1:end)),'+')
f= gca;
f.TickLabelInterpreter= 'latex';
xlabel('$\log\left(Ha\right)$')
ylabel('$\log\left(r_{root}\right)$')
f.FontSize= 12;      
f.XLabel.FontSize= 14;
f.YLabel.FontSize= 14;
title(sprintf(['Radial location of $\\alpha \\left( r, z=1\\right)'...
    '= 0$ versus $Ha$ \n ($R_\\eta= %4.2f$, $R_\\nu= %4.2e$)'],Reta_vect(1),Rnu_vect(1)))

