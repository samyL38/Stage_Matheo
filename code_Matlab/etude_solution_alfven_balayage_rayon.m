clear all
clc
%close all


% Ce programme fournit les solutions complètes pour les perturbations
% magnetique (b) et mécanique (u) azimutales, le gradient de potentiel 
% \partial_r phi, le champs electrique radial Er ainsi que les densités de
% radiale et axial (resp. jr et jz) à partir des équations MHD résistives
% linéarisées. 

% Le forçage de l'onde est effectué en forçant magnetiquement l'ondes au
% niveau de la plaque inférieure, à la pulsation omega

% Ce programme fournit les solutions complètes pour les perturbations
% magnetique (b) et mécanique (u) azimutales, le gradient de potentiel 
% \partial_r phi, le champs electrique radial Er ainsi que les densités de
% radiale et axial (resp. jr et jz)

% Adimensionnement de :
% - jr/jz par I0/L_carac^2
% - b par b0= mu_0*I0/(L_carac)
% - u par u0= eta/L_carac* b0/B0
% - E_0= eta/L_carac * b0
%( - TF_b_theta_kt par L_carac^2*b0)

% Les calculs sont effectués dans le plan (r,z) ==> permet d'extraire des
% cartographies dans ce même plan pour les variables précédemment citées
%%%%
%% toolbox flottant haute précision
%addpath('C:\Users\lalloz\Documents\Multiprecision Computing Toolbox')
precision= 10*32; % Setup default precision to 40 decimal digits (quadruple).
digits(precision);  

set(0,'defaultTextInterpreter','latex')
opengl software

precis_b= 32;

save_video= 0;
show_fig = 1;


%% ajout des chemins utiles au fonctionnement du code
% addpath(['C:\Users\lalloz\Documents\these\onde alfven\etude theorique numerique\'...
%     'programme_calcul_ratio_ampl_phase_shift']);
path_model =[]; % you can write here your own sought path
if ~exist('path_model') || isempty(path_model)
    path_model= uigetdir([],'Select the folder which contains the model scripts');
    path_model= fullfile(path_model,'\');
end
cd(path_model)
addpath(path_model)
path_com_function= [];
i=0;
while isempty(path_com_function) || i == 10
    cd ../
    path_com_function= dir('*\common_functions');
    i=i+1;
end
addpath(fullfile(path_com_function(1).folder,'\'));

% fit files downloading : allows the user to set values of delta kt and
% nk_kt regarding the sought error (in %)

%For r_investigation <0.2
%addpath(['C:\Users\lalloz\Documents\these\onde alfven\etude theorique numerique\'...
%    'modélisation_b_induit\fichier_fit_approx_b\radius_elec_0p5mm']);

%For r_investigation > 0.2
%addpath(['C:\Users\lalloz\Documents\these\onde alfven\etude theorique numerique\'...
%    'modélisation_b_induit\fichier_fit_approx_b\radius_elec_0p5mm\r_inv_max_0p8\bis'])
%% fitting files loading
%load('fit_r_investigation__delta_kt__to__err_relative_pc.mat')
%load('fit_r_investigation__delta_kt__to__nb_kt')
%load('fit_r_investigation__err_relative_pc__to__delta_kt')

%% Paramètre du problème
syms mu_0 sigma_gal visco_gal rho_gal omega_0 B0 h I0 kt_max epsilon

%initialisation grandeur exp
B0= [];%champ magnetique uniforme en Tesla


%Paramètres physiques
mu_0= vpa(1.25*sym(1e-6));
sigma_gal= vpa(3.46*sym(1e6));
visco_gal= vpa(0.0024); %Viscosité dynamique du Galstn en Pa.s
eta= 1/(sigma_gal * mu_0); %diffusivité mhd
rho_gal= vpa(6440);% masse volumique Galinstan kg/m3
nu_gal= visco_gal/rho_gal;

% Paramètres géométriques/expérimentaux
if ~exist('B0') || isempty(B0)
    B0= vpa(input('Write the value of the magnetic field B0 (in T):\n'));
end
diam_elec= vpa(1e-3); % en mètre
rayon_elec= diam_elec/2;
dist_inter_elec= 15*10^-3; % en mètre
h= vpa(0.1); % Distance between the top and bottom plate en meter

% Grandeurs caractéristiques
L_carac=h; % Characteristic length ==> height of the vessel
tau_eta= L_carac^2/eta; % Diffusive penetration time of the resistive skin effect
tau_nu= L_carac^2/nu_gal;  % Diffusive penetration time of the viscous skin effect
I0=3.18; % Current injected at by the electrode
j0= I0/(L_carac)^2; % Form of jz = j0/(pi*sigma_r^2)*exp(-(r/sigma_r)^2)
b0= mu_0*I0/(L_carac); % form of b_theta= b0*(1-exp(-((r)/sigma_r).^2))./(r);
E0= eta/L_carac * b0; % Characteristic electric field
joule_time= rho_gal/(sigma_gal*B0^2); % Joule time
tau_alfven=h*sqrt(rho_gal*mu_0)/B0; % Alfven time, time of propagation to cross the height h
% 2D time
TwoD_time= joule_time*10.^2;
TwoD_time_norm= TwoD_time*0.5;

% Paramètres adimensionnées
sigma_r= rayon_elec/L_carac;% diam_adm de l'électrode
hadm=h/L_carac;

% Nombre adm
Pm_num= nu_gal/(eta); % magnetic Prandtl number
va_num= B0/sqrt(rho_gal*mu_0); % Alfven speed
S_num=va_num*L_carac/eta; % Lundquist number based of the height of the vessel
Ha_num= B0*L_carac*sqrt(sigma_gal/visco_gal); % Hartmann number based on the height of the vessel

Freq_Alfven= B0/sqrt(rho_gal*mu_0)/L_carac; % Alfven frequency ==> inverse of the Alfven time

%%%% Controle parameter: frequency %%%%
frequence_forcage= [];% ; %en Hz %0.01*tau^-1
% frequence
if ~exist('frequence_forcage') || isempty(frequence_forcage)
    frequence_forcage= input('Write the frequency (in Hz) to test : \n');
end
omega_inf= vpa(2*sym(pi)*frequence_forcage);%390);
epsilon= tau_eta*omega_inf; % Dimensionless number equals to R_eta ==> w*tau_eta

% New control parameter
Rnu_fixed= 0;
Rnu_varying_only= 0;
Reta= tau_eta*2*vpa(pi)*frequence_forcage;%
Rnu= tau_nu*2*vpa(pi)*frequence_forcage;
if Rnu_fixed == 1
    Pm_num= Reta/Rnu;
    S_num= Ha_num*sqrt(Reta/Rnu);
    epsilon= Reta;
end
if Rnu_varying_only == 1
    Rnu= 5*10^6;
    Pm_num= Reta/Rnu;
    S_num= Ha_num*sqrt(Reta/Rnu);
    epsilon= Reta;
end
    %depth_penetration_eta= 1/sqrt(2*vpa('pi')*frequence_forcage/eta);
%depth_penetration_eta_adm= depth_penetration_eta/h;
%% Champ à calculer
% Entrer ici les champs à calculer
calculate_U= [];
calculate_V= [];
calculate_B= [];
calculate_E= [];
calculate_J= [];
calculate_time_evolution= 0;
save_data= [];
show_fig= [];

if ~exist('calculate_U') || isempty(calculate_U)
    calculate_U=input('Do you want to calculate the velocity field ? Yes (1), No (0) \n');
end
if ~exist('calculate_V') || isempty(calculate_V)
    calculate_V=input('Do you want to calculate the potential gradient ? Yes (1), No (0) \n');
end
if ~exist('calculate_B') || isempty(calculate_B)
    calculate_B=input('Do you want to calculate the magnetic field ? Yes (1), No (0) \n');
end
if ~exist('calculate_E') || isempty(calculate_E)
    calculate_E=input('Do you want to calculate the electric field ? Yes (1), No (0) \n');
end
if ~exist('calculate_J') || isempty(calculate_J)
    calculate_J=input('Do you want to calculate the radial current density ? Yes (1), No (0) \n');
end
if ~exist('calculate_time_evolution') || isempty(calculate_time_evolution)
    calculate_time_evolution=input(...
        'Do you want to calculate the time evolution along z ? Yes (1), No (0) \n');
end
if ~exist('save_data') || isempty(save_data)
    save_data=input(...
        'Do you want to save data ? Yes (1), No (0) \n');
end
if ~exist('show_fig') || isempty(show_fig)
    show_fig=input(...
        'Do you want to display figures ? Yes (1), No (0) \n');
end



%% Calcul des solutions de la relation de dispersion full MHD 
%(avec solution de Hartmann et Alfven)

% entrée des paramètres d'intérêt : 
% r_min : distance minimale en mètre depuis le centre de l'électrode
% r_max : distance maximale en mètre depuis le centre de l'électrode
% nb_point_r= nombre de points de calcul
% nb_kt= nombre de mode propre transversaux à utiliser pour l'approximation
%       du forcage magnetique (nb_kt est adimensionné par L_carac= 0.1m)
% delta_kt= pas entre deux modes transversaux
r_min= [];%0.2;%
r_max= [];%0.2;
nb_point_r= [];%150;
h_max= hadm;
nb_point_z= [];
nb_kt= 1200;%1800;%3200;%1000;%2100;
%delta_kt= 2.25; %% %0.5;
R= 1.45;
%type_forcing= 'Gaussian_deriv';%'linear'

if ~exist('r_min') || isempty(r_min)
    r_min= input('Write the minimum radial distance scaled from the electrode, scaled with the height to test: \n');
end
if ~exist('r_max') || isempty(r_max)
    r_max= input('Write the maximal radial distance scaled from the electrode, scaled with the height to test: \n');
end
if ~exist('nb_point_r') || isempty(nb_point_r)
    nb_point_r= input('Write the number of points along r: \n');
end
if ~exist('nb_point_z') || isempty(nb_point_z)
    nb_point_z= input('Write the number of points along z: \n');
end
if ~exist('nb_kt') || isempty(nb_kt)
nb_kt= input('Write the number of transverse mode (k_t) to test: \n');
end
if ~exist('R') || isempty(R)
R= input('Write the radial size of the vessel (scaled with h) ?\n'); % exemple, 355   
end

% if ~exist('delta_kt') || isempty(delta_kt)
% delta_kt= input('Value of the discretisation for kt ?\n'); % exemple, 355   
% end

ordre_B= 1;
% location and transverse mode vectors
J1roots= vpa(besselzero(1,nb_kt,1)); %determine the roots of J1 in order to find the transverse mode
%R= (J1roots(2)-J1roots(1))/vpa(delta_kt); % taille du domaine de définition du forçage magnétique
kt_adm= J1roots/vpa(R); %mode transversaux
r_min_adm= vpa(r_min);
r_max_adm= vpa(r_max);
%r_investigation= vpa(logspace(log10(r_min_adm),log10(r_max_adm),nb_point_r));
r_investigation= vpa(linspace(r_min_adm,r_max_adm,nb_point_r));
fprintf('scaled step of location \\Delta_r = %4.3f \n', ...
    (r_max_adm-r_min_adm)/(nb_point_r-1));

%initialisation matrice
ratio_ampl_mat= vpa(zeros(1,1));
phase_up_mat= vpa(zeros(1,1));
phase_bot_mat= vpa(zeros(1,1));
dephasage_mat= vpa(zeros(1,1));


%% ampl/phase shift coefficient calculation

 disp("calculing MHD solution for a radial sweeping...")

% initialisation vecteur signal max top/bottom
E_env= [];
E_phase= [];
J_env= [];
J_phase= [];
gradV_env= [];
gradV_phase= [];
B_env= [];
B_phase= [];
A_env= [];
A_phase= [];
s1_struc= struct;
k1_struc= struct;
s2_struc= struct;
k2_struc= struct;
V_struc= struct;
E_struc= struct;
B_struc= struct;
A_struc= struct;
J_struc= struct;

% Calcul des coefficients Bb_ri, projections du forcage magnetique suivant
% les modes propres transversaux

%TF_b_theta= Hankel_T_type_forcing(type_forcing,kt_adm,R,sigma_r,precision);
TF_b_theta_dot_kt= ((1- 0.5*kt_adm*sigma_r*sqrt(vpa(pi)).*exp(-sigma_r^2*kt_adm.^2/8)...
           .*besseli(0.5,sigma_r^2*kt_adm.^2/8))); %adimensionné par L_carac*b0
Bb_ri=(kt_adm)./(vpa(pi)*(J1roots.*besselj(2,J1roots)).^2).*TF_b_theta_dot_kt; 


disp('Calculating MHD solutions...')
%DigitsOld= digits(precision);
[s1,k1,s2,k2,A,B,V,J,E]= alfven_non_homogene_vpa(Pm_num,S_num,hadm,...
    kt_adm(1),Bb_ri(1),epsilon,precision,'bottom'); % 0 == bottom
if length(kt_adm) > 1
parfor i =2:length(kt_adm)-1
    [s1i,k1i,s2i,k2i,Ai,Bi,Vi,Ji,Ei]= alfven_non_homogene_vpa(vpa(Pm_num),vpa(S_num),vpa(hadm),...
        vpa(kt_adm(i)),Bb_ri(i),epsilon,precision,'bottom');

    A= [A ; Ai];
    B= [B ; Bi];
    V= [V ; Vi];
    J= [J ; Ji];
    E= [E ; Ei];
    s1= [s1 ; s1i];
    k1= [k1 ; k1i];
    s2= [s2 ; s2i];
    k2= [k2 ; k2i];
    i
end
else
end
disp('end solution calculation')

%
% diminution précision des paramètres
    disp('diminution of digit number')
    s1_less= vpa(s1);%,precis_b);
    k1_less= vpa(k1);%,precis_b);
    s2_less= vpa(s2);%,precis_b);
    k2_less= vpa(k2);%,precis_b);
    hadm_less= vpa(hadm);%,precis_b);
    disp('end diminution of digit number')


parfor ind= 1:length(r_investigation)  
        
 %  ind   
    fprintf('progressing %4.2f%% ..\n',ind*100/(length(r_investigation)))
     
    %ri= vpa(r_investigation(ind));
    s1_struc(ind).value= s1_less;
    k1_struc(ind).value= k1_less;
    s2_struc(ind).value= s2_less;
    k2_struc(ind).value= k2_less;

    % Ajout des fonctions propres aux coefficient et diminution de la
    % précisison
    if calculate_U == 1 
        A_less= vpa(A.*besselj(ordre_B,kt_adm(1:end-1)'*r_investigation(ind)));%,precis_b);
        disp('Start calculation envelope & phase velocity field ')
        [A_env_j,A_phase_j]=amplitudes_somme_onde(A_less,s1_less,...
        k1_less,s2_less,k2_less,h_max,nb_point_z,precision,epsilon);       
        disp('End calculation envelope & phase velocity field ')
        A_env= [A_env ; A_env_j];
        A_phase= [A_phase ; A_phase_j];
        A_struc(ind).value= A_less;
    end
    if calculate_B == 1
        B_less= vpa(B.*besselj(ordre_B,kt_adm(1:end-1)'*r_investigation(ind)));%,precis_b);
        disp('Start calculation envelope & phase magnetic field ')
        [B_env_j,B_phase_j]=amplitudes_somme_onde(B_less,s1_less,...
        k1_less,s2_less,k2_less,h_max,nb_point_z,precision,epsilon);       
        disp('End calculation envelope & phase magnetique field ')
        B_env= [B_env ; B_env_j];
        B_phase= [B_phase ; B_phase_j];
        B_struc(ind).value= B_less;
    end
    if calculate_V == 1
        V_less= vpa(V.*besselj(ordre_B,kt_adm(1:end-1)'*r_investigation(ind)));%,precis_b);
        disp('Start calculation envelope & phase gradient phi')
        [gradV_env_j,gradV_phase_j]=amplitudes_somme_onde(V_less,s1_less,...
        k1_less,s2_less,k2_less,h_max,nb_point_z,precision,epsilon);    
        disp('End calculation envelope & phase gradient phi')
        gradV_env= [gradV_env ; gradV_env_j];
        gradV_phase= [gradV_phase ; gradV_phase_j];
        V_struc(ind).value= V_less;
    end
    if calculate_J == 1
        J_less= vpa(J.*besselj(ordre_B,kt_adm(1:end-1)'*r_investigation(ind)));%,precis_b);
        disp('Start calculation envelope & phase radial current ')
        [J_env_j,J_phase_j]=amplitudes_somme_onde(J_less,s1_less,...
        k1_less,s2_less,k2_less,h_max,nb_point_z,precision,epsilon);       
        disp('End calculation envelope & phase radial current ')
        J_env= [J_env ; J_env_j];
        J_phase= [J_phase ; J_phase_j];
        J_struc(ind).value= J_less;
    end
    if calculate_E ==1
        E_less= vpa(E.*besselj(ordre_B,kt_adm(1:end-1)'*r_investigation(ind)));%,precis_b);
        disp('Start calculation envelope & phase electric field ')
        [E_env_j,E_phase_j]=amplitudes_somme_onde(E_less,s1_less,...
        k1_less,s2_less,k2_less,h_max,nb_point_z,precision,epsilon);       
        disp('End calculation envelope & phase electric field ')
        E_env= [E_env ; E_env_j];
        E_phase= [E_phase ; E_phase_j];
        E_struc(ind).value= E_less;
    end
   
end

%% Saving section 
if save_data == 1
    disp('choose the folder to save data')
    selpath = uigetdir(pwd);
    save([selpath,'\donnees.mat'],'epsilon','Ha_num','S_num','B0',"frequence_forcage",'Rnu',...
        "Reta","Pm_num","r_investigation","r_min","r_max","R","nb_point_z","nb_point_r","h",...
        "h_max","nb_kt","E_env","E_phase", "J_env","J_phase","gradV_env","gradV_phase",...
        "B_env","B_phase","A_env","A_phase","s1_struc","k1_struc","s2_struc",...
        "k2_struc","V_struc","E_struc","B_struc","A_struc","J_struc","precision",...
        "calculate_U","calculate_V","calculate_B","calculate_E","calculate_J")

end

%% Parameter of amplitude and phase mapping (r,z)
exp_Ha= floor(log10(double(Ha_num)));
mantisse_Ha= double(Ha_num)/(10^(exp_Ha));
exp_Rnu= floor(log10(double(Rnu)));
mantisse_Rnu= double(Rnu)/(10^(exp_Rnu));
z_vect= double(0:h_max/(nb_point_z-1):h_max);
[R_mat , Z_mat]= meshgrid(double(r_investigation),z_vect);      
nb_line_env= 40;
nb_line_phase= 20;

phase_tick_vect= [-6*pi  -5*pi -4*pi -3*pi -5*pi/2 -2*pi -7*pi/4 -3*pi/2 -5*pi/4 -pi...
    -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi];
phase_tick_lbl_vect= {'$-6\pi$', '$-5\pi$','$-4\pi$', '$-3\pi$', '$-5\pi/2$',...
    '$-2\pi$', '$-7\pi/4$', '$-3\pi/2$','$-5\pi/4$','$-\pi$', '$-3\pi/4$',...
    '$-\pi/2$', '$-\pi/4$', '$0$', '$\pi/4$','$\pi/2$', '$3\pi/4$', '$\pi$'};
font_sz_lbl_lgd= 14;
font_sz_x= 14;
font_sz_y= 14;
font_sz_glb= 12;
%% Figure: enveloppe, phase, attenuation and phase shift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comments on the parameters displayed be the model.
% For all field (u_\theta,b_\theta,partial_r_phi,Er_jr) are show,:
% - the envelope written "|\hat{a_i(r,z}|", where \hat represent the fourier
% transform on the time variable and a_i a given component of the field
% - the phase written "\varphi[ a_i(r,z)]" = \arg(\hat{a_i(r,z})
% For Er and \partial_r_phi are also displayed :
% - the attenuation coefficient written "\alpha(r,z)"=
% \ln(|\hat{a_i(r,z}|/|\hat{a_i(r,z=0}|)
% - the phase shift written "\varphi_\Delta(r,z)"= \varphi[ a_i(r,z)] - 
% \varphi[ a_i(r,z=0)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if show_fig == 1

if calculate_B == 1

    % figure  for b
    % ampl mapping
    figure('Name','b amplitude vs r and z')
    contourf(R_mat,Z_mat,B_env',30);
    %hold on
    %plot(1/sqrt(Ha_num)*ones(nb_point_z,1),z_vect,'--k','linewidth',1.5)
    % hold on
    % plot(double(r_investigation),sqrt(2*pi*eta/(frequence_forcage))/h*ones(nb_point_r,1),':k','linewidth',1.5)
    colormap(WR2)%colormap_vect)
    c= colorbar;
    c.Label.String = '$|\hat{b}_{\theta}\left(r,z\right)|$';
    c.Label.Interpreter= 'latex';
    c.TickLabelInterpreter='latex';
    c.Label.FontSize= font_sz_lbl_lgd;
    f= gca;
    f.ColorScale= 'linear';
    caxis([min(abs(double(B_env(:,1:end-1))),[],'all'),...
        max(abs(double(B_env(:,1:end-1))),[],'all')])
    %c.Limits= ([min(double(B_env(:,1:end-1)),[],'all'),max(double(B_env),[],'all')]);
    f.TickLabelInterpreter= 'latex';
    xlabel('$r$')
    ylabel('$z$')
    f.FontSize= font_sz_glb;
    f.XLabel.FontSize= font_sz_x;
    f.YLabel.FontSize= font_sz_y;
    title(sprintf(['Mapping of the $b_\\theta$ amplitude for\n ' ...
        ' $Ha =%4.2f \\times 10^%i$, $R_\\eta = %4.0f$ and $R_\\nu =%4.2f \\times 10^%i$'],mantisse_Ha,exp_Ha,...
            double(Reta),mantisse_Rnu,exp_Rnu))
     
    % phase mapping
    figure('Name','b phase vs r and z');
    contourf(R_mat(1:end-1,:),Z_mat(1:end-1,:),transpose(B_phase(:,1:end-1)),20)%linecontour_B_phase);
    colormap(BWR2)
    c= colorbar;
    c.Label.String = '$\varphi\left[b_{\theta}\left(r,z\right)\right]$';
    caxis([-max(abs(double(B_phase(:,1:end-1))),[],'all'), max(abs(double(B_phase(:,1:end-1))),[],'all')]);
    c.Limits= ([min(double(B_phase(:,1:end-1)),[],'all'),max(double(B_phase(:,1:end-1)),[],'all')]);
    c.Label.Interpreter= 'latex';
    c.TickLabelInterpreter='latex';
    c.Label.FontSize= font_sz_lbl_lgd;
    c.Ticks= phase_tick_vect;
    c.TickLabels= phase_tick_lbl_vect;
    f= gca;
    f.TickLabelInterpreter= 'latex';
    f.Color
    xlabel('$r$')
    ylabel('$z$')
    f.FontSize= font_sz_glb;
    f.XLabel.FontSize= font_sz_x;
    f.YLabel.FontSize= font_sz_y;
    title(sprintf(['Mapping of the $b_\\theta$ phase for\n ' ...
        '$Ha =%4.2f \\times 10^%i$, $R_\\eta = %4.0f$ and $R_\\nu =%4.2f \\times 10^%i$'],mantisse_Ha,exp_Ha,...
            double(Reta),mantisse_Rnu,exp_Rnu))
    ylim([0, 1])

end

if calculate_U == 1

    ctr_line_vect= logspace(log10(min(A_env(:,2:end-1),[],'all')),...
        log10(max(A_env(:,2:end-1),[],'all')),30);
    % figure for u
    % ampl mapping
    figure('Name','u amplitude vs r and z')
    contourf(R_mat,Z_mat,A_env',ctr_line_vect);
    %hold on
    %plot(1/sqrt(Ha_num)*ones(nb_point_z,1),z_vect,'--k','linewidth',1.5)
    % hold on
    % plot(double(r_investigation),sqrt(2*pi*eta/(frequence_forcage))/h*ones(nb_point_r,1),':k','linewidth',1.5)
    colormap(WR2)%colormap_vect)
    c= colorbar;
    c.Label.String = '$|\hat{u}_{\theta}\left(r,z\right)|$';
    c.Label.Interpreter= 'latex';
    c.TickLabelInterpreter='latex';
    c.Label.FontSize= font_sz_lbl_lgd;
    f= gca;
    f.ColorScale= 'log';
    caxis([min(abs(double(A_env(:,2:end-1))),[],'all'),...
        max(abs(double(A_env(:,2:end-1))),[],'all')])
    %c.Limits= ([min(double(B_env(:,1:end-1)),[],'all'),max(double(B_env),[],'all')]);
    f.TickLabelInterpreter= 'latex';
    xlabel('$r$')
    ylabel('$z$')
    f.FontSize= font_sz_glb;
    f.XLabel.FontSize= font_sz_x;
    f.YLabel.FontSize= font_sz_y;
    title(sprintf(['Mapping of the $u_\\theta$ amplitude for\n ' ...
        ' $Ha =%4.2f \\times 10^%i$, $R_\\eta = %4.0f$ and $R_\\nu =%4.2f \\times 10^%i$'],mantisse_Ha,exp_Ha,...
            double(Reta),mantisse_Rnu,exp_Rnu))
     
    % phase mapping
    figure('Name','u phase vs r and z');
    contourf(R_mat(1:end-1,:),Z_mat(1:end-1,:),transpose(A_phase(:,1:end-1)),20)%linecontour_B_phase);
    colormap(BWR2)
    c= colorbar;
    c.Label.String = '$\varphi\left[u_{\theta}\left(r,z\right)\right]$';
    caxis([-max(abs(double(A_phase(:,1:end-1))),[],'all'), max(abs(double(A_phase(:,1:end-1))),[],'all')]);
    c.Limits= ([min(double(A_phase(:,1:end-1)),[],'all'),max(double(A_phase(:,1:end-1)),[],'all')]);
    c.Label.Interpreter= 'latex';
    c.TickLabelInterpreter='latex';
    c.Label.FontSize= font_sz_lbl_lgd;
    c.Ticks= phase_tick_vect;
    c.TickLabels= phase_tick_lbl_vect;
    f= gca;
    f.TickLabelInterpreter= 'latex';
    f.Color
    xlabel('$r$')
    ylabel('$z$')
    f.FontSize= font_sz_glb;
    f.XLabel.FontSize= font_sz_x;
    f.YLabel.FontSize= font_sz_y;
    title(sprintf(['Mapping of the $u_\\theta$ phase for\n ' ...
        '$Ha =%4.2f \\times 10^%i$, $R_\\eta = %4.0f$ and $R_\\nu =%4.2f \\times 10^%i$'],mantisse_Ha,exp_Ha,...
            double(Reta),mantisse_Rnu,exp_Rnu))
    ylim([0, 1])

end

if calculate_J == 1
    % figure  for b
    % ampl mapping
    figure('Name','u amplitude vs r and z')
    contourf(R_mat,Z_mat,J_env',30);
    %hold on
    %plot(1/sqrt(Ha_num)*ones(nb_point_z,1),z_vect,'--k','linewidth',1.5)
    % hold on
    % plot(double(r_investigation),sqrt(2*pi*eta/(frequence_forcage))/h*ones(nb_point_r,1),':k','linewidth',1.5)
    colormap(WR2)%colormap_vect)
    c= colorbar;
    c.Label.String = '$|\hat{j}_{r}\left(r,z\right)|$';
    c.Label.Interpreter= 'latex';
    c.TickLabelInterpreter='latex';
    c.Label.FontSize= font_sz_lbl_lgd;
    f= gca;
    f.ColorScale= 'log';
    caxis([min(abs(double(J_env(:,1:end-1))),[],'all'),...
        max(abs(double(J_env(:,1:end-1))),[],'all')])
    %c.Limits= ([min(double(B_env(:,1:end-1)),[],'all'),max(double(B_env),[],'all')]);
    f.TickLabelInterpreter= 'latex';
    xlabel('$r$')
    ylabel('$z$')
    f.FontSize= font_sz_glb;
    f.XLabel.FontSize= font_sz_x;
    f.YLabel.FontSize= font_sz_y;
    title(sprintf(['Mapping of the $j_r$ amplitude for\n ' ...
        ' $Ha =%4.2f \\times 10^%i$, $R_\\eta = %4.0f$ and $R_\\nu =%4.2f \\times 10^%i$'],mantisse_Ha,exp_Ha,...
            double(Reta),mantisse_Rnu,exp_Rnu))
     
    % phase mapping
    figure('Name','j phase vs r and z');
    contourf(R_mat(1:end-1,:),Z_mat(1:end-1,:),transpose(J_phase(:,1:end-1)),20)%linecontour_B_phase);
    colormap(BWR2)
    c= colorbar;
    c.Label.String = '$\varphi\left[j_{r}\left(r,z\right)\right]$';
    caxis([-max(abs(double(J_phase(:,1:end-1))),[],'all'), max(abs(double(J_phase(:,1:end-1))),[],'all')]);
    c.Limits= ([min(double(J_phase(:,1:end-1)),[],'all'),max(double(J_phase(:,1:end-1)),[],'all')]);
    c.Label.Interpreter= 'latex';
    c.TickLabelInterpreter='latex';
    c.Label.FontSize= font_sz_lbl_lgd;
    c.Ticks= phase_tick_vect;
    c.TickLabels= phase_tick_lbl_vect;
    f= gca;
    f.TickLabelInterpreter= 'latex';
    f.Color
    xlabel('$r$')
    ylabel('$z$')
    f.FontSize= font_sz_glb;
    f.XLabel.FontSize= font_sz_x;
    f.YLabel.FontSize= font_sz_y;
    title(sprintf(['Mapping of the $j_r$ phase for\n ' ...
        '$Ha =%4.2f \\times 10^%i$, $R_\\eta = %4.0f$ and $R_\\nu =%4.2f \\times 10^%i$'],mantisse_Ha,exp_Ha,...
            double(Reta),mantisse_Rnu,exp_Rnu))
    ylim([0, 1])

end
%
if calculate_V == 1
    att_coef_V= log(transpose(gradV_env./gradV_env(1:end,1)));
    phase_diff_V= (gradV_phase-gradV_phase(:,1))';
    ctr_line_vect= logspace(log10(min(gradV_env(:,:),[],'all')),...
        log10(max(gradV_env(:,:),[],'all')),30);
    % figure  for grad phi
    % ampl mapping
    figure('Name','grad phi amplitude vs r and z')
    contourf(R_mat,Z_mat,gradV_env',ctr_line_vect);
    %hold on
    %plot(1/sqrt(Ha_num)*ones(nb_point_z,1),z_vect,'--k','linewidth',1.5)
    % hold on
    % plot(double(r_investigation),sqrt(2*pi*eta/(frequence_forcage))/h*ones(nb_point_r,1),':k','linewidth',1.5)
    colormap(WR2)%colormap_vect)
    c= colorbar;
    c.Label.String = '$|\widehat{\partial_{r}\phi}\left(r,z\right)|$';
    c.Label.Interpreter= 'latex';
    c.TickLabelInterpreter='latex';
    c.Label.FontSize= font_sz_lbl_lgd;
    f= gca;
    f.ColorScale= 'log';
    caxis([min(abs(double(gradV_env(:,1:end-1))),[],'all'),...
        max(abs(double(gradV_env(:,1:end-1))),[],'all')])
    %c.Limits= ([min(double(B_env(:,1:end-1)),[],'all'),max(double(B_env),[],'all')]);
    f.TickLabelInterpreter= 'latex';
    xlabel('$r$')
    ylabel('$z$')
    f.FontSize= font_sz_glb;
    f.XLabel.FontSize= font_sz_x;
    f.YLabel.FontSize= font_sz_y;
    title(sprintf(['Mapping of the $\\partial_r \\phi$ amplitude for\n ' ...
        ' $Ha =%4.2f \\times 10^%i$, $R_\\eta = %4.0f$ and $R_\\nu =%4.2f \\times 10^%i$'],mantisse_Ha,exp_Ha,...
            double(Reta),mantisse_Rnu,exp_Rnu))
     
    % phase mapping
    figure('Name','grad phi phase vs r and z');
    contourf(R_mat(1:end-1,:),Z_mat(1:end-1,:),transpose(gradV_phase(:,1:end-1)),20)%linecontour_B_phase);
    colormap(BWR2)
    c= colorbar;
    c.Label.String = '$\varphi\left[\partial_{r}\phi\left(r,z\right)\right]$';
    caxis([-max(abs(double(gradV_phase(:,1:end-1))),[],'all'), max(abs(double(gradV_phase(:,1:end-1))),[],'all')]);
    c.Limits= ([min(double(gradV_phase(:,1:end-1)),[],'all'),max(double(gradV_phase(:,1:end-1)),[],'all')]);
    c.Label.Interpreter= 'latex';
    c.TickLabelInterpreter='latex';
    c.Label.FontSize= font_sz_lbl_lgd;
    c.Ticks= phase_tick_vect;
    c.TickLabels= phase_tick_lbl_vect;
    f= gca;
    f.TickLabelInterpreter= 'latex';
    f.Color
    xlabel('$r$')
    ylabel('$z$')
    f.FontSize= font_sz_glb;
    f.XLabel.FontSize= font_sz_x;
    f.YLabel.FontSize= font_sz_y;
    title(sprintf(['Mapping of the $\\partial_r\\phi$ phase for\n ' ...
        '$Ha =%4.2f \\times 10^%i$, $R_\\eta = %4.0f$ and $R_\\nu =%4.2f \\times 10^%i$'],mantisse_Ha,exp_Ha,...
            double(Reta),mantisse_Rnu,exp_Rnu))
    ylim([0, 1])


    % att coef mapping
    figure('Name','grad phi att coef vs r and z')
    contourf(R_mat,Z_mat,att_coef_V,30);
    %hold on
    %plot(1/sqrt(Ha_num)*ones(nb_point_z,1),z_vect,'--k','linewidth',1.5)
    % hold on
    % plot(double(r_investigation),sqrt(2*pi*eta/(frequence_forcage))/h*ones(nb_point_r,1),':k','linewidth',1.5)
    colormap(BWR2)%colormap_vect)
    c= colorbar;
    c.Label.String = '$|\alpha \left(r,z\right)|$';
    c.Label.Interpreter= 'latex';
    c.TickLabelInterpreter='latex';
    c.Label.FontSize= font_sz_lbl_lgd;
    f= gca;
    f.ColorScale= 'linear';
    caxis([-max(abs(double(att_coef_V(:,1:end-1))),[],'all'), max(abs(double(att_coef_V(:,1:end-1))),[],'all')])
    c.Limits= ([min(double(att_coef_V(:,1:end-1)),[],'all'),max(double(att_coef_V(:,1:end-1)),[],'all')]);
    f.TickLabelInterpreter= 'latex';
    xlabel('$r$')
    ylabel('$z$')
    f.FontSize= font_sz_glb;
    f.XLabel.FontSize= font_sz_x;
    f.YLabel.FontSize= font_sz_y;
    title(sprintf(['Mapping of the $\\alpha$ amplitude for\n ' ...
        ' $Ha =%4.2f \\times 10^%i$, $R_\\eta = %4.0f$ and $R_\\nu =%4.2f \\times 10^%i$'],mantisse_Ha,exp_Ha,...
            double(Reta),mantisse_Rnu,exp_Rnu))

    % phase shift mapping
    figure('Name','grad phi phase shift vs r and z');
    contourf(R_mat,Z_mat,phase_diff_V,20)%linecontour_B_phase);
    colormap(BWR2)
    c= colorbar;
    c.Label.String = '$\varphi_\Delta \left\{\partial_r\phi \left(r,z\right)\right\}$';
    caxis([-max(abs(double(phase_diff_V(:,1:end-1))),[],'all'), max(abs(double(phase_diff_V(:,1:end-1))),[],'all')]);
    c.Limits= ([min(double(phase_diff_V(:,1:end-1)),[],'all'),max(double(phase_diff_V(:,1:end-1)),[],'all')]);
    c.Label.Interpreter= 'latex';
    c.TickLabelInterpreter='latex';
    c.Label.FontSize= font_sz_lbl_lgd;
    c.Ticks= phase_tick_vect;
    c.TickLabels= phase_tick_lbl_vect;
    f= gca;
    f.TickLabelInterpreter= 'latex';
    f.Color
    xlabel('$r$')
    ylabel('$z$')
    f.FontSize= font_sz_glb;
    f.XLabel.FontSize= font_sz_x;
    f.YLabel.FontSize= font_sz_y;
    title(sprintf(['Mapping of the $\\partial_r\\phi$ phase difference for\n ' ...
        '$Ha =%4.2f \\times 10^%i$, $R_\\eta = %4.0f$ and $R_\\nu =%4.2f \\times 10^%i$'],mantisse_Ha,exp_Ha,...
            double(Reta),mantisse_Rnu,exp_Rnu))
    ylim([0, 1])

end


if calculate_E == 1
    att_coef_E= log(transpose(E_env./E_env(1:end,1)));
    phase_diff_E= (E_phase-E_phase(:,1))';
    ctr_line_vect= logspace(log10(min(E_env(:,:),[],'all')),...
        log10(max(E_env(:,:),[],'all')),30);
    % figure  for grad phi
    % ampl mapping
    figure('Name','E amplitude vs r and z')
    contourf(R_mat,Z_mat,E_env',ctr_line_vect);
    %hold on
    %plot(1/sqrt(Ha_num)*ones(nb_point_z,1),z_vect,'--k','linewidth',1.5)
    % hold on
    % plot(double(r_investigation),sqrt(2*pi*eta/(frequence_forcage))/h*ones(nb_point_r,1),':k','linewidth',1.5)
    colormap(WR2)%colormap_vect)
    c= colorbar;
    c.Label.String = '$|\hat{E}_{r} \left(r,z\right)|$';
    c.Label.Interpreter= 'latex';
    c.TickLabelInterpreter='latex';
    c.Label.FontSize= font_sz_lbl_lgd;
    f= gca;
    f.ColorScale= 'log';
    caxis([min(abs(double(E_env(:,1:end-1))),[],'all'),...
        max(abs(double(E_env(:,1:end-1))),[],'all')])
    %c.Limits= ([min(double(B_env(:,1:end-1)),[],'all'),max(double(B_env),[],'all')]);
    f.TickLabelInterpreter= 'latex';
    xlabel('$r$')
    ylabel('$z$')
    f.FontSize= font_sz_glb;
    f.XLabel.FontSize= font_sz_x;
    f.YLabel.FontSize= font_sz_y;
    title(sprintf(['Mapping of the $E_r$ amplitude for\n ' ...
        ' $Ha =%4.2f \\times 10^%i$, $R_\\eta = %4.0f$ and $R_\\nu =%4.2f \\times 10^%i$'],mantisse_Ha,exp_Ha,...
            double(Reta),mantisse_Rnu,exp_Rnu))
     
    % phase mapping
    figure('Name','E phase vs r and z');
    contourf(R_mat,Z_mat,transpose(E_phase),20)%linecontour_B_phase);
    colormap(BWR2)
    c= colorbar;
    c.Label.String = '$\varphi\left[E_{r}\left(r,z\right)\right]$';
    caxis([-max(abs(double(E_phase)),[],'all'), max(abs(double(E_phase)),[],'all')]);
    c.Limits= ([min(double(E_phase),[],'all'),max(double(E_phase),[],'all')]);
    c.Label.Interpreter= 'latex';
    c.TickLabelInterpreter='latex';
    c.Label.FontSize= font_sz_lbl_lgd;
    c.Ticks= phase_tick_vect;
    c.TickLabels= phase_tick_lbl_vect;
    f= gca;
    f.TickLabelInterpreter= 'latex';
    f.Color
    xlabel('$r$')
    ylabel('$z$')
    f.FontSize= font_sz_glb;
    f.XLabel.FontSize= font_sz_x;
    f.YLabel.FontSize= font_sz_y;
    title(sprintf(['Mapping of the $E_r$ phase for\n ' ...
        '$Ha =%4.2f \\times 10^%i$, $R_\\eta = %4.0f$ and $R_\\nu =%4.2f \\times 10^%i$'],mantisse_Ha,exp_Ha,...
            double(Reta),mantisse_Rnu,exp_Rnu))
    ylim([0, 1])


    % att coef mapping
    figure('Name','E att coef vs r and z')
    contourf(R_mat,Z_mat,att_coef_E,30);
    %hold on
    %plot(1/sqrt(Ha_num)*ones(nb_point_z,1),z_vect,'--k','linewidth',1.5)
    % hold on
    % plot(double(r_investigation),sqrt(2*pi*eta/(frequence_forcage))/h*ones(nb_point_r,1),':k','linewidth',1.5)
    colormap(BWR2)%colormap_vect)
    c= colorbar;
    c.Label.String = '$|\alpha \left(r,z\right)|$';
    c.Label.Interpreter= 'latex';
    c.TickLabelInterpreter='latex';
    c.Label.FontSize= font_sz_lbl_lgd;
    f= gca;
    f.ColorScale= 'linear';
    caxis([-max(abs(double(att_coef_E(:,1:end-1))),[],'all'), max(abs(double(att_coef_E(:,1:end-1))),[],'all')])
    c.Limits= ([min(double(att_coef_E(:,1:end-1)),[],'all'),max(double(att_coef_E(:,1:end-1)),[],'all')]);
    f.TickLabelInterpreter= 'latex';
    xlabel('$r$')
    ylabel('$z$')
    f.FontSize= font_sz_glb;
    f.XLabel.FontSize= font_sz_x;
    f.YLabel.FontSize= font_sz_y;
    title(sprintf(['Mapping of the $\\alpha$ amplitude for\n ' ...
        ' $Ha =%4.2f \\times 10^%i$, $R_\\eta = %4.0f$ and $R_\\nu =%4.2f \\times 10^%i$'],mantisse_Ha,exp_Ha,...
            double(Reta),mantisse_Rnu,exp_Rnu))

    % phase shift mapping
    figure('Name','E phase shift vs r and z');
    contourf(R_mat,Z_mat,phase_diff_E,20)%linecontour_B_phase);
    colormap(BWR2)
    c= colorbar;
    c.Label.String = '$\varphi_\Delta \left\{E_r \left(r,z\right)\right\}$';
    caxis([-max(abs(double(phase_diff_E(:,1:end-1))),[],'all'), max(abs(double(phase_diff_E(:,1:end-1))),[],'all')]);
    c.Limits= ([min(double(phase_diff_E(:,1:end-1)),[],'all'),max(double(phase_diff_E(:,1:end-1)),[],'all')]);
    c.Label.Interpreter= 'latex';
    c.TickLabelInterpreter='latex';
    c.Label.FontSize= font_sz_lbl_lgd;
    c.Ticks= phase_tick_vect;
    c.TickLabels= phase_tick_lbl_vect;
    f= gca;
    f.TickLabelInterpreter= 'latex';
    f.Color
    xlabel('$r$')
    ylabel('$z$')
    f.FontSize= font_sz_glb;
    f.XLabel.FontSize= font_sz_x;
    f.YLabel.FontSize= font_sz_y;
    title(sprintf(['Mapping of the $E_r$ phase shift for\n ' ...
        '$Ha =%4.2f \\times 10^%i$, $R_\\eta = %4.0f$ and $R_\\nu =%4.2f \\times 10^%i$'],mantisse_Ha,exp_Ha,...
            double(Reta),mantisse_Rnu,exp_Rnu))
    ylim([0, 1])

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evolution temporelle a coder !!                    %
% rq: peut être utilisé manuellement mais PAS EN RUN % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % calculation of the spatial time evolution
% 
% % calculation signal at a given time outside the Ha layer (10 time delta_Ha)
% Ha_num= S_num*sqrt(Pm_num^-1);                      
% z_Op= 10*Ha_num^-1;
% z_1m= 1-10*Ha_num^-1;
% t= [0:0.02*2*pi/epsilon:2*pi/(epsilon)];
% z_vect=0:h_max/(40-1):h_max;
% 
% 
% 
% gradV_time_struc= struct;   
% B_time_struc= struct;
% gradV_time_z0p= zeros(length(double(r_investigation)),length(double(t)));
% gradV_time_z1m= zeros(length(double(r_investigation)),length(double(t)));
% gradV_time_zi= zeros(length(double(r_investigation)),length(double(t)));
% B_time_z0p= zeros(length(double(r_investigation)),length(double(t)));
% B_time_z1m= zeros(length(double(r_investigation)),length(double(t)));
% B_time_zi= zeros(length(double(r_investigation)),length(double(t)));
% %E_time_z0p1= zeros(length(r_investigation),length(t));
% %E_time_z0p9= zeros(length(r_investigation),length(t));
% 
% ind_zi= 1;
% 
% s1=  s1_struc(1).value;           
% k1= k1_struc(1).value;
% s2=  s2_struc(1).value;
% k2= k2_struc(1).value;
% 
% 
% for zi= z_vect
%     parfor ind_r= 1:length(r_investigation)
%     
%     
%         V=  V_struc(1,ind_r).value;
%         B=  B_struc(1,ind_r).value;
%       %    E= E_struc(ind).value;
%         phase_num= vpa('0');
%         if ind_zi == 1
%             ind_r
% %             gradV_time_z0p(ind,:)= calcul_time_signal_given_z(V,s1,k1,s2,k2,epsilon,...
% %                 z_Op,t,precision,phase_num);
% %             gradV_time_z1m(ind,:)= calcul_time_signal_given_z(V,s1,k1,s2,k2,epsilon,...
% %                 z_1m,t,precision,phase_num);
% %             B_time_z0p(ind,:)= calcul_time_signal_given_z(B,s1,k1,s2,k2,epsilon,...
% %                 z_Op,t,precision,phase_num);
% %             B_time_z1m(ind,:)= calcul_time_signal_given_z(B,s1,k1,s2,k2,epsilon,...
% %                 z_1m,t,precision,phase_num);
%             
%         end
%         ind_r
%         gradV_time_zi(ind_r,:)= calcul_time_signal_given_z(V,s1,k1,s2,k2,epsilon,...
%             zi,t,[],phase_num);
%         B_time_zi(ind_r,:)= calcul_time_signal_given_z(B,s1,k1,s2,k2,epsilon,...
%             zi,t,[],phase_num);
%         ind_r
%     end
%     gradV_time_struc(ind_zi).z_value= gradV_time_zi;
%     B_time_struc(ind_zi).z_value= B_time_zi;
%     fprintf('progressing %4.2f%% ..\n',100*ind_zi/nb_point_z)
%     ind_zi= ind_zi +1;
% end
%     
% %% save time data
% selpath = uigetdir(['C:\Users\lalloz\Documents\these\onde alfven\etude theorique numerique'...
%     '\1_electrode_model_Results\one plate forcing case\location_sweeping']);
% save([selpath,'\time_data_evolution.mat'],'z_Op','z_1m','Ha_num',...
%     't','gradV_time_z0p','gradV_time_z1m','B_time_z0p','B_time_z1m','r_investigation',...
%     'gradV_time_struc','B_time_struc','z_vect')
%   
%    
% %% save
% selpath = uigetdir(['C:\Users\lalloz\Documents\these\onde alfven\etude theorique numerique'...
%     '\1_electrode_model_Results\one plate forcing case']);
% %%
% save(['donnees_nb_z_50.mat'],'Bb_ri','gradV_env',... %_nb_z_2 %selpath,'\
%     'gradV_phase','E_env','E_env','B_env','B_phase',...
%     'A_env','A_phase','Freq_Alfven','frequence_forcage','Ha_num',...
%     'R','S_num','B0','r_investigation','r_min','r_max','nb_point_r',...
%     'sigma_r','epsilon','phase_bot_mat','hadm','dist_inter_elec',...
%    'J1roots','dist_inter_elec','kt_adm','delta_kt','nb_kt','s1_struc',...
%     'k1_struc','s2_struc','k2_struc','V_struc','E_struc','B_struc','A_struc',...
%     'Rnu','Reta','nb_point_z','precision')
% 
% 
% 
% 
% %% figure showing the spatial time evolution of the B and gradV
% 
% figure('Name','B_spatial_time_evolution')
% B_max= max(double(B_time_struc(1).z_value),[],'all');
% B_min= min(double(B_time_struc(1).z_value),[],'all');
% linecontour_B_time= linspace(B_min/2,B_max/2,60);
% for i_time= 1:length(t)
%     B_time_i= [];
%     for j= 1:length(z_vect)
%         B_time_i= [B_time_i, B_time_struc(j).z_value(:,i_time)];
%     end
%     contourf(R_mat,Z_mat,B_time_i',linecontour_B_time)
%     colormap(BWR2)
%     c= colorbar;
%     caxis([B_min, B_max]);
%     c.Limits= ([B_min,B_max]);
%     pause(0.2)
% end
% %%
% fig_gradV_sp_time_evolution= figure('Name','gradV_spatial_time_evolution');
% gradV_max= 500;%max(double(gradV_time_struc(1).z_value),[],'all');
% gradV_min= -500;%min(double(gradV_time_struc(1).z_value),[],'all');
% linecontour_gradV_time= linspace(-500,500,100);
% for i_time= 1:length(t)
%     gradV_time_i= [];
%     for j= 1:length(z_vect)
%         gradV_time_i= [gradV_time_i, gradV_time_struc(j).z_value(:,i_time)];
%     end
%     gradV_time_i(gradV_time_i>=gradV_max)= gradV_max;
%     gradV_time_i(gradV_time_i<=gradV_min)= gradV_min;
%     contourf(R_mat,Z_mat,gradV_time_i',linecontour_gradV_time)
%     colormap(BWR2);
%     c= colorbar;
%     caxis([-500, 500]);
%     %c.Limits= ([-500,500]);
%     c.Label.String = '$\nabla\phi_r\left(r,z\right)$';
%     c.Label.Interpreter= 'latex';
%     c.Label.FontSize= 14;
%     c.TickLabelInterpreter='latex';
%      f= gca;
%      f.ColorScale= 'linear';
%     f.TickLabelInterpreter= 'latex';
%     xlabel('$r$')
%     ylabel('$z$')
%     f.FontSize= 12;
%     f.XLabel.FontSize= 14;
%     f.YLabel.FontSize= 14;
%     pause(0.2)
%    video_struc(i_time)= getframe(fig_gradV_sp_time_evolution);
%     hold off
% end
% 
% 
% 
% %% sauvegarde video
% 
%     for i=1:length(video_struc)%/2-1
%         video_struc_bis(i)= video_struc(i);
%         i
%     end
%     v= VideoWriter('gradV_along_z_and_r_time_evolution.avi','Uncompressed AVI');
%     v.FrameRate= 5;
%     %v.Quality= 100;
%     %v.CompressionRatio= 1;
%     open(v)
%     writeVideo(v,video_struc_bis)
%     close(v)
    clear v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fin signaux temporels 
%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autres mapping non temporels pouvant être tracés %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% B_test= B_env;
% B_test(B_env == 0)=NaN;
% 
% transit_time= (gradV_phase-gradV_phase(:,1))'/double(epsilon);
% phase_diff= (gradV_phase-gradV_phase(:,1))';
% %A_test= A_env;
% %A_test(A_env == 0)=NaN;
% linecontour_B_env= logspace(-4,... %log10(min(double(B_test),[],'all',"omitnan"))
%     log10(max(double(B_test),[],'all')),nb_line_env);
% linecontour_gradV_env= logspace(log10(1e-10),...
%     log10(5e5),nb_line_env);
% linecontour_gradV_env_norm = logspace(log10(1e-11),...
%     log10(1),nb_line_env);
% linecontour_B_phase= linspace(min(double(B_phase),[],'all',"omitnan"),...
%     max(double(B_phase),[],'all'),nb_line_phase);
% linecontour_gradV_phase= linspace(min(double(gradV_phase),[],'all',"omitnan"),...
%     max(double(gradV_phase),[],'all'),nb_line_phase);
% 
% linecontour_attenuation_coef_env= linspace(-9,1,20);%real([linspace(-12,...
% %-0.1,15), linspace(0.1,1,10)]);%4
% 
% %linecontour_A_env= logspace(-2,...
% %    2,nb_line_env);
% 
% % limits for the colorbar
% b_ampl_inf= 100*min(abs(double(B_env(:,1:end-1))),[],'all');
% b_ampl_sup= max(abs(double(B_env(:,1:end-1))),[],'all');
% gradV_ampl_inf= 1;%min(abs(double(gradV_env(:,1:end-1))),[],'all');
% gradV_ampl_sup= 10^4;%max(abs(double(gradV_env(:,1:end-1))),[],'all');
% %A_ampl_inf= min(abs(double(A_env(:,1:end-1))),[],'all');
% %A_ampl_sup= max(abs(double(A_env(:,1:end-1))),[],'all');
% 
% depth_penetration_eta= 1/sqrt(vpa(pi)*frequence_forcage/eta);
% depth_penetration_eta_adm= depth_penetration_eta/h;
% 
% 
% %%% max gradient phi, B0= 10T
% % max_f_56Hz= 3.70e5
% % max_f_56Hz= 6.44e4
% % max_f_560Hz= 2.17e4
% % max_f_1114Hz= 1.55e4
% 
% %%% min gradient phi, B0= 10T
% % min_f_56Hz= 3.70e5
% % min_f_56Hz= 6.44e4
% % min_f_560Hz= 0.34
% % min_f_1114Hz= 1.55e4
% 
%        
% %%  figure gradV
% figure('Name','gradV amplitude vs r and z')
% s=contourf(R_mat(:,1:end),Z_mat(:,1:end),transpose(gradV_env(1:end,:)),linecontour_gradV_env);
% hold on
% plot(1/sqrt(Ha_num)*ones(nb_point_z,1),z_vect,'--k','linewidth',1.5)
% hold on
% plot(double(r_investigation),depth_penetration_eta_adm*ones(nb_point_r,1),':w','linewidth',1.5)
% colormap(WR2(40))
% c= colorbar;
% caxis([min(gradV_env,[],"all"),max(gradV_env,[],"all")]);
% %c.Limits= ([min(double(B_phase(:,1:end-1)),[],'all'),max(double(B_phase(:,1:end-1)),[],'all')]);
% c.Label.String = '$|\nabla\phi_r\left(r,z\right)|$';
% c.Label.Interpreter= 'latex';
% c.Label.FontSize= 14;
% c.TickLabelInterpreter='latex';
% f= gca;
% f.ColorScale= 'log';
% f.TickLabelInterpreter= 'latex';
% xlabel('$r$')
% ylabel('$z$')
% f.FontSize= 12;
% f.XLabel.FontSize= 14;
% f.YLabel.FontSize= 14;
% title(sprintf(['Mapping of the radial potential gradient amplitude for\n $Ha= %4.2e,\\, R_{\\eta}/2\\pi = %4.2f,',...
%            ' \\,R_\\nu/2\\pi = %4.2e$'],Ha_num,epsilon/(2*vpa(pi)),Rnu/(2*vpa(pi))))
%  
% 
% % figure gradV normalised max
% figure('Name','gradV amplitude vs r and z')
% s=contourf(R_mat(:,1:end),Z_mat(:,1:end),transpose(gradV_env(1:end,:)/max(gradV_env,[],"all")),linecontour_gradV_env_norm);
% hold on
% plot(1/sqrt(Ha_num)*ones(nb_point_z,1),z_vect,'--k','linewidth',1.5)
% hold on
% plot(double(r_investigation),depth_penetration_eta_adm*ones(nb_point_r,1),':w','linewidth',1.5)
% colormap(WR2(40))
% c= colorbar;
% caxis([1e-4,1]);
% %c.Limits= ([min(double(B_phase(:,1:end-1)),[],'all'),max(double(B_phase(:,1:end-1)),[],'all')]);
% c.Label.String = '$|\nabla\phi_r\left(r,z\right)|/\max\left(|\nabla\phi_r|\right)$';
% c.Label.Interpreter= 'latex';
% c.Label.FontSize= 14;
% c.TickLabelInterpreter='latex';
% f= gca;
% f.ColorScale= 'log';
% f.TickLabelInterpreter= 'latex';
% xlabel('$r$')
% ylabel('$z$')
% f.FontSize= 12;
% f.XLabel.FontSize= 14;
% f.YLabel.FontSize= 14;
% title(sprintf(['Mapping of the radial potential gradient amplitude for\n $Ha= %4.2e,\\, R_{\\eta}/2\\pi = %4.2f,',...
%            ' \\,R_\\nu/2\\pi = %4.2e$'],Ha_num,epsilon/(2*vpa(pi)),Rnu/(2*vpa(pi))))
%  
%        
% %     
% figure('Name','gradV phase vs r and z')
% contourf(R_mat,Z_mat,gradV_phase',linecontour_gradV_phase);
% colormap(BWR2(50))
% c= colorbar;
% c.Label.String = '$\varphi\left[b_{\theta}\left(r,z\right)\right]$';
% caxis([-max(abs(double(gradV_phase)),[],'all'), max(abs(double(gradV_phase)),[],'all')]);
% c.Limits= ([min(double(gradV_phase),[],'all'),max(double(gradV_phase),[],'all')]);
% c.Label.String = '$\varphi\left[\nabla\phi_r\left(r,z\right)\right]$';
% c.Label.Interpreter= 'latex';
% c.Label.FontSize= 14;
% c.Ticks= [-2*pi -7*pi/4 -3*pi/2 -5*pi/4 -pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi];
% c.TickLabels= {'$-2\pi$', '$-7\pi/4$', '$-3\pi/2$', '$-5\pi/4$','$-\pi$', '$-3\pi/4$',...
%     '$-\pi/2$', '$-\pi/4$', '$0$', '$\pi/4$', '$\pi/2$', '$3\pi/4$', '$\pi$'};
% c.TickLabelInterpreter='latex';
% f= gca;
% f.TickLabelInterpreter= 'latex';
% xlabel('$r$')
% ylabel('$z$')
% f.FontSize= 12;
% f.XLabel.FontSize= 14;
% f.YLabel.FontSize= 14;
% title(sprintf(['Mapping of the radial potential gradient phase for\n $Ha= %4.2e,\\, R_{\\eta}/2\\pi = %4.2f,',...
%            ' \\,R_\\nu/2\\pi = %4.2e$'],Ha_num,epsilon/(2*vpa(pi)),Rnu/(2*vpa(pi))))
%  
% 
% % transit time 
% figure('Name','gradV transit time vs r and z')
% contourf(R_mat,Z_mat,transit_time,15);%linecontour_gradV_phase);
% hold on
% %contour(R_mat,Z_mat,transit_time,[0 0],'k','linewidth',2);
% colormap(BWR2(100))
% c= colorbar;
% 
% caxis([-max(abs(transit_time),[],'all'), max(abs(transit_time),[],'all')]);
% c.Limits= ([min(double(transit_time),[],'all'),max(double(transit_time),[],'all')]);
% %c.Limits= ([min(double(transit_time),[],'all'),0.01])
% c.Label.String = '$\tau\left(r,z\right)$';
% c.Label.Interpreter= 'latex';
% c.Label.FontSize= 14;
% %c.Ticks= [-2*pi -7*pi/4 -3*pi/2 -5*pi/4 -pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi];
% %c.TickLabels= {'$-2\pi$', '$-7\pi/4$', '$-3\pi/2$', '$-5\pi/4$','$-\pi$', '$-3\pi/4$',...
% %    '$-\pi/2$', '$-\pi/4$', '$0$', '$\pi/4$', '$\pi/2$', '$3\pi/4$', '$\pi$'};
% c.TickLabelInterpreter='latex';
% f= gca;
% f.TickLabelInterpreter= 'latex';
% xlabel('$r$')
% ylabel('$z$')
% f.FontSize= 12;
% f.XLabel.FontSize= 14;
% f.YLabel.FontSize= 14;
% title(sprintf(['Mapping of the transit time for \n $Ha= %4.2e,\\, R_{\\eta}/2\\pi = %4.2f,',...
%            ' \\,R_\\nu/2\\pi = %4.2e$'],Ha_num,epsilon/(2*vpa(pi)),Rnu/(2*vpa(pi))))
%  
% % phase diff for gradV
% 
% figure('Name','gradV phase diff vs r and z')
% contourf(R_mat,Z_mat,phase_diff,15);%linecontour_gradV_phase);
% hold on
% contour(R_mat,Z_mat,transit_time,[0 0],'k','linewidth',2);
% colormap(BWR2(60))
% c= colorbar;
% 
% caxis([-max(abs(phase_diff),[],'all'), max(abs(phase_diff),[],'all')]);
% c.Limits= ([min(double(phase_diff),[],'all'),max(double(phase_diff),[],'all')]);
% %c.Limits= ([min(double(transit_time),[],'all'),0.01])
% c.Label.String= ['$\varphi\left[\nabla\phi_r\left(r,z\right)\right]',...
%     ' - \varphi\left[\nabla\phi_r\left(r,z=0\right)\right]$'];
% c.Label.Interpreter= 'latex';
% c.Label.FontSize= 14;
% c.Ticks= [-6*pi  -5*pi -4*pi -3*pi -5*pi/2 -2*pi -7*pi/4 -3*pi/2 -5*pi/4 -pi...
%     -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi];
% c.TickLabels= {'$-6\pi$', '$-5\pi$','$-4\pi$', '$-3\pi$', '$-5\pi/2$','$-2\pi$', '$-7\pi/4$', '$-3\pi/2$',...
%     '$-5\pi/4$','$-\pi$', '$-3\pi/4$','$-\pi/2$', '$-\pi/4$', '$0$', '$\pi/4$',...
%     '$\pi/2$', '$3\pi/4$', '$\pi$'};
% c.TickLabelInterpreter='latex';
% f= gca;
% f.TickLabelInterpreter= 'latex';
% xlabel('$r$')
% ylabel('$z$')
% f.FontSize= 12;
% f.XLabel.FontSize= 14;
% f.YLabel.FontSize= 14;
% title(sprintf(['Mapping of the phase difference for \n $Ha= %4.2e,\\, R_{\\eta}/2\\pi = %4.2f,',...
%            ' \\,R_\\nu/2\\pi = %4.2e$'],Ha_num,epsilon/(2*vpa(pi)),Rnu/(2*vpa(pi))))
%  
%            
% % attenauation coeff
% nb_pt_r_bis= 20;
% figure('Name','attenuation coef vs r and z');
% s=contourf(R_mat(:,1:end),Z_mat(:,1:end),log(transpose(gradV_env./gradV_env(1:end,1))),linecontour_attenuation_coef_env);%linecontour_attenuation_coef_env);
% hold on
% plot(1/sqrt(Ha_num)*ones(nb_point_z,1),z_vect,'--w','linewidth',1.5)
% hold on
% plot(linspace(r_min,r_max,nb_pt_r_bis),...
%     double(depth_penetration_eta_adm*ones(nb_pt_r_bis,1)),'.k','linewidth',1.5)
% hold on
% contour(R_mat,Z_mat,log(transpose(gradV_env./gradV_env(1:end,1))),[0 0],'k','linewidth',2);
% colormap(BWR2(100))
% c= colorbar;
% caxis([-14,14]);
% c.Limits=([min(double(log(gradV_env./gradV_env(1:end,1))),[],'all'),...
%     max(double(log(gradV_env./gradV_env(1:end,1))),[],'all')]); %[linecontour_attenuation_coef_env(1),linecontour_attenuation_coef_env(end)] ;%
% c.Label.String = '$\alpha \left( r, z\right)$';
% c.Label.Interpreter= 'latex';
% c.Label.FontSize= 14;
% c.TickLabelInterpreter='latex';
% ax= gca;
% ax.ColorScale= 'linear';
% ax.TickLabelInterpreter= 'latex';
% ax.PlotBoxAspectRatio= [0.6000, 1.0000, 1.0000];
% ax.Children
% xlabel('$r$')
% ylabel('$z$')
% ax.FontSize= 12;
% ax.XLabel.FontSize= 14;
% ax.YLabel.FontSize= 14;
% title(sprintf(['Mapping of the radial potential gradient amplitude for \n $Ha= %4.2e,\\, R_{\\eta}/2\\pi = %4.2f,',...
%            ' \\,R_\\nu/2\\pi = %4.2e$'],Ha_num,epsilon/(2*vpa(pi)),Rnu/(2*vpa(pi))))
%        
%        
% % Amplitude ratio
% 
% i_lim= length(r_investigation);
% gradV_max_top= gradV_env(1:i_lim,end); 
% gradV_max_bot= gradV_env(1:i_lim,1);
% E_max_top= E_env(1:i_lim,end);  
% E_max_bot= E_env(1:i_lim,1); 
% 
% gradV_phase_top= unwrap(gradV_phase(1:i_lim,end)); %unwrap
% gradV_phase_bot= gradV_phase(1:i_lim,1); 
% E_phase_top= unwrap(E_phase(1:i_lim,end));  %unwrap
% E_phase_bot= (E_phase(1:i_lim,1)); %unwrap
% 
% gradV_ampl_ratio= double(gradV_max_top./gradV_max_bot); 
% E_ampl_ratio= double(E_max_top./E_max_bot);
% 
% E_dephasage= E_phase_top - E_phase_bot;
% gradV_dephasage= gradV_phase_top - gradV_phase_bot;
% 
% 
% %gradV_ampl_ratio= double(gradV_max_top./gradV_max_bot);
% %E_ampl_ratio= double(E_max_bot./E_max_top);
% 
% % transit time and velocity
% joule_time= rho_gal/(sigma_gal*B0^2);
% gradV_time_delay= -double(gradV_dephasage./(2*vpa('pi')*double(frequence_forcage)));
% E_time_delay= -double(E_dephasage./(2*vpa('pi')*double(frequence_forcage)));
% gradV_velocity= double(h./gradV_time_delay);
% E_velocity= double(h./E_time_delay);
% 
% gradV_velocity_adm= double(gradV_velocity/va_num);
% E_velocity_adm= double(E_velocity/va_num);
% 
% % 2D time
% TwoD_time= joule_time*1.^2;
% TwoD_time_norm= TwoD_time*0.5;
% 
% 
% % Figure spatial attenuation coefficient
% figure
% %yyaxis left
% plot(r_investigation,log(gradV_ampl_ratio))
% hold on
% xlabel('$r$')
% ylabel('$s_{\nabla\phi}$')
% title(sprintf(['Spatial attenuation coefficient for the potential gradient\n'...
% ' ($S= %4.1f, \\Omega= %4.2f$)'],double(S_num),2*pi*double(frequence_forcage/Freq_Alfven)))
% %txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
% %t= text(700,0.4,txt);
% %t.Units= 'normalized';
% f_att_coef_ax= gca;
% f_att_coef_ax.TickLabelInterpreter= 'latex';
% f_att_coef_ax.FontSize= 10;
% f_att_coef_ax.YLabel.FontSize = 15;
% f_att_coef_ax.XLabel.FontSize = 12;
% grid on
% 
% % Figure transit time scaled with the alfven time coefficient
% figure
% %yyaxis left
% plot(r_investigation,gradV_time_delay*Freq_Alfven)
% hold on
% xlabel('$r$')
% ylabel('$\tau_{\nabla\phi}/\tau_A$')
% title(sprintf(['scaled transit time for the potential gradient\n'...
% ' ($S= %4.1f, \\Omega= %4.2f$)'],double(S_num),2*pi*double(frequence_forcage/Freq_Alfven)))
% %txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
% %t= text(700,0.4,txt);
% %t.Units= 'normalized';
% f_transit_time_ax= gca;
% f_transit_time_ax.TickLabelInterpreter= 'latex';
% f_transit_time_ax.FontSize= 10;
% f_transit_time_ax.YLabel.FontSize = 15;
% f_transit_time_ax.XLabel.FontSize = 12;
% grid on

%%

% %% figure Evolution along r at a specific time
% time_loc= 3;
% 
% fig_animation= figure;
% 
% for i= 1:length(t)
% %yyaxis left
% plot(r_investigation,gradV_time_z0p(:,i))
% hold on
% plot(r_investigation,gradV_time_z1m(:,i))
% hold on
% xlabel('$r$')
% ylabel('$\nabla\phi_r (r,z,t)$') %'$b_\theta (r,z,t) $'
% title(sprintf(['Magnetic field regarding the location r at a given time\n'...
% ' ($S= %4.1f ,\\, \\left( 2 \\pi \\right)^{-1}\\Omega= %4.2f$)'],double(S_num),double(frequence_forcage/Freq_Alfven)))
% %txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
% %t= text(700,0.4,txt);
% %t.Units= 'normalized';
% legend('$z=0^+$','$z=1^-$','Interpreter','latex') %
% f_transit_time_ax= gca;
% %f_transit_time_ax.YScale= 'log';
% f_transit_time_ax.TickLabelInterpreter= 'latex';
% f_transit_time_ax.FontSize= 10;
% f_transit_time_ax.YLabel.FontSize = 15;
% f_transit_time_ax.XLabel.FontSize = 12;
% xlim([double(r_investigation(1)) double(r_investigation(end))])
% ylim([-max(double(gradV_time_z0p),[],'all'),max(double(gradV_time_z0p),[],'all')])
% grid on
% pause(0.1)
% video_struc(i)= getframe(fig_animation);
% hold off
% end
% 
% 
% 
% %% sauvegarde video
% 
%     for i=1:length(video_struc)%/2-1
%         video_struc_bis(i)= video_struc(i);
%         i
%     end
%     v= VideoWriter('b_z0p_and_z1m_vs_r.avi');
%     v.FrameRate= 10;
%     open(v)
%     writeVideo(v,video_struc_bis)
%     close(v)
%     clear v
% 
% 
% 
% 
% %% grad phi
% size_ft= 12;
% figure
% subplot(2,1,1)
%  plot(r_investigation,gradV_max_top)
%     xlabel('$\tilde{r}/h$')
%     ylabel('$|\nabla\phi|/\nabla\phi_0$')
%     title('top plate')
%     
% subplot(2,1,2)
%   plot(r_investigation,gradV_max_bot)
%     xlabel('$\tilde{r}/h$')
%     ylabel('$|\nabla\phi|/\nabla\phi_0$')
%     title('bottom plate')
%     
% sgtitle(sprintf(['Evolution of the potential gradient amplitude \n regarding '...
%     'the distance from the electrode \n(S= %4.1f, Frequency= %4.1f Hz, $\\mathrm{E}_0$= %4.1e V/m)'],...
%     double(S_num),frequence_forcage,double(E0))) 
% 
% 
% %% amplitude ratio
% 
% figure
% 
% %subplot(2,1,1)  %yyaxis left
% plot(r_investigation,gradV_ampl_ratio)
% hold on 
% %p= plot(r_investigation,mean_gradV_ratio,'b');
% %p.LineWidth = 1;
% xlabel('r/$\sigma_r$')
% ylabel(sprintf('amplitude ratio between \n the top and bottom plate'))
% 
% legend('radial \nabla{\phi} ratio')
% title(sprintf(['Potentiel gradient amplitude ratio in function of the distance'...
%     ' from\n the electrode (S= %4.1f, Frequency = %4.1f Hz)'],double(S_num),frequence_forcage)) 
% 
% %yyaxis right
% %subplot(2,1,2) 
% % %hold on
% % plot(r_investigation,gradV_max_bot)
% % hold on
% % plot(r_investigation,gradV_max_top)
% % ylabel(sprintf('amplitude of the \n non scaled potential gradient '))
% % xlabel('non scaled distance from the centre of the electrode')
% % 
% % legend('Amplitude of the bottom gradV','Amplitude of the top gradV')
% % title(sprintf(['Potentiel gradient amplitude in function of the distance'...
% %     ' from\n the electrode (S= %4.1f)'],double(S_num))) 
% % 
% % txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0)
% % t= text(0.14,400,txt);
% % t.Units= 'normalized';
% % set(gca, 'YScale', 'log')
% 
% % figure
% % plot(epsilon_mat,gradV_ampl_ratio)
% % figure
% % plot(epsilon_mat,gradV_ampl_ratio)
% % xlabel('relative error of $b_{approx}$ en percentage')
% % ylabel('amplitude ratio between the top and bottom plate')
% % legend('potential gradient ratio')



