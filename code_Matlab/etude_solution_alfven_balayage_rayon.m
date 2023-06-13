clear all
clc
%close all


% Calcul des amplitudes pour le potentiel, vitesse, champ magnétique
% induits au niveau des plaques inférieures et supérieures. Le programme 
% utilise la fonction alfven_non_homogene qui permet de résoudre la 
% relation de dispersion full
% MHD pour un cas de forçage non_homogène 
% (existance d'un nombre d'onde transverse k_t)
% On considère ici que seul le courant injecté produit de la vorticité (la
% plaque inférieure n'est pas oscillante)

% Le forcage peut être produit par les deux plaques a deux frequences données
% différentes. Néanmoins la distribution spatiale du forcage électrique est
% le même sur la plaque supérieure que sur la plaque inférieure

% Ce programme fournit la valeur du ratio d'amplitude pour un champ donné
% lors q'un balayage de la distance du point de mesure depuis le centre de
% l'electrode. 
%L_carac est la hauteur du récipient

% Adimensionnement de :
% - j par I0/L_carac^2
% - b par b0= mu_0*I0/(L_carac)
% - TF_b_theta_kt par L_carac^2*b0
% - u par u0= eta/L_carac* b0/B0
% - tau= L_carac^2/eta;
% - E_0= eta/L_carac * b_0
%%%%
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
    path_com_function= dir('common_functions');
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
%Paramètres physiques
mu_0= vpa(1.25*sym(1e-6));
sigma_gal= vpa(3.46*sym(1e6));
visco_gal= vpa(0.0024); %Viscosité dynamique du Galstn en Pa.s
eta= 1/(sigma_gal * mu_0); %diffusivité mhd
rho_gal= vpa(6440);% masse volumique Galinstan kg/m3
nu_gal= visco_gal/rho_gal;

% Paramètres géométriques/expérimentaux
B0= vpa(10); %champ magnetique uniforme en Tesla
diam_elec= vpa(1e-3); % en mètre
rayon_elec= diam_elec/2;
dist_inter_elec= 15*10^-3; % en mètre
h= vpa(0.1); % distance entre la plaque inf et sup en mètre

% Grandeurs caractéristiques
L_carac=h;%(eta/omega_0)^(1/2);%longueur caractéristique ==> hauteur de la cuve
tau_eta= L_carac^2/eta;
tau_nu= L_carac^2/nu_gal;
I0=3.18; % Injection courant en Ampère
j0= I0/(L_carac)^2; %%% forme de jz = j0/(pi*sigma_r^2)*exp(-(r/sigma_r)^2)
b0= mu_0*I0/(L_carac); %%% forme de b_theta= b0*(1-exp(-((r)/sigma_r).^2))./(r);
E0= eta/L_carac * b0;

% Paramètres adimensionnées
sigma_r= rayon_elec/L_carac;% diam_adm de l'électrode
hadm=h/L_carac;


% Paramètre tracé
Rmax= hadm;%r_exp_adm;
nr= 5*10^2;%floor(Rmax);
%nz=500;

% Nombre adm
va_num= vpa(B0/sqrt(rho_gal*mu_0)); %Vitesse d'Alfven 
S_num= vpa(va_num*L_carac/eta); %Nombre de Lundquist basé sur la hauteur récipient
Ha_num= vpa(B0*L_carac*sqrt(sigma_gal/visco_gal));

Tau_Alfven= L_carac/va_num;
Freq_Alfven= 1/Tau_Alfven;

%%%% Controle parameter: frequency %%%%
frequence_forcage= Freq_Alfven;% ; %en Hz %0.01*tau^-1
omega= vpa(2*sym(pi)*frequence_forcage);
%epsilon= 100;%tau_eta*omega;
Pm_num= vpa(nu_gal/(eta)); %nb de Prandtl magnétique 

%%% New control parameter
Rnu_fixed= 0;
Rnu_varying_only= 0;
Reta= tau_eta*2*vpa(pi)*frequence_forcage;%
frequence_forcage= Reta/(tau_eta*2*vpa(pi));
Rnu= tau_nu*2*vpa(pi)*frequence_forcage;
epsilon= Reta;
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
%% Modélisation du champ magnétique induit par une série de Fourier-Bessel
% ordre_B= 1;
% % champ induit pour un r positif
% b_theta= 1/(2*pi)*(1-exp(-((r)/sigma_r).^2))./(r); %adm par b0
% b_theta(1)= 0;
% 
% % Calcul de la transformée de Hankel du champ induit
% kt_adm= linspace(kt_min,kt_max,nb_kt);
% TF_b_theta_kt=(1- 0.5*sigma_r*sqrt(pi)*kt_adm.*exp(-sigma_r^2*kt_adm.^2/8)...
%     .*besseli(0.5+ordre_B,sigma_r^2*kt_adm.^2/8)); %adimensionnement par L_carac^2/b0* 


%% Calcul des solutions de la relation de dispersion full MHD 
%(avec solution de Hartmann et Alfven)

% entrée des paramètres d'intérêt : 
% r_min : distance minimale en mètre depuis le centre de l'électrode
% r_max : distance maximale en mètre depuis le centre de l'électrode
% nb_point_r= nombre de points de calcul
% nb_kt= nombre de mode propre transversaux à utiliser pour l'approximation
%       du forcage magnetique (nb_kt est adimensionné par L_carac= 0.1m)
% delta_kt= pas entre deux modes transversaux
r_min= 0.005;%0.2;%
r_max= 0.6;%0.2;
nb_point_r= 150;%150;
h_max= hadm;
nb_point_z= 50  ;
nb_kt= 1800;%1800;%3200;%1000;%2100;
delta_kt=1.2; %2.25; %% %0.5;
type_forcing= 'Gaussian_deriv';%'linear'

if ~exist('r_min') || isempty(r_min)
    r_min= input('Write the minimum scaled location from the electrode to test : \n');
end
if ~exist('r_max') || isempty(r_max)
    r_max= input('Write the maximal scaled location from the electrode to test : \n');
end
if ~exist('nb_point_r') || isempty(nb_point_r)
    nb_point_r= input('Write the number of points along r: \n');
end
if ~exist('nb_point_z') || isempty(nb_point_z)
    nb_point_z= input('Write the number of points along z: \n');
end
if ~exist('nb_kt') || isempty(nb_kt)
nb_kt= input('Write the number of transverse mode to test: \n');
end
if ~exist('delta_kt') || isempty(delta_kt)
delta_kt= input('Value of the discretisation for kt ?\n'); % exemple, 355   
end

ordre_B= 1;
% location and transverse mode vectors
J1roots= vpa(besselzero(1,nb_kt,1)); %determine the roots of J1 in order to find the transverse mode
R= (J1roots(2)-J1roots(1))/vpa(delta_kt); % taille du domaine de définition du forçage magnétique
kt_adm= J1roots/R; %mode transversaux
r_min_adm= vpa(r_min);
r_max_adm= vpa(r_max);
%r_investigation= vpa(logspace(log10(r_min_adm),log10(r_max_adm),nb_point_r));
r_investigation= vpa(linspace(r_min_adm,r_max_adm,nb_point_r));
fprintf('scaled step of location Delta_r = %4.3f \n', (r_max_adm-r_min_adm)/(nb_point_r-1));

% plage frequence
if ~exist('frequence_forcage') || isempty(frequence_forcage)
    frequence_forcage= input('Write the forcing frequency (Hz) ?\n');
    omega= vpa(2*sym(pi)*frequence_forcage);
    epsilon= tau_eta*omega;
end



%initialisation matrice
ratio_ampl_mat= vpa(zeros(1,1));
phase_up_mat= vpa(zeros(1,1));
phase_bot_mat= vpa(zeros(1,1));
dephasage_mat= vpa(zeros(1,1));


%% ampl/phase shift coefficient calculation


% solution for the bottom plate
disp("début du calcul des solutions d'alfven pour un balayage de localisation r") 

% initialisation vecteur signal max top/bottom
% E_max_top= [];
% E_max_bot= [];
E_env= [];
E_phase= [];
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

% Calcul des coefficients Bb_ri, projections du forcage magnetique suivant
% les modes propres transversaux

% syms zeta
% Hankel_Transform_approx= vpaintegral((1-exp(-((zeta)/sigma_r).^2))... 
%     .*besselj(1,J1roots*zeta/R), 0, R,'ArrayValued',true); %Hankel transform de 0 à R
% TF_b_theta= Hankel_Transform_approx;
TF_b_theta= Hankel_T_type_forcing(type_forcing,kt_adm,R,sigma_r,precision);
%TF_b_theta= (2*vpa('pi'))^-1*(1./kt_adm- 0.5*sigma_r*sqrt(vpa('pi')).*exp(-sigma_r^2*kt_adm.^2/8)...
%    .*besseli(0.5,sigma_r^2*kt_adm.^2/8)); %adimensionné par L_carac*b0

Bb_ri= 2./(R*besselj(2,J1roots)).^2 .*TF_b_theta;

disp('Calculating alfven wave solutions for the bottom plate forcing...')
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
disp('end solution searching for bottom plate forcing')

%
% diminution précision des paramètres
    disp('diminution de la précision des paramètre')
    s1_less= vpa(s1);%,precis_b);
    k1_less= vpa(k1);%,precis_b);
    s2_less= vpa(s2);%,precis_b);
    k2_less= vpa(k2);%,precis_b);
    hadm_less= vpa(hadm);%,precis_b);
    disp('fin diminution précision')


parfor ind= 1:length(r_investigation)  
        
 %  ind   
    ind
    fprintf('progressing %4.2f%% ..\n',ind*100/(length(r_investigation)))
     
    %ri= vpa(r_investigation(ind));
    
    % Ajout des fonctions propres aux coefficient et diminution de la
% précisison
%     A_less= vpa(A.*besselj(ordre_B,kt_adm(1:end-1)'*r_investigation(ind)));%,precis_b);
    B_less= vpa(B.*besselj(ordre_B,kt_adm(1:end-1)'*r_investigation(ind)));%,precis_b);
    V_less= vpa(V.*besselj(ordre_B,kt_adm(1:end-1)'*r_investigation(ind)));%,precis_b);
%    J_less= vpa(J.*besselj(ordre_B,kt_adm(1:end-1)'*r_investigation(ind)));%,precis_b);
%    E_less= vpa(E.*besselj(ordre_B,kt_adm(1:end-1)'*r_investigation(ind)));%,precis_b);
%     Fr_phi= vpa((1- besselj(ordre_B-1,kt_adm(1:end-1)'.*r_investigation(ind)))./kt_adm(1:end-1)');
%     Fr_phi(1)= vpa(sym(0));
%     Phi_less= vpa(V.*Fr_phi,precis_b);
    
% calcul des signaux aux limites (top and bottom)

%     t=[0:0.05*2*pi/epsilon:5*2*pi/(epsilon)];
%     lim= 0.5*length(A(:,1));

%     % champ vitesse
%     [asignal_top,asignal_bottom]=calcul_signal_boundaries(A_less,s1_less,k1_less,...
%            s2_less,k2_less,epsilon_inf, hadm_less,t,precision);
    % gradient de potentiel
%     [vsignal_top_j,vsignal_bottom_j]=calcul_signal_boundaries(V_less,s1_less,k1_less,...
%                s2_less,k2_less,epsilon, hadm_less,t,precision);
    % champ magnétique
%     [bsignal_top,bsignal_bottom]=calcul_signal_boundaries(B_less,s1_less,k1_less,...
%                s2_less,k2_less,epsilon_inf, hadm_less,t,precision);
%     % courant électrique
%     [jsignal_top,jsignal_bottom]=calcul_signal_boundaries(J_less,s1_less,k1_less,...
%                s2_less,k2_less,epsilon_inf, hadm_less,t,precision);
    % champ électrique loi d'ohm
%     [esignal_top_j,esignal_bottom_j]=calcul_signal_boundaries(E_less,s1_less,k1_less,...
%                s2_less,k2_less,epsilon_inf, hadm_less,t,precision);
    
%     % potentiel electrique
%         [phisignal_top_j,phisignal_bottom_j]=calcul_signal_boundaries(Phi_less,s1_less,k1_less,...
%                s2_less,k2_less,epsilon, hadm_less,t,precision);
%    
    disp('debut calcul enveloppe & phase gradient phi')
    [gradV_env_j,gradV_phase_j]=amplitudes_somme_onde(V_less,s1_less,...
    k1_less,s2_less,k2_less,h_max,nb_point_z,precision,epsilon);    
    disp('fin calcul enveloppe & phase gradient phi')

%     disp('debut calcul enveloppe & phase electric field ')
%     [E_env_j,E_phase_j]=amplitudes_somme_onde(E_less,s1_less,...
%     k1_less,s2_less,k2_less,h_max,2,precision,epsilon);       
%     disp('fin calcul enveloppe & phase electric field ')
     
    disp('debut calcul enveloppe & phase magnetic field ')
    [B_env_j,B_phase_j]=amplitudes_somme_onde(B_less,s1_less,...
    k1_less,s2_less,k2_less,h_max,nb_point_z,precision,epsilon);       
    disp('fin calcul enveloppe & phase magnetique field ')

%     disp('debut calcul enveloppe & phase velocity field ')
%     [A_env_j,A_phase_j]=amplitudes_somme_onde(A_less,s1_less,...
%     k1_less,s2_less,k2_less,h_max,nb_point_z,precision,epsilon);       
%     disp('fin calcul enveloppe & phase velocity field ')
% 
    % 
% 

%    gradV_max_top_j= max(vsignal_top_j);
%    gradV_max_bot_j= max(vsignal_bottom_j);
%     E_max_top_j= max(esignal_top_j);
%     E_max_bot_j= max(esignal_bottom_j);
%     phi_max_top_j= max(phisignal_top_j);
%     phi_max_bot_j= max(phisignal_bo"ttom_j);    
    

%potential gradient
gradV_env= [gradV_env ; gradV_env_j];
gradV_phase= [gradV_phase ; gradV_phase_j];

% %electric field
%E_env= [E_env ; E_env_j];
%E_phase= [E_phase ; E_phase_j];

%magnetic field
B_env= [B_env ; B_env_j];
B_phase= [B_phase ; B_phase_j];

% %velocity field
% A_env= [A_env ; A_env_j];
% A_phase= [A_phase ; A_phase_j];


s1_struc(ind).value= s1_less;
k1_struc(ind).value= k1_less;
s2_struc(ind).value= s2_less;
k2_struc(ind).value= k2_less;
V_struc(ind).value= V_less;
% E_struc(ind).value= E_less;
B_struc(ind).value= B_less;
% A_struc(ind).value= A_less;

%     gradV_top= [gradV_top ; vsignal_top_j];
%     gradV_bot= [gradV_bot ; vsignal_bottom_j];
% %     E_max_top= [E_max_top ; E_max_top_j];
% %     E_max_bot= [E_max_bot ; E_max_bot_j];
%      phi_top= [phi_top ; phisignal_top_j];
%      phi_bot= [phi_bot ; phisignal_bottom_j];

end

%% calculation signal at a given time outside the Ha layer (10 time delta_Ha)
Ha_num= S_num*sqrt(Pm_num^-1);                      
z_Op= 10*Ha_num^-1;
z_1m= 1-10*Ha_num^-1;
t= [0:0.02*2*pi/epsilon:2*pi/(epsilon)];
z_vect=0:h_max/(40-1):h_max;

gradV_time_struc= struct;   
B_time_struc= struct;
%% calculation of the spatial time evolution
gradV_time_z0p= zeros(length(double(r_investigation)),length(double(t)));
gradV_time_z1m= zeros(length(double(r_investigation)),length(double(t)));
gradV_time_zi= zeros(length(double(r_investigation)),length(double(t)));
B_time_z0p= zeros(length(double(r_investigation)),length(double(t)));
B_time_z1m= zeros(length(double(r_investigation)),length(double(t)));
B_time_zi= zeros(length(double(r_investigation)),length(double(t)));
%E_time_z0p1= zeros(length(r_investigation),length(t));
%E_time_z0p9= zeros(length(r_investigation),length(t));

ind_zi= 1;

s1=  s1_struc(1).value;           
k1= k1_struc(1).value;
s2=  s2_struc(1).value;
k2= k2_struc(1).value;


for zi= z_vect
    parfor ind_r= 1:length(r_investigation)
    
    
        V=  V_struc(1,ind_r).value;
        B=  B_struc(1,ind_r).value;
      %    E= E_struc(ind).value;
        phase_num= vpa('0');
        if ind_zi == 1
            ind_r
%             gradV_time_z0p(ind,:)= calcul_time_signal_given_z(V,s1,k1,s2,k2,epsilon,...
%                 z_Op,t,precision,phase_num);
%             gradV_time_z1m(ind,:)= calcul_time_signal_given_z(V,s1,k1,s2,k2,epsilon,...
%                 z_1m,t,precision,phase_num);
%             B_time_z0p(ind,:)= calcul_time_signal_given_z(B,s1,k1,s2,k2,epsilon,...
%                 z_Op,t,precision,phase_num);
%             B_time_z1m(ind,:)= calcul_time_signal_given_z(B,s1,k1,s2,k2,epsilon,...
%                 z_1m,t,precision,phase_num);
            
        end
        ind_r
        gradV_time_zi(ind_r,:)= calcul_time_signal_given_z(V,s1,k1,s2,k2,epsilon,...
            zi,t,[],phase_num);
        B_time_zi(ind_r,:)= calcul_time_signal_given_z(B,s1,k1,s2,k2,epsilon,...
            zi,t,[],phase_num);
        ind_r
    end
    gradV_time_struc(ind_zi).z_value= gradV_time_zi;
    B_time_struc(ind_zi).z_value= B_time_zi;
    fprintf('progressing %4.2f%% ..\n',100*ind_zi/nb_point_z)
    ind_zi= ind_zi +1;
end
    
%% save time data
selpath = uigetdir(['C:\Users\lalloz\Documents\these\onde alfven\etude theorique numerique'...
    '\1_electrode_model_Results\one plate forcing case\location_sweeping']);
save([selpath,'\time_data_evolution.mat'],'z_Op','z_1m','Ha_num',...
    't','gradV_time_z0p','gradV_time_z1m','B_time_z0p','B_time_z1m','r_investigation',...
    'gradV_time_struc','B_time_struc','z_vect')
  
   
%% save
selpath = uigetdir(['C:\Users\lalloz\Documents\these\onde alfven\etude theorique numerique'...
    '\1_electrode_model_Results\one plate forcing case']);
%%
save(['donnees_nb_z_50.mat'],'Bb_ri','gradV_env',... %_nb_z_2 %selpath,'\
    'gradV_phase','E_env','E_env','B_env','B_phase',...
    'A_env','A_phase','Freq_Alfven','frequence_forcage','Ha_num',...
    'R','S_num','B0','r_investigation','r_min','r_max','nb_point_r',...
    'sigma_r','epsilon','phase_bot_mat','hadm','dist_inter_elec',...
   'J1roots','dist_inter_elec','kt_adm','delta_kt','nb_kt','s1_struc',...
    'k1_struc','s2_struc','k2_struc','V_struc','E_struc','B_struc','A_struc',...
    'Rnu','Reta','nb_point_z','precision')


%% amplitude mapping (r,z)
z_vect= double(0:h_max/(nb_point_z-1):h_max);
[R_mat , Z_mat]= meshgrid(double(r_investigation),z_vect);      
nb_line_env= 40;
nb_line_phase= 20;

%% figure showing the spatial time evolution of the B and gradV

figure('Name','B_spatial_time_evolution')
B_max= max(double(B_time_struc(1).z_value),[],'all');
B_min= min(double(B_time_struc(1).z_value),[],'all');
linecontour_B_time= linspace(B_min/2,B_max/2,60);
for i_time= 1:length(t)
    B_time_i= [];
    for j= 1:length(z_vect)
        B_time_i= [B_time_i, B_time_struc(j).z_value(:,i_time)];
    end
    contourf(R_mat,Z_mat,B_time_i',linecontour_B_time)
    colormap(BWR2)
    c= colorbar;
    caxis([B_min, B_max]);
    c.Limits= ([B_min,B_max]);
    pause(0.2)
end
%%
fig_gradV_sp_time_evolution= figure('Name','gradV_spatial_time_evolution');
gradV_max= 500;%max(double(gradV_time_struc(1).z_value),[],'all');
gradV_min= -500;%min(double(gradV_time_struc(1).z_value),[],'all');
linecontour_gradV_time= linspace(-500,500,100);
for i_time= 1:length(t)
    gradV_time_i= [];
    for j= 1:length(z_vect)
        gradV_time_i= [gradV_time_i, gradV_time_struc(j).z_value(:,i_time)];
    end
    gradV_time_i(gradV_time_i>=gradV_max)= gradV_max;
    gradV_time_i(gradV_time_i<=gradV_min)= gradV_min;
    contourf(R_mat,Z_mat,gradV_time_i',linecontour_gradV_time)
    colormap(BWR2);
    c= colorbar;
    caxis([-500, 500]);
    %c.Limits= ([-500,500]);
    c.Label.String = '$\nabla\phi_r\left(r,z\right)$';
    c.Label.Interpreter= 'latex';
    c.Label.FontSize= 14;
    c.TickLabelInterpreter='latex';
     f= gca;
     f.ColorScale= 'linear';
    f.TickLabelInterpreter= 'latex';
    xlabel('$r$')
    ylabel('$z$')
    f.FontSize= 12;
    f.XLabel.FontSize= 14;
    f.YLabel.FontSize= 14;
    pause(0.2)
   video_struc(i_time)= getframe(fig_gradV_sp_time_evolution);
    hold off
end



%% sauvegarde video

    for i=1:length(video_struc)%/2-1
        video_struc_bis(i)= video_struc(i);
        i
    end
    v= VideoWriter('gradV_along_z_and_r_time_evolution.avi','Uncompressed AVI');
    v.FrameRate= 5;
    %v.Quality= 100;
    %v.CompressionRatio= 1;
    open(v)
    writeVideo(v,video_struc_bis)
    close(v)
    clear v



%%
B_test= B_env;
B_test(B_env == 0)=NaN;

transit_time= (gradV_phase-gradV_phase(:,1))'/double(epsilon);
phase_diff= (gradV_phase-gradV_phase(:,1))';
%A_test= A_env;
%A_test(A_env == 0)=NaN;
linecontour_B_env= logspace(-4,... %log10(min(double(B_test),[],'all',"omitnan"))
    log10(max(double(B_test),[],'all')),nb_line_env);
linecontour_gradV_env= logspace(log10(1e-10),...
    log10(5e5),nb_line_env);
linecontour_gradV_env_norm = logspace(log10(1e-11),...
    log10(1),nb_line_env);
linecontour_B_phase= linspace(min(double(B_phase),[],'all',"omitnan"),...
    max(double(B_phase),[],'all'),nb_line_phase);
linecontour_gradV_phase= linspace(min(double(gradV_phase),[],'all',"omitnan"),...
    max(double(gradV_phase),[],'all'),nb_line_phase);

linecontour_attenuation_coef_env= linspace(-9,1,20);%real([linspace(-12,...
%-0.1,15), linspace(0.1,1,10)]);%4

%linecontour_A_env= logspace(-2,...
%    2,nb_line_env);

% limits for the colorbar
b_ampl_inf= 100*min(abs(double(B_env(:,1:end-1))),[],'all');
b_ampl_sup= max(abs(double(B_env(:,1:end-1))),[],'all');
gradV_ampl_inf= 1;%min(abs(double(gradV_env(:,1:end-1))),[],'all');
gradV_ampl_sup= 10^4;%max(abs(double(gradV_env(:,1:end-1))),[],'all');
%A_ampl_inf= min(abs(double(A_env(:,1:end-1))),[],'all');
%A_ampl_sup= max(abs(double(A_env(:,1:end-1))),[],'all');

depth_penetration_eta= 1/sqrt(vpa(pi)*frequence_forcage/eta);
depth_penetration_eta_adm= depth_penetration_eta/h;


%%% max gradient phi, B0= 10T
% max_f_56Hz= 3.70e5
% max_f_56Hz= 6.44e4
% max_f_560Hz= 2.17e4
% max_f_1114Hz= 1.55e4

%%% min gradient phi, B0= 10T
% min_f_56Hz= 3.70e5
% min_f_56Hz= 6.44e4
% min_f_560Hz= 0.34
% min_f_1114Hz= 1.55e4
%% figure B
figure('Name','b amplitude vs r and z')
s=contourf(R_mat,Z_mat,B_env',linecontour_B_env);
%hold on
%plot(1/sqrt(Ha_num)*ones(nb_point_z,1),z_vect,'--k','linewidth',1.5)
% hold on
% plot(double(r_investigation),sqrt(2*pi*eta/(frequence_forcage))/h*ones(nb_point_r,1),':k','linewidth',1.5)
colormap(WR2)%colormap_vect)
c= colorbar;
c.Label.String = '$|b_{\theta}\left(r,z\right)|$';
c.Label.Interpreter= 'latex';
c.TickLabelInterpreter='latex';
c.Label.FontSize= 14;
f= gca;
f.ColorScale= 'log';
caxis([100*min(abs(double(B_env(:,1:end-1))),[],'all'), max(abs(double(B_env(:,1:end-1))),[],'all')])
%c.Limits= ([min(double(B_env(:,1:end-1)),[],'all'),max(double(B_env),[],'all')]);
f.TickLabelInterpreter= 'latex';
xlabel('$r$')
ylabel('$z$')
f.FontSize= 12;
f.XLabel.FontSize= 14;
f.YLabel.FontSize= 14;
title(sprintf(['Mapping of the magnetic perturbation amplitude for\n $Ha= %4.2e,\\, R_{\\eta}/2\\pi = %4.2f,',...
           ' \\,R_\\nu/2\\pi = %4.2e$'],Ha_num,epsilon/(2*vpa(pi)),Rnu/(2*vpa(pi))))
 
%%
figure('Name','b phase vs r and z');
s=contourf(R_mat(1:end-1,:),Z_mat(1:end-1,:),transpose(B_phase(:,1:end-1)),20)%linecontour_B_phase);
colormap(BWR2)
c= colorbar;
c.Label.String = '$\varphi\left[b_{\theta}\left(r,z\right)\right]$';
caxis([-max(abs(double(B_phase(:,1:end-1))),[],'all'), max(abs(double(B_phase(:,1:end-1))),[],'all')]);
c.Limits= ([min(double(B_phase(:,1:end-1)),[],'all'),max(double(B_phase(:,1:end-1)),[],'all')]);
c.Label.Interpreter= 'latex';
c.TickLabelInterpreter='latex';
c.Label.FontSize= 14;
c.Ticks= [-6*pi  -5*pi -4*pi -3*pi -5*pi/2 -2*pi -7*pi/4 -3*pi/2 -5*pi/4 -pi...
    -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi];
c.TickLabels= {'$-6\pi$', '$-5\pi$','$-4\pi$', '$-3\pi$', '$-5\pi/2$','$-2\pi$', '$-7\pi/4$', '$-3\pi/2$',...
    '$-5\pi/4$','$-\pi$', '$-3\pi/4$','$-\pi/2$', '$-\pi/4$', '$0$', '$\pi/4$',...
    '$\pi/2$', '$3\pi/4$', '$\pi$'};
f= gca;
f.TickLabelInterpreter= 'latex';
f.Color
xlabel('$r$')
ylabel('$z$')
f.FontSize= 12;
f.XLabel.FontSize= 14;
f.YLabel.FontSize= 14;
title(sprintf(['Mapping of the magnetic perturbation phase for\n $Ha= %4.2e,\\, R_{\\eta}/2\\pi = %4.2f,',...
           ' \\,R_\\nu/2\\pi = %4.2e$'],Ha_num,epsilon/(2*vpa(pi)),Rnu/(2*vpa(pi))))
 
ylim([0, 1])
       
%% figure gradV
figure('Name','gradV amplitude vs r and z')
s=contourf(R_mat(:,1:end),Z_mat(:,1:end),transpose(gradV_env(1:end,:)),linecontour_gradV_env);
hold on
plot(1/sqrt(Ha_num)*ones(nb_point_z,1),z_vect,'--k','linewidth',1.5)
hold on
plot(double(r_investigation),depth_penetration_eta_adm*ones(nb_point_r,1),':w','linewidth',1.5)
colormap(WR2(40))
c= colorbar;
caxis([min(gradV_env,[],"all"),max(gradV_env,[],"all")]);
%c.Limits= ([min(double(B_phase(:,1:end-1)),[],'all'),max(double(B_phase(:,1:end-1)),[],'all')]);
c.Label.String = '$|\nabla\phi_r\left(r,z\right)|$';
c.Label.Interpreter= 'latex';
c.Label.FontSize= 14;
c.TickLabelInterpreter='latex';
f= gca;
f.ColorScale= 'log';
f.TickLabelInterpreter= 'latex';
xlabel('$r$')
ylabel('$z$')
f.FontSize= 12;
f.XLabel.FontSize= 14;
f.YLabel.FontSize= 14;
title(sprintf(['Mapping of the radial potential gradient amplitude for\n $Ha= %4.2e,\\, R_{\\eta}/2\\pi = %4.2f,',...
           ' \\,R_\\nu/2\\pi = %4.2e$'],Ha_num,epsilon/(2*vpa(pi)),Rnu/(2*vpa(pi))))
 

%% figure gradV normalised max
figure('Name','gradV amplitude vs r and z')
s=contourf(R_mat(:,1:end),Z_mat(:,1:end),transpose(gradV_env(1:end,:)/max(gradV_env,[],"all")),linecontour_gradV_env_norm);
hold on
plot(1/sqrt(Ha_num)*ones(nb_point_z,1),z_vect,'--k','linewidth',1.5)
hold on
plot(double(r_investigation),depth_penetration_eta_adm*ones(nb_point_r,1),':w','linewidth',1.5)
colormap(WR2(40))
c= colorbar;
caxis([1e-4,1]);
%c.Limits= ([min(double(B_phase(:,1:end-1)),[],'all'),max(double(B_phase(:,1:end-1)),[],'all')]);
c.Label.String = '$|\nabla\phi_r\left(r,z\right)|/\max\left(|\nabla\phi_r|\right)$';
c.Label.Interpreter= 'latex';
c.Label.FontSize= 14;
c.TickLabelInterpreter='latex';
f= gca;
f.ColorScale= 'log';
f.TickLabelInterpreter= 'latex';
xlabel('$r$')
ylabel('$z$')
f.FontSize= 12;
f.XLabel.FontSize= 14;
f.YLabel.FontSize= 14;
title(sprintf(['Mapping of the radial potential gradient amplitude for\n $Ha= %4.2e,\\, R_{\\eta}/2\\pi = %4.2f,',...
           ' \\,R_\\nu/2\\pi = %4.2e$'],Ha_num,epsilon/(2*vpa(pi)),Rnu/(2*vpa(pi))))
 
       
  %%     
figure('Name','gradV phase vs r and z')
contourf(R_mat,Z_mat,gradV_phase',linecontour_gradV_phase);
colormap(BWR2(50))
c= colorbar;
c.Label.String = '$\varphi\left[b_{\theta}\left(r,z\right)\right]$';
caxis([-max(abs(double(gradV_phase)),[],'all'), max(abs(double(gradV_phase)),[],'all')]);
c.Limits= ([min(double(gradV_phase),[],'all'),max(double(gradV_phase),[],'all')]);
c.Label.String = '$\varphi\left[\nabla\phi_r\left(r,z\right)\right]$';
c.Label.Interpreter= 'latex';
c.Label.FontSize= 14;
c.Ticks= [-2*pi -7*pi/4 -3*pi/2 -5*pi/4 -pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi];
c.TickLabels= {'$-2\pi$', '$-7\pi/4$', '$-3\pi/2$', '$-5\pi/4$','$-\pi$', '$-3\pi/4$',...
    '$-\pi/2$', '$-\pi/4$', '$0$', '$\pi/4$', '$\pi/2$', '$3\pi/4$', '$\pi$'};
c.TickLabelInterpreter='latex';
f= gca;
f.TickLabelInterpreter= 'latex';
xlabel('$r$')
ylabel('$z$')
f.FontSize= 12;
f.XLabel.FontSize= 14;
f.YLabel.FontSize= 14;
title(sprintf(['Mapping of the radial potential gradient phase for\n $Ha= %4.2e,\\, R_{\\eta}/2\\pi = %4.2f,',...
           ' \\,R_\\nu/2\\pi = %4.2e$'],Ha_num,epsilon/(2*vpa(pi)),Rnu/(2*vpa(pi))))
 
%% transit time 

figure('Name','gradV transit time vs r and z')
contourf(R_mat,Z_mat,transit_time,15);%linecontour_gradV_phase);
hold on
%contour(R_mat,Z_mat,transit_time,[0 0],'k','linewidth',2);
colormap(BWR2(100))
c= colorbar;

caxis([-max(abs(transit_time),[],'all'), max(abs(transit_time),[],'all')]);
c.Limits= ([min(double(transit_time),[],'all'),max(double(transit_time),[],'all')]);
%c.Limits= ([min(double(transit_time),[],'all'),0.01])
c.Label.String = '$\tau\left(r,z\right)$';
c.Label.Interpreter= 'latex';
c.Label.FontSize= 14;
%c.Ticks= [-2*pi -7*pi/4 -3*pi/2 -5*pi/4 -pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi];
%c.TickLabels= {'$-2\pi$', '$-7\pi/4$', '$-3\pi/2$', '$-5\pi/4$','$-\pi$', '$-3\pi/4$',...
%    '$-\pi/2$', '$-\pi/4$', '$0$', '$\pi/4$', '$\pi/2$', '$3\pi/4$', '$\pi$'};
c.TickLabelInterpreter='latex';
f= gca;
f.TickLabelInterpreter= 'latex';
xlabel('$r$')
ylabel('$z$')
f.FontSize= 12;
f.XLabel.FontSize= 14;
f.YLabel.FontSize= 14;
title(sprintf(['Mapping of the transit time for \n $Ha= %4.2e,\\, R_{\\eta}/2\\pi = %4.2f,',...
           ' \\,R_\\nu/2\\pi = %4.2e$'],Ha_num,epsilon/(2*vpa(pi)),Rnu/(2*vpa(pi))))
 
%% phase diff for gradV

figure('Name','gradV phase diff vs r and z')
contourf(R_mat,Z_mat,phase_diff,15);%linecontour_gradV_phase);
hold on
contour(R_mat,Z_mat,transit_time,[0 0],'k','linewidth',2);
colormap(BWR2(60))
c= colorbar;

caxis([-max(abs(phase_diff),[],'all'), max(abs(phase_diff),[],'all')]);
c.Limits= ([min(double(phase_diff),[],'all'),max(double(phase_diff),[],'all')]);
%c.Limits= ([min(double(transit_time),[],'all'),0.01])
c.Label.String= ['$\varphi\left[\nabla\phi_r\left(r,z\right)\right]',...
    ' - \varphi\left[\nabla\phi_r\left(r,z=0\right)\right]$'];
c.Label.Interpreter= 'latex';
c.Label.FontSize= 14;
c.Ticks= [-6*pi  -5*pi -4*pi -3*pi -5*pi/2 -2*pi -7*pi/4 -3*pi/2 -5*pi/4 -pi...
    -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi];
c.TickLabels= {'$-6\pi$', '$-5\pi$','$-4\pi$', '$-3\pi$', '$-5\pi/2$','$-2\pi$', '$-7\pi/4$', '$-3\pi/2$',...
    '$-5\pi/4$','$-\pi$', '$-3\pi/4$','$-\pi/2$', '$-\pi/4$', '$0$', '$\pi/4$',...
    '$\pi/2$', '$3\pi/4$', '$\pi$'};
c.TickLabelInterpreter='latex';
f= gca;
f.TickLabelInterpreter= 'latex';
xlabel('$r$')
ylabel('$z$')
f.FontSize= 12;
f.XLabel.FontSize= 14;
f.YLabel.FontSize= 14;
title(sprintf(['Mapping of the phase difference for \n $Ha= %4.2e,\\, R_{\\eta}/2\\pi = %4.2f,',...
           ' \\,R_\\nu/2\\pi = %4.2e$'],Ha_num,epsilon/(2*vpa(pi)),Rnu/(2*vpa(pi))))
 
           
%% attenauation coeff
nb_pt_r_bis= 20;
figure('Name','attenuation coef vs r and z');
s=contourf(R_mat(:,1:end),Z_mat(:,1:end),log(transpose(gradV_env./gradV_env(1:end,1))),linecontour_attenuation_coef_env);%linecontour_attenuation_coef_env);
hold on
plot(1/sqrt(Ha_num)*ones(nb_point_z,1),z_vect,'--w','linewidth',1.5)
hold on
plot(linspace(r_min,r_max,nb_pt_r_bis),...
    double(depth_penetration_eta_adm*ones(nb_pt_r_bis,1)),'.k','linewidth',1.5)
hold on
contour(R_mat,Z_mat,log(transpose(gradV_env./gradV_env(1:end,1))),[0 0],'k','linewidth',2);
colormap(BWR2(100))
c= colorbar;
caxis([-14,14]);
c.Limits=([min(double(log(gradV_env./gradV_env(1:end,1))),[],'all'),...
    max(double(log(gradV_env./gradV_env(1:end,1))),[],'all')]); %[linecontour_attenuation_coef_env(1),linecontour_attenuation_coef_env(end)] ;%
c.Label.String = '$\alpha \left( r, z\right)$';
c.Label.Interpreter= 'latex';
c.Label.FontSize= 14;
c.TickLabelInterpreter='latex';
ax= gca;
ax.ColorScale= 'linear';
ax.TickLabelInterpreter= 'latex';
ax.PlotBoxAspectRatio= [0.6000, 1.0000, 1.0000];
ax.Children
xlabel('$r$')
ylabel('$z$')
ax.FontSize= 12;
ax.XLabel.FontSize= 14;
ax.YLabel.FontSize= 14;
title(sprintf(['Mapping of the radial potential gradient amplitude for \n $Ha= %4.2e,\\, R_{\\eta}/2\\pi = %4.2f,',...
           ' \\,R_\\nu/2\\pi = %4.2e$'],Ha_num,epsilon/(2*vpa(pi)),Rnu/(2*vpa(pi))))
       
       
%% Amplitude ratio

i_lim= length(r_investigation);
gradV_max_top= gradV_env(1:i_lim,end); 
gradV_max_bot= gradV_env(1:i_lim,1);
E_max_top= E_env(1:i_lim,end);  
E_max_bot= E_env(1:i_lim,1); 

gradV_phase_top= unwrap(gradV_phase(1:i_lim,end)); %unwrap
gradV_phase_bot= gradV_phase(1:i_lim,1); 
E_phase_top= unwrap(E_phase(1:i_lim,end));  %unwrap
E_phase_bot= (E_phase(1:i_lim,1)); %unwrap

gradV_ampl_ratio= double(gradV_max_top./gradV_max_bot); 
E_ampl_ratio= double(E_max_top./E_max_bot);

E_dephasage= E_phase_top - E_phase_bot;
gradV_dephasage= gradV_phase_top - gradV_phase_bot;


%gradV_ampl_ratio= double(gradV_max_top./gradV_max_bot);
%E_ampl_ratio= double(E_max_bot./E_max_top);

%% transit time and velocity
joule_time= rho_gal/(sigma_gal*B0^2);
gradV_time_delay= -double(gradV_dephasage./(2*vpa('pi')*double(frequence_forcage)));
E_time_delay= -double(E_dephasage./(2*vpa('pi')*double(frequence_forcage)));
gradV_velocity= double(h./gradV_time_delay);
E_velocity= double(h./E_time_delay);

gradV_velocity_adm= double(gradV_velocity/va_num);
E_velocity_adm= double(E_velocity/va_num);

% 2D time
TwoD_time= joule_time*1.^2;
TwoD_time_norm= TwoD_time*0.5;


%% Figure spatial attenuation coefficient
figure
%yyaxis left
plot(r_investigation,log(gradV_ampl_ratio))
hold on
xlabel('$r$')
ylabel('$s_{\nabla\phi}$')
title(sprintf(['Spatial attenuation coefficient for the potential gradient\n'...
' ($S= %4.1f, \\Omega= %4.2f$)'],double(S_num),2*pi*double(frequence_forcage/Freq_Alfven)))
%txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
%t= text(700,0.4,txt);
%t.Units= 'normalized';
f_att_coef_ax= gca;
f_att_coef_ax.TickLabelInterpreter= 'latex';
f_att_coef_ax.FontSize= 10;
f_att_coef_ax.YLabel.FontSize = 15;
f_att_coef_ax.XLabel.FontSize = 12;
grid on

%% Figure transit time scaled with the alfven time coefficient
figure
%yyaxis left
plot(r_investigation,gradV_time_delay*Freq_Alfven)
hold on
xlabel('$r$')
ylabel('$\tau_{\nabla\phi}/\tau_A$')
title(sprintf(['scaled transit time for the potential gradient\n'...
' ($S= %4.1f, \\Omega= %4.2f$)'],double(S_num),2*pi*double(frequence_forcage/Freq_Alfven)))
%txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
%t= text(700,0.4,txt);
%t.Units= 'normalized';
f_transit_time_ax= gca;
f_transit_time_ax.TickLabelInterpreter= 'latex';
f_transit_time_ax.FontSize= 10;
f_transit_time_ax.YLabel.FontSize = 15;
f_transit_time_ax.XLabel.FontSize = 12;
grid on

%% figure Evolution along r at a specific time
time_loc= 3;

fig_animation= figure;

for i= 1:length(t)
%yyaxis left
plot(r_investigation,gradV_time_z0p(:,i))
hold on
plot(r_investigation,gradV_time_z1m(:,i))
hold on
xlabel('$r$')
ylabel('$\nabla\phi_r (r,z,t)$') %'$b_\theta (r,z,t) $'
title(sprintf(['Magnetic field regarding the location r at a given time\n'...
' ($S= %4.1f ,\\, \\left( 2 \\pi \\right)^{-1}\\Omega= %4.2f$)'],double(S_num),double(frequence_forcage/Freq_Alfven)))
%txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
%t= text(700,0.4,txt);
%t.Units= 'normalized';
legend('$z=0^+$','$z=1^-$','Interpreter','latex') %
f_transit_time_ax= gca;
%f_transit_time_ax.YScale= 'log';
f_transit_time_ax.TickLabelInterpreter= 'latex';
f_transit_time_ax.FontSize= 10;
f_transit_time_ax.YLabel.FontSize = 15;
f_transit_time_ax.XLabel.FontSize = 12;
xlim([double(r_investigation(1)) double(r_investigation(end))])
ylim([-max(double(gradV_time_z0p),[],'all'),max(double(gradV_time_z0p),[],'all')])
grid on
pause(0.1)
video_struc(i)= getframe(fig_animation);
hold off
end



%% sauvegarde video

    for i=1:length(video_struc)%/2-1
        video_struc_bis(i)= video_struc(i);
        i
    end
    v= VideoWriter('b_z0p_and_z1m_vs_r.avi');
    v.FrameRate= 10;
    open(v)
    writeVideo(v,video_struc_bis)
    close(v)
    clear v




%% grad phi
size_ft= 12;
figure
subplot(2,1,1)
 plot(r_investigation,gradV_max_top)
    xlabel('$\tilde{r}/h$')
    ylabel('$|\nabla\phi|/\nabla\phi_0$')
    title('top plate')
    
subplot(2,1,2)
  plot(r_investigation,gradV_max_bot)
    xlabel('$\tilde{r}/h$')
    ylabel('$|\nabla\phi|/\nabla\phi_0$')
    title('bottom plate')
    
sgtitle(sprintf(['Evolution of the potential gradient amplitude \n regarding '...
    'the distance from the electrode \n(S= %4.1f, Frequency= %4.1f Hz, $\\mathrm{E}_0$= %4.1e V/m)'],...
    double(S_num),frequence_forcage,double(E0))) 


%% amplitude ratio

figure

%subplot(2,1,1)  %yyaxis left
plot(r_investigation,gradV_ampl_ratio)
hold on 
%p= plot(r_investigation,mean_gradV_ratio,'b');
%p.LineWidth = 1;
xlabel('r/$\sigma_r$')
ylabel(sprintf('amplitude ratio between \n the top and bottom plate'))

legend('radial \nabla{\phi} ratio')
title(sprintf(['Potentiel gradient amplitude ratio in function of the distance'...
    ' from\n the electrode (S= %4.1f, Frequency = %4.1f Hz)'],double(S_num),frequence_forcage)) 

%yyaxis right
%subplot(2,1,2) 
% %hold on
% plot(r_investigation,gradV_max_bot)
% hold on
% plot(r_investigation,gradV_max_top)
% ylabel(sprintf('amplitude of the \n non scaled potential gradient '))
% xlabel('non scaled distance from the centre of the electrode')
% 
% legend('Amplitude of the bottom gradV','Amplitude of the top gradV')
% title(sprintf(['Potentiel gradient amplitude in function of the distance'...
%     ' from\n the electrode (S= %4.1f)'],double(S_num))) 
% 
% txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0)
% t= text(0.14,400,txt);
% t.Units= 'normalized';
% set(gca, 'YScale', 'log')

% figure
% plot(epsilon_mat,gradV_ampl_ratio)
% figure
% plot(epsilon_mat,gradV_ampl_ratio)
% xlabel('relative error of $b_{approx}$ en percentage')
% ylabel('amplitude ratio between the top and bottom plate')
% legend('potential gradient ratio')

%% Figure ratio d'amplitude de GradV et E en fonction de la localisation de la mesure
% Par rapport au centre de l'électrode

figure
plot(r_investigation,gradV_max_bot./gradV_max_top)
xlabel('non scaled location from the electrode')
ylabel('ampltiude ratio between the bottom and the top value')
title(sprintf(['Potentiel gradient amplitude ratio in function of \n the distance'...
    ' from the electrode (S= %4.1f , nb kt= %i, kt_{max}= %i)'],double(S_num),int16(nb_kt),int16(kt_max)))

figure
plot(r_investigation,E_max_bot./E_max_top)
xlabel('non scaled location from the electrode')
ylabel('ampltiude ratio between the bottom and the top value')
title(sprintf(['Electric field amplitude ratio in function of \n the distance'...
    ' from the electrode (S= %4.1f , nb kt= %i, kt_{max}= %i)'],double(S_num),int16(nb_kt),int16(kt_max)))

%% Calculation of the wave enveloppe and the phases for each field
    z=0:hadm_less/(nz-1):hadm_less;
    disp('debut calcul enveloppe')
        if length(A_less(:,1))==1
%             [venv,vphase,z2]=amplitudes(V_less,s1_less,k1_less,...
%                s2_less,k2_less,hadm_less,nz);
           [benv,bphase,z2]=amplitudes(B_less,s1_less,k1_less,...
               s2_less,k2_less,hadm_less,nz);
%            [aenv,aphase,z2]=amplitudes(A_less,s1_less,k1_less,...
%                s2_less,k2_less,hadm_less,nz);
           [jenv,jphase,z2]=amplitudes(J_less,s1_less,k1_less,...
               s2_less,k2_less,hadm_less,nz);
           [jenv_ohm,jphase_ohm,z2]=amplitudes(J_ohm_less,s1_less,k1_less,...
                s2_less,k2_less,hadm_less,nz);
           disp('fin')
        else
            tic
%             disp('debut calcul pour vitesse')
%             [aenv,aphase,z2]=amplitudes_somme_onde(A_less,s1_less,k1_less,...
%                 s2_less,k2_less,hadm_less,nz,precision,sym(0),epsilon_inf);
            disp('debut calcul pour gradient potentiel')
            [venv,vphase,z2]=amplitudes_somme_onde(V_less,s1_less,k1_less,...
                s2_less,k2_less,hadm_less,nz,precision,sym(0),epsilon);
%             disp('debut calcul pour champ magnetic')
%             [benv,bphase,z2]=amplitudes_somme_onde(B_less,s1_less,k1_less,...
%                 s2_less,k2_less,hadm_less,nz,precision,sym(0),epsilon_inf);
%             disp('debut calcul pour courant')
%             [jenv,jphase,z2]=amplitudes_somme_onde(J_less,s1_less,k1_less,...
%                 s2_less,k2_less,hadm_less,nz,precision,sym(0),epsilon_inf);
            disp('debut calcul champ électrique')
            [eenv,jehase,z2]=amplitudes_somme_onde(E_less,s1_less,k1_less,...
                s2_less,k2_less,hadm_less,nz,precision,sym(0),epsilon);

            toc
        end
    disp('fin calcul enveloppe')
    
% %% calcul de j par calcul numérique de db_theta/dz
% z=0:hadm_less/(nz-1):hadm_less;
% time= [0:0.01*2*vpa(pi)/epsilon_inf:6*vpa(pi)/epsilon_inf]';
% grad_f_b_theta= [];
% grad_g_b_theta= [];
% f_B_grad_ana= [];
% g_B_grad_ana= [];
% f_B= [];
% g_B= [];
% 
% for j=1:nb_kt
%     f_b_theta_j= cos(k1(j).*z).*(B_less(j,1).*exp(s1(j).*z)+B_less(j,3).*exp(-s1(j).*z))...
%            +sin(k1(j).*z).*(-B_less(j,2).*exp(s1(j).*z)+B_less(j,4).*exp(-s1(j).*z))...
%            +cos(k2(j).*z).*(B_less(j,5).*exp(s2(j).*z)+B_less(j,7).*exp(-s2(j).*z))...
%            +sin(k2(j).*z).*(-B_less(j,6).*exp(s2(j).*z)+B_less(j,8).*exp(-s2(j).*z));
% 
%     g_b_theta_j= cos(k1(j).*z).*(-B_less(j,2).*exp(s1(j).*z)-B_less(j,4).*exp(-s1(j).*z))...
%            +sin(k1(j).*z).*(-B_less(j,1).*exp(s1(j).*z)+B_less(j,3).*exp(-s1(j).*z))...
%            +cos(k2(j).*z).*(-B_less(j,6).*exp(s2(j).*z)-B_less(j,8).*exp(-s2(j).*z))...
%            +sin(k2(j).*z).*(-B_less(j,5).*exp(s2(j).*z)+B_less(j,7).*exp(-s2(j).*z));
%             
%        
%     grad_f_b_j= gradient(double(f_b_theta_j),double(hadm_less/(nz-1)));    
%     grad_g_b_j= gradient(double(g_b_theta_j),double(hadm_less/(nz-1))); 
%       
%      
%     f_B= [f_B; f_b_theta_j];
%     g_B= [g_B; g_b_theta_j];
% 
%     grad_f_b_theta= [grad_f_b_theta; grad_f_b_j];
%     grad_g_b_theta= [grad_g_b_theta; grad_g_b_j];
%     
%     f_B_grad_ana= [f_B_grad_ana; f_B_grad_j];
%     g_B_grad_ana= [g_B_grad_ana; g_B_grad_j];
% 
% j
% end

% %% calcul enveloppe pour chaque mode propre
% env_grad_b_mode_num= [];
% env_grad_b_mode_ana= [];
% signal_ana= struct;
% signal_num= struct;
% signal_b= struct;
% 
% for j=1:nb_kt
%     signal_num_j= S_num^2*grad_f_b_theta(j,:).*cos(epsilon_inf*time)+S_num^2*grad_g_b_theta(j,:).*sin(epsilon_inf *time);
%     env_grad_b_num_j= max(abs(hilbert(double(signal_num_j))));
%     
%     signal_ana_j= S_num^2*f_B_grad_ana(j,:).*cos(epsilon_inf*time)+S_num^2*g_B_grad_ana(j,:).*sin(epsilon_inf *time);
%     env_grad_b_ana_j= max(abs(hilbert(double(signal_ana_j))));
%     
%     signal_b_j= S_num^2*f_B(j,:).*cos(epsilon_inf*time)+S_num^2*g_B(j,:).*sin(epsilon_inf *time);
%     env_b_ana_j= max(abs(hilbert(double(signal_ana_j))));
%         
%     signal_num(j).valeur= signal_num_j;
%     signal_ana(j).valeur= signal_ana_j;
%     signal_b(j).valeur= signal_b_j;
%     signal_num(j).kti= kt_adm(j);
%     signal_ana(j).kti= kt_adm(j);
%     
%     env_grad_b_mode_num= [env_grad_b_mode_num ; env_grad_b_num_j];
%     env_grad_b_mode_ana= [env_grad_b_mode_ana ; env_grad_b_ana_j];
% 
% j
% end
% % comparaison évolution temporelle RotB num/ana 
% rot_b_num= signal_num(4).valeur;
% rot_b_ana= signal_ana(4).valeur;
% b_ana= signal_b(4).valeur;

% %% Calcul s^2*grad(b)
% S= double(S_num);
% grad_b_num= S^2*sum(grad_f_b_theta).*cos(epsilon_inf*time)+S^2*sum(grad_g_b_theta).*sin(epsilon_inf *time);
%     env_grad_b_num= max(abs(hilbert(double(grad_b_num))));


% %% comparaison des coeffecients entre les méthodes numériques et analytiques
% figure
% for j= 1:nb_kt
% subplot(4,2,j)
% plot(-env_grad_b_mode_ana(j,2:end-1),z(2:end-1),'r',-env_grad_b_mode_num(j,2:end-1),z(2:end-1),'b')
% hold on
% plot(env_grad_b_mode_ana(j,2:end-1),z(2:end-1),'r',env_grad_b_mode_num(j,2:end-1),z(2:end-1),'b')
% xlabel('z axis')
% ylabel('coefficient')
% 
% title(sprintf('kt= %4.1f',kt_adm(j)))
% end
% sgtitle('enveloppe par mode, forçage inférieur')
% legend('grad b analytique','grad b numérique')
%% regard mode transversal par mode transversal    
  benv_par_mode= [];
  jenv_par_mode= [];
for j = 1:nb_kt
            [benv_par_mode_j,bphase,z2]=amplitudes(B_less(j,:),s1_less(j),k1_less(j),...
            s2_less(j),k2_less(j),hadm_less,nz);
%            [aenv,aphase,z2]=amplitudes(A_less,s1_less,k1_less,...
%                s2_less,k2_less,hadm_less,nz);
            [jenv_par_mode_j,jphase,z2]=amplitudes(J_less(j,:),s1_less(j),k1_less(j),...
               s2_less(j),k2_less(j),hadm_less,nz);
            benv_par_mode= [benv_par_mode ; benv_par_mode_j];
            jenv_par_mode= [jenv_par_mode ;jenv_par_mode_j];
           j
           
end   

%% plot current with the analytical method
figure
 for j= 1:nb_kt
    subplot(2,4,j)
    plot(-jenv_par_mode(j,2:end-1),z(2:end-1),'r',jenv_par_mode(j,2:end-1),z(2:end-1),'r')
    xlabel('current j')
    ylabel('z axis')
    legend(sprintf('kt= %4.1f',kt_adm(j)))
end
sgtitle(sprintf('Current envelope simple forcing (top) \n analytical calculation of rot(b)'))


%% 
figure
 for j= 1:nb_kt
subplot(2,4,j)
plot(-benv_par_mode(j,2:end-1),z(2:end-1),'r',benv_par_mode(j,2:end-1),z(2:end-1),'r')
xlabel('magnetic field')
ylabel('z axis')
legend(sprintf('kt= %4.1f',kt_adm(j)))
end
sgtitle('magnetic field envelope simple forcing (top)')


%% courant j
figure 
plot(-jenv(2:end-1),z(2:end-1),'r',jenv(2:end-1),z(2:end-1),'r')
%hold on
%plot(-jenv_forcage_bot(end-1:-1:2),z(2:end-1),'b',jenv_forcage_bot(end-1:-1:2),z(2:end-1),'b')
title('Current envelope simple forcing (top)')
xlabel('current j')
ylabel('z axis')

% figure 
% plot(-jenv_forcage_bot(2:end-1),z(2:end-1),jenv_forcage_bot(2:end-1),z(2:end-1))
% title('Current envelope')
% xlabel('current j')
% ylabel('z axis')
% legend('coefficent vitesse (solution 1)')
% 
% figure 
% plot(-jenv_forcage_top,z,jenv_forcage_top,z)
% title('Current envelope')
% xlabel('current j')
% ylabel('z axis')
% legend('coefficent champ mgnetique (solution 2)')

%% champ magnetique b
figure 
plot(-benv(1:end),z(1:end),benv(1:end),z(1:end))
title('magnetic field envelope simple forcing (top)')
xlabel('magnetic field b')
ylabel('z axis')

% figure 
% plot(-benv_forcage_bot(2:end-1),z(2:end-1),benv_forcage_bot(2:end-1),z(2:end-1))
% title('magnetic field envelope')
% xlabel('magnetic field b')
% ylabel('z axis')
% legend('coefficent vitesse (solution 1)')
% 
% figure 
% plot(-benv_forcage_top(2:end-1),z(2:end-1),benv_forcage_top(2:end-1),z(2:end-1))
% title('magnetic field envelope')
% xlabel('magnetic field b')
% ylabel('z axis')
% legend('coefficent champ mgnetique (solution 2)')

%% calcul des signaux aux limites (top and bottom)
t=[0:0.05*2*pi/epsilon:5*2*pi/(epsilon)];
lim= 0.5*length(A(:,1));
tic
%     % champ vitesse
%     [asignal_top,asignal_bottom]=calcul_signal_boundaries(A_less,s1_less,k1_less,...
%            s2_less,k2_less,epsilon_inf, hadm_less,t,precision);
    % gradient de potentiel
    [vsignal_top,vsignal_bottom]=calcul_signal_boundaries(V_less,s1_less,k1_less,...
               s2_less,k2_less,epsilon, hadm_less,t,precision);
%     % champ magnétique
%     [bsignal_top,bsignal_bottom]=calcul_signal_boundaries(B_less,s1_less,k1_less,...
%                s2_less,k2_less,epsilon_inf, hadm_less,t,precision);
%     % courant électrique
%     [jsignal_top,jsignal_bottom]=calcul_signal_boundaries(J_less,s1_less,k1_less,...
%                s2_less,k2_less,epsilon_inf, hadm_less,t,precision);
    % champ électrique loi d'ohm
    [esignal_top,esignal_bottom]=calcul_signal_boundaries(E_less,s1_less,k1_less,...
               s2_less,k2_less,epsilon, hadm_less,t,precision);
toc

%%

figure 
plot(t,jsignal_top)
hold on
plot(t,jsignal_bottom)
xlabel('non-scaled time')
ylabel('current')
legend('signal Jtop','signal Jbottom')

figure 
plot(t,asignal_top)
hold on
plot(t,asignal_bottom)
xlabel('non-scaled time')
ylabel('current')
legend('signal velocity top','signal velocity down')

figure 
plot(t,vsignal_top)
hold on
plot(t,vsignal_bottom)
xlabel('non-scaled time')
ylabel('potential gradient')
legend('signal grad(phi) top','signal grad(phi) bottom')

figure 
plot(t,esignal_top)
hold on
plot(t,esignal_bottom)
xlabel('non-scaled time')
ylabel('electric field')
legend('signal E top','signal E bottom')
%%

figure 
plot(t,jsignal_bottom)
hold on
plot(t,jsignal_top)
legend('signal B bottom','signal B top')


%%  
time=[0:0.2*2*vpa(pi)/epsilon:1*2*vpa(pi)/(epsilon)];
length(time)
i=1
asignal= [];
vsignal= [];
bsignal= [];
jsignal= [];
esignal= [];
for t = time
% asignal_inf= calcul_signal_along_z(A_less,s1_less,k1_less,...
%                s2_less,k2_less,epsilon_inf,z,t,precision,0);
vsignal_inf= calcul_signal_along_z(V_less,s1_less,k1_less,...
               s2_less,k2_less,epsilon,z,t,precision,0);
            
% bsignal_inf= calcul_signal_along_z(B_less,s1_less,k1_less,...
%                s2_less,k2_less,epsilon_inf,z,t,precision,0);
% 
% jsignal_inf= calcul_signal_along_z(J_less,s1_less,k1_less,...
%                s2_less,k2_less,epsilon_inf,z,t,precision,0);
           
esignal_inf= calcul_signal_along_z(E_less,s1_less,k1_less,...
               s2_less,k2_less,epsilon,z,t,precision,0);      
           
% j_ohm_signal_inf= calcul_signal_along_z(J_ohm_less,s1_less,k1_less,...
%                s2_less,k2_less,epsilon_inf,z,t,precision,0);

% asignal= [asignal ; asignal_inf ];           
vsignal= [vsignal ; vsignal_inf ];
% bsignal= [bsignal ; bsignal_inf ];
% jsignal= [jsignal ; jsignal_inf ];
esignal= [esignal ; esignal_inf ];

i= i+1
end

%%

set(0,'defaultTextInterpreter','latex')
opengl software
figure
clf


%for i = [1:1:length(time)]
% subplot(1,5,1)
% plot(venv,z,'r',-venv,z,'r')
% hold on
% plot(vsignal(i,:),z,'k');
% hold off
% title('Electric potential')
% xlabel('$V/B_0\eta$')
% ylabel('$z (\omega/\eta)^{1/2}$')

subplot(1,2,1)
plot(venv(2:end-1),z(2:end-1),'r',-venv(2:end-1),z(2:end-1),'r')
hold on
%plot(vsignal(i,:),z,'k');
%hold off
title(sprintf('potential gradient'))
xlabel('$\nabla\,\phi / \nabla\,\phi_0$')
ylabel('$z/L_{carac}$')

subplot(1,2,2)
plot(eenv,z,'r',-eenv,z,'r')
hold on
%plot(esignal(i,:),z,'k');
title('electric field')
xlabel('$E/ E_0$')
ylabel('$z/L_{carac}$')

% subplot(1,5,3)
% plot(benv,z,'r',-benv,z,'r')
% hold on
% % plot(bsignal(i,:),z,'k');
% hold off
% title('magnetic field')
% xlabel('$b/ b_0$')
% ylabel('$z/L_{carac}$')
% 
% subplot(1,5,4)
% plot(jenv(2:end-1),z(2:end-1),'r',-jenv(2:end-1),z(2:end-1),'r')
% hold on
% % plot(jsignal(i,2:end-1),z(2:end-1),'k');
% hold off
% title(sprintf('electric current'))
% xlabel('$j/ j_0$')
% ylabel('$z/L_{carac}$')
% 
% 
% subplot(1,5,5)
% plot(aenv,z,'r',-aenv,z,'r')
% hold on
% % plot(asignal(i,:),z,'k');
% title('velicty field')
% xlabel('$u/ u_0$')
% ylabel('$z/L_{carac}$')

sgtitle(sprintf('Enveloppe of different parameters, \n simple bottom forcing'))
pause(1.0)
i
hold off
%end

          
%% figure enveloppe et phase

% for the electric potential
    if show_fig == 1
    Fig_enveloppe_phase_onde= figure('name','fig_enveloppe_onde');
    subplot(1,2,1)
    %axis([-Bmax Bmax 0 hadm])
    title('enveloppe of Electric potential')
    xlabel('$V/ \left(b_0\,L_{carac}\,\omega\right)$')
    ylabel('$z/L_{carac}$')
    hold on
    plot(venv,z,'r',-venv,z,'r')

    subplot(1,2,2)
    %axis([-Vmax Vmax 0 hadm])
    plot(vphase,z,'k')
    ylabel('$z/L_{carac}$')
    xlabel('Phase to forcing (rad)')
    title('Phase of Electric potential')
    xticks(-pi/2:pi/4:pi/2) 
    xticklabels({'-\pi/2','-\pi/4','0','\pi/4','\pi/2'})
    sgtitle('Electric potential behaviours')
    end 

% for the magnetic field    
    if show_fig == 1
    Fig_enveloppe_phase_onde= figure('name','fig_enveloppe_onde');
    subplot(1,2,1)
    %axis([-Bmax Bmax 0 hadm])
    title('enveloppe of magnetic field')
    xlabel('$b/ \left(b_0)$')
    ylabel('$z/L_{carac}$')
    hold on
    plot(benv,z,'r',-benv,z,'r')

    subplot(1,2,2)
    %axis([-Vmax Vmax 0 hadm])
    plot(bphase,z,'k')
    ylabel('$z/L_{carac}$')
    xlabel('Phase to forcing (rad)')
    title('Phase of magnetic field')
    xticks(-pi/2:pi/4:pi/2) 
    xticklabels({'-\pi/2','-\pi/4','0','\pi/4','\pi/2'})
    sgtitle('magnetic field behaviours')
    end   
    
% extraction ratio ampl / phase
    disp('debut calcul ration ampl et dephasage')
    ratio_ampl= venv(end)/venv(1);
    phase_bot= vphase(1);
    phase_up= vphase(end);
    dephasage= vphase(end)-vphase(1);
    disp('fin calcul ration ampl et dephasage')




