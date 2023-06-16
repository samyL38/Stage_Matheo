qclear all
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


% Les calculs sont effectués à une pulsation et une distance radiale fixe
%% toolbox flottant haute précision
%addpath('C:\Users\lalloz\Documents\Multiprecision Computing Toolbox')
precision= 8*32;% % Setup default precision to 40 decimal digits (quadruple).
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
diam_elec= 1e-3; % en mètre
rayon_elec= diam_elec/2;
dist_inter_elec= 2*10^-3; % en mètre
h= vpa(0.1); % distance entre la plaque inf et sup en mètre
phase_inf= sym(0);

% Grandeurs caractéristiques
L_carac=h;%(eta/omega_0)^(1/2);%longueur caractéristique ==> hauteur de la cuve
tau_eta= L_carac^2/eta;
tau_nu= L_carac^2/nu_gal;
I0=sqrt(2)*2; % Injection courant en Ampère
j0= I0/(L_carac)^2; %%% forme de jz = j0/(pi*sigma_r^2)*exp(-(r/sigma_r)^2)
b0= mu_0*I0/(L_carac); %%% forme de b_theta= b0*(1-exp(-((r)/sigma_r).^2))./(r);
E0= eta/L_carac * b0;


% Paramètres adimensionnées
sigma_r= rayon_elec/L_carac;% diam_adm de l'électrode
hadm=h/L_carac;
L_trans= rayon_elec;

% Paramètre tracé
Rmax= hadm;%r_exp_adm;
nr= 5*10^2;%floor(Rmax);
nz= 20;


% Nombre adm
S_num= B0/sqrt(rho_gal*mu_0)*L_carac/eta; %Nombre de Lundquist basé sur la hauteur récipient
Ha_num= B0*L_carac*sqrt(sigma_gal/visco_gal); %Nombre d'Hartmann basé sur la hauteur récipient
Pm_num= vpa(nu_gal/(eta)); %nb de Prandtl magnétique 

% grandeurs caractéristiques
va_num= B0/sqrt(rho_gal*mu_0); %Vitesse d'Alfven 
Tau_Alfven= L_carac/va_num;
Tau_Ha= h^2/(2*nu_gal*Ha_num); 
Tau_joule= rho_gal/(sigma_gal*B0^2);
Tau_2D= Tau_joule*(h/L_trans)^2;
Freq_Alfven= 1/Tau_Alfven;


%%%% Controle parameter: frequency %%%%
frequence_forcage= [];% ; %en Hz %0.01*tau^-1
% frequence
if ~exist('frequence_forcage') || isempty(frequence_forcage)
    frequence_forcage= input('Write the frequency (in Hz) to test : \n');
end
omega_inf= vpa(2*sym(pi)*frequence_forcage);%390);
epsilon_inf= tau_eta*omega_inf;


%%% New control parameter
Rnu_fixed= 0;
Reta=  tau_eta*2*vpa(pi)*frequence_forcage;
Rnu= tau_nu*2*vpa(pi)*frequence_forcage;
if Rnu_fixed == 1
    Pm_num= Reta/Rnu;
    S_num= Ha_num*sqrt(Reta/Rnu);
    epsilon_inf= Reta;
end 

%% Champ à calculer
% Entrer ici les champs à calculer
calculate_U= [];
calculate_V= [];
calculate_B= [];
calculate_E= [];
calculate_J= [];
calculate_boundary_signal= [];
calculate_time_evolution_along_z= [];
save_data= [];
save_video_time_evol_along_z= [];

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
if ~exist('calculate_boundary_signal') || isempty(calculate_boundary_signal)
    calculate_boundary_signal=input(...
        'Do you want to calculate boundary signals ? Yes (1), No (0) \n');
end
if ~exist('calculate_time_evolution_along_z') || isempty(calculate_time_evolution_along_z)
    calculate_time_evolution_along_z=input(...
        'Do you want to calculate the time evolution along z ? Yes (1), No (0) \n');
end
if ~exist('save_data') || isempty(save_data)
    save_data=input(...
        'Do you want to save data ? Yes (1), No (0) \n');
end
if ~exist('save_video_time_evol_along_z') || isempty(save_video_time_evol_along_z)
    save_video_time_evol_along_z=input(...
        'Do you want to save figure on time evolution along z ? Yes (1), No (0) \n');
end


%% Modélisation du champ magnétique induit par une série de Fourier-Bessel (obsolete)
% ordre_B= 1;
% % champ induit pour un r positif
% b_theta= 1/(2*pi)*(1-exp(-((r)/sigma_r).^2))./(r); %adm par b0
% b_theta(1)= 0;
% 
% % Calcul de la transformée de Hankel du champ induit
% kt_adm= linspace(kt_min,kt_max,nb_kt);
% TF_b_theta_kt=(1- 0.5*sigma_r*sqrt(pi)*kt_adm.*exp(-sigma_r^2*kt_adm.^2/8)...

%     .*besseli(0.5*ordre_B,sigma_r^2*kt_adm.^2/8)); %adimensionnement par L_carac^2/b0* 

%% obsolete
% % Calcul approx champ induit par série de Bessel
% loc_sigma_L= length(kt_adm);
% delta_kt= (kt_max-kt_min)/(nb_kt-1)*vpa(ones(1,loc_sigma_L));
% 
% %%% Calcul des termes du champ magnetique imposé pour chaque mode transversal %%%
% %%% Attention les fonctions propres (fonction de Bessel) ne sont pas          %%%  
% %%% introduites dans le calcul terme à terme : cela en vu de calculer les     %%% 
% %%% coefficient de la vitesse, du champ etc, libéré des fonctions propres.    %%% 
% %%% Les fonctions propres sont introduites après le calcul des coeffcient     %%% 
% %%% A,B, etc (voir partie "Ajout des fonctions propres")                      %%% 
% %%% Cela permet de calculer le potentiel électrique directement dans la       %%% 
% %%% fonction "alfven non homogene". En effet le potentiel électrique ne       %%% 
% %%% s'exprime pas dans la même base que les autres champs !!                  %%% 
% 
% B= 1/(2*pi)*(TF_b_theta_kt.*delta_kt)'.*ones(nb_kt,nz);%besselj(ordre_B,kt_adm'.*r);
% % B_bessel= 1/(2*pi)*(TF_b_theta_kt.*delta_kt)'.*besselj(ordre_B,kt_adm'.*r);
% % Bb_tot= sum(B_bessel);
% % coef_corr= max(b_theta)/max(Bb_tot);
% % Bb_c= coef_corr*B;
% % Bb_tot_c= coef_corr*Bb_tot;
% disp('fin')

%% Définition of the boundary magnetic perturbation
%(avec solution de Hartmann et Alfven)
% i_ri= find(r >= r_exp_adm,1);%r_exp_adm; % longueur par rapport au centre de l'électrode
% ri=r(i_ri);
% Bb_c_ri= vpa(1);%B(:,i_ri);%Bb_c(:,i_ri);
%Bb_ri_bessel = B_bessel(:,i_ri);

% parameters by default
nb_point_z= [];
nb_kt= 20;%1200;%2500;%3200;%1000;%2100;%
r_investigation= []; %location scaled with h
kt_adm= [];
R= 1.45;
%delta_kt= 2.25;%0.5;%



cas_single_mode= input('Study of a single transverse mode case ? Yes (1), No (0) \n');

switch cas_single_mode
    
    case 0
        if ~exist('r_investigation') || isempty(r_investigation)
            r_investigation= vpa(input('Write the radial distance from the electrode, scaled with the height: \n'));
        end
        if ~exist('nb_kt') || isempty(nb_kt)
            nb_kt= input('Write the number of transverse mode to test: \n');
        end
        if ~exist('R') || isempty(R)
            R= input('Write the radial size of the vessel (scaled with h) ?\n'); % exemple, 355   
        end
        if ~exist('nb_point_z') || isempty(nb_point_z)
            nb_point_z= input('Write the number of points along the z axis ?\n'); %1 exemple, 355   
        end
%         if ~exist('delta_kt') || isempty(delta_kt)
%             delta_kt= input('Value of the discretisation for kt ?\n'); %1 exemple, 355   
%         end
        ordre_B= 1;
        % location and transverse mode vectors
        J1roots= vpa(besselzero(1,nb_kt,1)); %determine the roots of J1 in order to find the transverse mode
        %R= (J1roots(2)-J1roots(1))/vpa(delta_kt); % taille du domaine de définition du forçage magnétique
        kt_adm= J1roots/vpa(R); %mode transversaux
        fprintf('Value of the domain R %4.2f \n',double(R))

        % Calcul des coefficients Bb_ri, projections du forcage magnetique suivant
        % les modes propres transversaux
        TF_b_theta_dot_kt= double((1- 0.5*kt_adm*sigma_r*sqrt(vpa(pi)).*exp(-sigma_r^2*kt_adm.^2/8)...
           .*besseli(0.5,sigma_r^2*kt_adm.^2/8))); %adimensionné par L_carac*b0
       Bb_ri=(kt_adm)./(pi*(double(J1roots.*besselj(2,J1roots))).^2).*TF_b_theta_dot_kt; 
       
    case 1
        if ~exist('kt_adm') || isempty(kt_adm)
        kt_adm=  vpa(input('write the mode to study: \n'));
        end
        type_forcing= input('Do you want to work with Bessel function/coefficent (1) or with a unit forcing (2) ?\n');
        switch type_forcing
            case 2
                Bb_ri= 1;
            case 1
                r_investigation= input('Write the radial distance from the electrode, scaled with the height: \n');
                J1roots= vpa(besselzero(1,2,1)); %determine the roots of J1 in order to find the transverse mode
                R= (J1roots(2)-J1roots(1))/kt_adm;
                TF_b_theta_dot_kt= (1- 0.5*sigma_r*kt_adm*sqrt(vpa(pi)).*exp(-sigma_r^2*kt_adm.^2/8)...
                    .*besseli(0.5,sigma_r^2*kt_adm.^2/8)); %adimensionné par L_carac*b0
                Bb_ri= kt_adm/(vpa(pi)*(J1roots(1).*besselj(2,J1roots(1))).^2)*TF_b_theta_dot_kt; 
                ordre_B= 1;
        end

end

ordre_B= 1;
%initialisation matrice
ratio_ampl_mat= vpa(zeros(1,1));
phase_up_mat= vpa(zeros(1,1));
phase_bot_mat= vpa(zeros(1,1));
dephasage_mat= vpa(zeros(1,1));

%% Calculation of the Full MHD solutions (Hartmann and Alfven solutions)

disp("Start of calculation of the full MHD solution")


[s1,k1,s2,k2,A,B,V,J,E,Ar,dtAr]= alfven_non_homogene_vpa(Pm_num,S_num,hadm,...
    kt_adm(1),Bb_ri(1),epsilon_inf,precision,'bottom'); % 0 == bottom
if length(kt_adm) > 1
parfor i =2:length(kt_adm)-1
    [s1i,k1i,s2i,k2i,Ai,Bi,Vi,Ji,Ei]= alfven_non_homogene_vpa(vpa(Pm_num),vpa(S_num),vpa(hadm),...
        kt_adm(i),Bb_ri(i),epsilon_inf,precision,'bottom');

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



% diminution précision des paramètres
disp('diminution of digit number')
s1_less= vpa(s1,precis_b);
k1_less= vpa(k1,precis_b);
s2_less= vpa(s2,precis_b);
k2_less= vpa(k2,precis_b);
hadm_less= vpa(hadm,precis_b);

if isempty( r_investigation)
    A_less= vpa(A,precis_b);
    B_less= vpa(B,precis_b);
    V_less= vpa(V,precis_b);
    J_less= vpa(J,precis_b);    
    E_less= vpa(E,precis_b); 
    Ar_less= vpa(Ar,precis_b);
    dtAr_less= vpa(dtAr,precis_b);

elseif ~isempty( r_investigation) && length(kt_adm) >1
    % Ajout des fonctions propres aux coefficient et diminution de la
    % précisison
    A_less= vpa(A.*besselj(ordre_B,kt_adm(1:end-1)'.*r_investigation),precis_b);
    disp('fin A')
    B_less= vpa(B.*besselj(ordre_B,kt_adm(1:end-1)'.*r_investigation),precis_b);
    disp('fin B')
    V_less= vpa(V.*besselj(ordre_B,kt_adm(1:end-1)'.*r_investigation),precis_b);
    disp('fin V')
    J_less= vpa(J.*besselj(ordre_B,kt_adm(1:end-1)'.*r_investigation),precis_b);
    disp('fin J')
    E_less= vpa(E.*besselj(ordre_B,kt_adm(1:end-1)'.*r_investigation),precis_b);
    % Ar_less= vpa(Ar,precis_b);
    % dtAr_less= vpa(dtAr,precis_b);
    % Fr_phi= vpa((1- besselj(ordre_B-1,kt_adm(1:end-1)'.*r_investigation))./kt_adm(1:end-1)');
    % Fr_phi(1)= vpa(sym(0));
    % Phi_less= vpa(V.*Fr_phi,precis_b);
else
    A_less= vpa(A.*besselj(ordre_B,kt_adm.*r_investigation),precis_b);
    disp('fin A')
    B_less= vpa(B.*besselj(ordre_B,kt_adm.*r_investigation),precis_b);
    disp('fin B')
    V_less= vpa(V.*besselj(ordre_B,kt_adm.*r_investigation),precis_b);
    disp('fin V')
    J_less= vpa(J.*besselj(ordre_B,kt_adm.*r_investigation),precis_b);
    disp('fin J')
    E_less= vpa(E.*besselj(ordre_B,kt_adm.*r_investigation),precis_b);
    disp('fin E')
    %Ar_less= vpa(Ar.*besselj(ordre_B,kt_adm.*r_investigation),precis_b);
    %disp('fin Ar')
    %dtAr_less= vpa(dtAr.*besselj(ordre_B,kt_adm.*r_investigation),precis_b);
    %disp('fin dtAr')
end
disp('end diminution of digit number')

% %%
% 
% benv_mat=[];
% bphase_mat= [];
% aenv_mat= [];
% aphase_mat= [];
% eenv_mat= [];
% ephase_mat= [];   
% 
% for i= 1:5:length(kt_adm)
%     [benv_i,bphase_i,z2]=amplitudes(B_less(i,:),s1_less(i),k1_less(i),...
%                s2_less(i),k2_less(i),hadm_less,nz);
%     [aenv_i,aphase_i,z2]=amplitudes(A_less(i,:),s1_less(i),k1_less(i),...
%                s2_less(i),k2_less(i),hadm_less,nz);
%     [eenv_i,ephase_i,z2]=amplitudes(E_less(i,:),s1_less(i),k1_less(i),...
%                s2_less(i),k2_less(i),hadm_less,nz);
%     benv_mat=[benv_mat; benv_i];
%     bphase_mat= [bphase_mat; bphase_i];
%     aenv_mat= [aenv_mat ; aenv_i];
%     aphase_mat= [aphase_mat ; aphase_i];
%     eenv_mat= [eenv_mat ; eenv_i];
%     ephase_mat= [ephase_mat ; ephase_i];  
%     i
% end
% %%
% figure
% for i= 1:5:length(kt_adm)
% sgtitle(sprintf('Phase of different parameters, \n simple bottom forcing'))
% subplot(1,3,1)
% plot(bphase_mat(i,:),z)
% xlabel('magnetic field phase')
% ylabel('z axis')
% subplot(1,3,2)
% plot(aphase_mat_b(i,:),z)
% xlabel('azimuthal velocity phase')
% ylabel('z axis')
% subplot(1,3,3)
% plot(ephase_mat(i,:),z)
% xlabel('electric field')
% ylabel('z axis')
% pause(0.5)
% end

% programme TEST a supprimer si non compris
% for i= 1:nb_kt-1
% sk1=vpa(sin(k1(i)*hadm));
% ck1=vpa(cos(k1(i)*hadm));
% ep1=vpa(exp(s1(i)*hadm));
% em1=vpa(exp(-s1(i)*hadm));
% sk2=vpa(sin(k2(i)*hadm));
% ck2=vpa(cos(k2(i)*hadm));
% ep2=vpa(exp(s2(i)*hadm));
% em2=vpa(exp(-s2(i)*hadm));
% 
% N1= [ 0 1 0 1;...
%     1 0 1 0;...
%     -sk1*ep1 -ck1*ep1 sk1*em1 -ck1*em1;...
%     ck1*ep1 -sk1*ep1 ck1*em1 sk1*em1];
% N2= [ 0 1 0 1;...
%     1 0 1 0;...
%     -sk2*ep2 -ck2*ep2 sk2*em2 -ck2*em2;...
%     ck2*ep2 -sk2*ep2 ck2*em2 sk2*em2];
% 
% vect_b_s= [N1(1,:) N2(1,:)];
% vect_b_c= [N1(2,:) N2(2,:)];
% vect_t_s= [N1(3,:) N2(3,:)];
% vect_t_c= [N1(4,:) N2(4,:)];
% 
% V_b_s_somme_kt(i)= vect_b_s*V(i,:)';
% V_b_c_somme_kt(i)= vect_b_c*V(i,:)';
% V_t_s_somme_kt(i)= vect_t_s*V(i,:)';
% V_t_c_somme_kt(i)= vect_t_c*V(i,:)';
% i
% end
% %%
% figure
% plot(double(kt_adm(1:end-1)),double(V_b_s_somme_kt))
% hold on
% plot(double(kt_adm(1:end-1)),double(V_t_s_somme_kt))
% legend('bottom','top','interpreter','latex')
% xlabel('$kt_i$')
% ylabel('$Az(kt_i)$')
% title(sprintf(['Amplitudes of the sine components of the bessel-fourier coefficients for the \n'...
% ' potential gradient (freq= %4.1f Hz, S= %2.1f, $r= %4.1f$)'],freq,S_num,r_investigation))
% 
% figure
% plot(double(kt_adm(1:end-1)),double(V_b_c_somme_kt))
% hold on
% plot(double(kt_adm(1:end-1)),double(V_t_c_somme_kt))
% legend('bottom','top','interpreter','latex')
% xlabel('$kt_i$')
% ylabel('$Az(kt_i)$')
% title(sprintf(['Amplitudes of the cosine components of the bessel-fourier coefficients for the \n'...
% ' potential gradient (freq= %4.1f Hz, S= %2.1f, $r= %4.1f$)'],freq,S_num,r_investigation))
% 
% 
% 
% figure
% plot(double(kt_adm(1:end-1)),(double(V_b_s_somme_kt).^2+double(V_b_c_somme_kt).^2).^0.5)
% hold on
% plot(double(kt_adm(1:end-1)),(double(V_t_s_somme_kt).^2+double(V_t_c_somme_kt).^2).^0.5)
% legend('bottom','top','interpreter','latex')
% xlabel('$kt_i$')
% ylabel('$(Az_{sin}(kt_i)^2 + Az_{cos}(kt_i)^2)^{0.5}$')
% title(sprintf(['Quadratic amplitude of the bessel-fourier coefficients for the \n'...
% ' potential gradient (freq= %4.1f Hz, S= %2.1f, $r= %4.1f$)'],freq,S_num,r_investigation))

% Calculation of the wave enveloppe and the phases for each field
    z=0:hadm_less/(nz-1):hadm_less;
    
    % initialisation
    venv= [];
    vphase= [];
    benv=[];
    bphase= [];
    aenv= [];
    aphase= [];
    jenv= [];
    jphase= [];
    eenv= [];
    ephase= [];
    disp('Start calculation envelope')
        if length(A_less(:,1))==1
            [venv,vphase,z2]=amplitudes(V_less,s1_less,k1_less,...
               s2_less,k2_less,hadm_less,nz);
           [benv,bphase,z2]=amplitudes(B_less,s1_less,k1_less,...
               s2_less,k2_less,hadm_less,nz);
           [aenv,aphase,z2]=amplitudes(A_less,s1_less,k1_less,...
               s2_less,k2_less,hadm_less,nz);
           [jenv,jphase,z2]=amplitudes(J_less,s1_less,k1_less,...
               s2_less,k2_less,hadm_less,nz);
           [eenv,ephase,z2]=amplitudes(E_less,s1_less,k1_less,...
               s2_less,k2_less,hadm_less,nz);
%            [arenv,arphase,z2]=amplitudes(Ar_less,s1_less,k1_less,...
%                s2_less,k2_less,hadm_less,nz);

           disp('End calculation envelope')
        else
            tic
            if calculate_U == 1
                disp('Start calculation for u_\theta')
                [aenv,aphase,z2]=amplitudes_somme_onde(A_less,s1_less,k1_less,...
                    s2_less,k2_less,hadm_less,nz,precision,epsilon_inf,phase_inf);
                aphase= mk_same_congruence_phase_shift(aphase);
            end
            if calculate_V == 1
                disp('Start calculation for \partial_r\phi')
                [venv,vphase,z2]=amplitudes_somme_onde(V_less,s1_less,k1_less,...
                    s2_less,k2_less,hadm_less,nz,precision,epsilon_inf,phase_inf);
                vphase= mk_same_congruence_phase_shift(vphase);
            end
            if calculate_B == 1
                disp('Start calculation for b_\theta')
               [benv,bphase,z2]=amplitudes_somme_onde(B_less,s1_less,k1_less,...
                   s2_less,k2_less,hadm_less,nz,precision,epsilon_inf,phase_inf);
               bphase= mk_same_congruence_phase_shift(bphase);
            end
            if calculate_J == 1
                disp('Start calculation for j_r')
                [jenv,jphase,z2]=amplitudes_somme_onde(J_less,s1_less,k1_less,...
                s2_less,k2_less,hadm_less,nz,precision,epsilon_inf,phase_inf);
                jphase= mk_same_congruence_phase_shift(jphase);
            end
            if calculate_E == 1
                disp('Start calculation for E_r')
                [eenv,ephase,z2]=amplitudes_somme_onde(E_less,s1_less,k1_less,...
                    s2_less,k2_less,hadm_less,nz,precision,epsilon_inf,phase_inf);
                ephase= mk_same_congruence_phase_shift(ephase);
            end
            toc
        end
    disp('End calculation envelope')
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Codes that were used to check parts of the previous code lines          % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% %% regard mode transversal par mode transversal    
%   benv_par_mode= [];
%   jenv_par_mode= [];
% for j = 1:nb_kt
%             [benv_par_mode_j,bphase,z2]=amplitudes(B_less(j,:),s1_less(j),k1_less(j),...
%             s2_less(j),k2_less(j),hadm_less,nz);
% %            [aenv,aphase,z2]=amplitudes(A_less,s1_less,k1_less,...
% %                s2_less,k2_less,hadm_less,nz);
%             [jenv_par_mode_j,jphase,z2]=amplitudes(J_less(j,:),s1_less(j),k1_less(j),...
%                s2_less(j),k2_less(j),hadm_less,nz);
%             benv_par_mode= [benv_par_mode ; benv_par_mode_j];
%             jenv_par_mode= [jenv_par_mode ;jenv_par_mode_j];
%            j
%            
% end   
% 
% %% plot current with the analytical method
% figure
%  for j= 1:nb_kt
%     subplot(2,4,j)
%     plot(-jenv_par_mode(j,2:end-1),z(2:end-1),'r',jenv_par_mode(j,2:end-1),z(2:end-1),'r')
%     xlabel('current j')
%     ylabel('z axis')
%     legend(sprintf('kt= %4.1f',kt_adm(j)))
% end
% sgtitle(sprintf('Current envelope simple forcing (top) \n analytical calculation of rot(b)'))
% 
% 
% %% 
% figure
%  for j= 1:nb_kt
% subplot(2,4,j)
% plot(-benv_par_mode(j,2:end-1),z(2:end-1),'r',benv_par_mode(j,2:end-1),z(2:end-1),'r')
% xlabel('magnetic field')
% ylabel('z axis')
% legend(sprintf('kt= %4.1f',kt_adm(j)))
% end
% sgtitle('magnetic field envelope simple forcing (top)')
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% manual unwrap

% % velocity
% i_shift= [];
% for i= 1:length(aphase)-1
% i_shift_i= find(abs(aphase(i)-aphase(i+1)) >= pi/2);
% i_shift= [i_shift , i*i_shift_i];
% end
% aphase(i_shift(1)+1:end)= aphase(i_shift(1)+1:end) - pi;
% 
% % potential gradient
% i_shift= [];
% for i= 1:length(vphase)-1
% i_shift_i= find(abs(vphase(i)-vphase(i+1)) >= pi/2);
% i_shift= [i_shift , i*i_shift_i];
% end
% vphase(i_shift(1)+1:end)= vphase(i_shift(1)+1:end) - pi;
% 
% % electric field
% i_shift= [];
% for i= 1:length(ephase)-1
% i_shift_i= find(abs(ephase(i)-ephase(i+1)) >= pi/2);
% i_shift= [i_shift , i*i_shift_i];
% end
% ephase(i_shift(1)+1:end)= ephase(i_shift(1)+1:end) - pi;
% 
% % magnetic field
% i_shift= [];
% for i= 1:length(bphase)-1
% i_shift_i= find(abs(bphase(i)-bphase(i+1)) >= pi/2);
% i_shift= [i_shift , i*i_shift_i];
% end
% bphase(i_shift(1)+1:end)= bphase(i_shift(1)+1:end) - pi;
%%enveloppe des signaux
% figure 
% plot(-benv(1:end),z(1:end),benv(1:end),z(1:end))
% title('magnetic field envelope simple forcing (top)')
% xlabel('potential gradient')
% ylabel('z axis')

%% save
if save_data ==1 
    disp('Choose the sought folder to save data')
    selpath = uigetdir(['C:\Users\lalloz\Documents\these\onde alfven\etude theorique numerique',...
        '\1_electrode_model_results\one plate forcing case\freq_fixe']);
    save([selpath,'\donnees.mat'],'Bb_ri','aenv','aphase','benv','bphase',...
        'venv','vphase','eenv','ephase','asignal','bsignal','esignal','vsignal',...
        'time','Freq_Alfven','Ha_num','R','S_num','frequence_forcage','va_num',...
    'kt_adm','kt_max','nb_kt','r_investigation','A_less','B_less','V_less',...
    'J_less','E_less','k1_less','k2_less','s1_less','s2_less')  
end
%% figure to plot
exp_Ha= floor(log10(double(Ha_num)));
mantisse_Ha= double(Ha_num)/(10^(exp_Ha));
exp_Rnu= floor(log10(double(Rnu)));
mantisse_Rnu= double(Rnu)/(10^(exp_Rnu));
exp_r_invest= floor(log10(double(r_investigation)));
mantisse_r_invest= double(r_investigation)/(10^(exp_r_invest));

if calculate_U == 1
    Cfig= figure('name','Envelope_and_phase_of_u_\theta');
    sgtitle(sprintf(['Envelope and phase of $u_\\theta$ at $Ha =%4.2f \\times 10^%i$,\n',...
        '$R_\\eta = %4.2f$, $R_\\nu =%4.2f \\times 10^%i$, and $r = %4.2f$'],mantisse_Ha,exp_Ha,...
        double(Reta),mantisse_Rnu,exp_Rnu,double(r_investigation)))
    subplot(1,2,1)
    plot(-aenv,z,aenv,z)
    xlabel('$\max|\tilde{u}_\theta|/u_0$')
    ylabel('z')
    ax1= gca;
    subplot(1,2,2)
    plot(aphase,z)
    xlabel('$\varphi u_\theta$')
    ylabel('z')
    ax2= gca; 
    Cfig.Children(1).TickLabelInterpreter= "latex";
    Cfig.Children(2).TickLabelInterpreter= "latex";
    Cfig.Children(1).FontSize= 14;
    Cfig.Children(2).FontSize= 14;
    ax1.XAxis.FontSize= 16;
    ax1.YAxis.FontSize= 16;
    ax1.XAxis.TickLabelInterpreter="latex";
    ax1.YAxis.TickLabelInterpreter="latex";
    ax2.XAxis.FontSize= 16;
    ax2.YAxis.FontSize= 16;
    ax2.XAxis.TickLabelInterpreter="latex";
    ax2.YAxis.TickLabelInterpreter="latex";
    %set

end

if calculate_V == 1
    Cfig= figure('name','Envelope_and_phase_of_grad_phi_r');
    sgtitle(sprintf(['Envelope and phase of $\\partial_r \\phi$ at $Ha =%4.2f \\times 10^%i$,\n',...
        '$R_\\eta = %4.2f$, $R_\\nu =%4.2f \\times 10^%i$, and $r = %4.2f$'],mantisse_Ha,exp_Ha,...
        double(Reta),mantisse_Rnu,exp_Rnu,double(r_investigation)))
    subplot(1,2,1)
    plot(-venv,z,venv,z)
    xlabel('$\max|\partial_r \tilde\phi|/E_0$')
    ylabel('z')
    subplot(1,2,2)
    ax1= gca;
    plot(vphase,z)
    xlabel('$\varphi_\Delta \partial_r \phi$')
    ylabel('z')
    ax2= gca; 
    Cfig.Children(1).TickLabelInterpreter= "latex";
    Cfig.Children(2).TickLabelInterpreter= "latex";
    Cfig.Children(1).FontSize= 14;
    Cfig.Children(2).FontSize= 14;
    ax1.XAxis.FontSize= 16;
    ax1.YAxis.FontSize= 16;
    ax1.XAxis.TickLabelInterpreter="latex";
    ax1.YAxis.TickLabelInterpreter="latex";
    ax2.XAxis.FontSize= 16;
    ax2.YAxis.FontSize= 16;
    ax2.XAxis.TickLabelInterpreter="latex";
    ax2.YAxis.TickLabelInterpreter="latex";
end

if calculate_B == 1
    Cfig= figure('name','Envelope_and_phase_of_b_\theta');
    sgtitle(sprintf(['Envelope and phase of $b_\\theta$ at $Ha =%4.2f \\times 10^%i$,\n',...
        '$R_\\eta = %4.2f$, $R_\\nu =%4.2f \\times 10^%i$, and $r = %4.2f$'],mantisse_Ha,exp_Ha,...
        double(Reta),mantisse_Rnu,exp_Rnu,double(r_investigation)))
    subplot(1,2,1)
    plot(-benv,z,benv,z)
    xlabel('$\max|\tilde{b}_\theta|/b_0$')
    ylabel('z')
    ax1= gca;
    subplot(1,2,2)
    plot(bphase,z)
    xlabel('$\varphi_\Delta b_\theta$')
    ylabel('z')
    ax2= gca; 
    Cfig.Children(1).TickLabelInterpreter= "latex";
    Cfig.Children(2).TickLabelInterpreter= "latex";
    Cfig.Children(1).FontSize= 14;
    Cfig.Children(2).FontSize= 14;
    ax1.XAxis.FontSize= 16;
    ax1.YAxis.FontSize= 16;
    ax1.XAxis.TickLabelInterpreter="latex";
    ax1.YAxis.TickLabelInterpreter="latex";
    ax2.XAxis.FontSize= 16;
    ax2.YAxis.FontSize= 16;
    ax2.XAxis.TickLabelInterpreter="latex";
    ax2.YAxis.TickLabelInterpreter="latex";
end

if calculate_E == 1
    Cfig= figure('name','Envelope_and_phase_of_E_r');
    sgtitle(sprintf(['Envelope and phase of $E_r$ at $Ha =%4.2f \\times 10^%i$,\n',...
        '$R_\\eta = %4.2f$, $R_\\nu =%4.2f \\times 10^%i$, and $r = %4.2f$'],mantisse_Ha,exp_Ha,...
        double(Reta),mantisse_Rnu,exp_Rnu,double(r_investigation)))
    subplot(1,2,1)
    plot(-eenv,z,eenv,z)
    xlabel('$\max|\tilde{E}_r|/E_0$')
    ylabel('z')
    ax1= gca;
    subplot(1,2,2)
    plot(ephase,z)
    xlabel('$\varphi_\Delta E_r$')
    ylabel('z')
    ax2= gca; 
    Cfig.Children(1).TickLabelInterpreter= "latex";
    Cfig.Children(2).TickLabelInterpreter= "latex";
    Cfig.Children(1).FontSize= 14;
    Cfig.Children(2).FontSize= 14;
    ax1.XAxis.FontSize= 16;
    ax1.YAxis.FontSize= 16;
    ax1.XAxis.TickLabelInterpreter="latex";
    ax1.YAxis.TickLabelInterpreter="latex";
    ax2.XAxis.FontSize= 16;
    ax2.YAxis.FontSize= 16;
    ax2.XAxis.TickLabelInterpreter="latex";
    ax2.YAxis.TickLabelInterpreter="latex";
end

if calculate_J == 1
    Cfig= figure('name','Envelope_and_phase_of_j_r');
    sgtitle(sprintf(['Envelope and phase of $j_r$ at $Ha =%4.2f \\times 10^%i$,\n',...
        '$R_\\eta = %4.2f$, $R_\\nu =%4.2f \\times 10^%i$, and $r = %4.2f$'],mantisse_Ha,exp_Ha,...
        double(Reta),mantisse_Rnu,exp_Rnu,double(r_investigation)))
    subplot(1,2,1)
    plot(-jenv,z,jenv,z)
    xlabel('$\max|\tilde{j}_r|/j_0$')
    ylabel('z')
    ax1= gca;
    subplot(1,2,2)
    plot(jphase,z)
    xlabel('$\varphi_\Delta j_r$')
    ylabel('z')
    ax2= gca; 
    Cfig.Children(1).TickLabelInterpreter= "latex";
    Cfig.Children(2).TickLabelInterpreter= "latex";
    Cfig.Children(1).FontSize= 14;
    Cfig.Children(2).FontSize= 14;
    ax1.XAxis.FontSize= 16;
    ax1.YAxis.FontSize= 16;
    ax1.XAxis.TickLabelInterpreter="latex";
    ax1.YAxis.TickLabelInterpreter="latex";
    ax2.XAxis.FontSize= 16;
    ax2.YAxis.FontSize= 16;
    ax2.XAxis.TickLabelInterpreter="latex";
    ax2.YAxis.TickLabelInterpreter="latex";
end

%% calcul des signaux aux limites (top and bottom)
if calculate_boundary_signal == 1
    t=[0:0.01*2*pi/epsilon_inf:2*pi/(epsilon_inf)];
    lim= 0.5*length(A(:,1));
   
    if calculate_U == 1
    % champ vitesse
        [asignal_top,asignal_bottom]=calcul_signal_boundaries(A_less,s1_less,k1_less,...
               s2_less,k2_less,epsilon_inf, hadm_less,t,precision,sym('0'));
           disp('fin u')
    end
    if calculate_V == 1
        % gradient de potentiel
        [vsignal_top,vsignal_bottom]=calcul_signal_boundaries(V_less,s1_less,k1_less,...
                   s2_less,k2_less,epsilon_inf, hadm_less,t,precision,sym('0'));
        disp('fin grad_r phi')
    end
    if calculate_B == 1
        % champ magnétique
        [bsignal_top,bsignal_bottom]=calcul_signal_boundaries(B_less,s1_less,k1_less,...
                   s2_less,k2_less,epsilon_inf, hadm_less,t,precision,sym('0'));
           disp('fin b')
    end
    if calculate_J == 1
    %     % courant électrique
        [jsignal_top,jsignal_bottom]=calcul_signal_boundaries(J_less,s1_less,k1_less,...
                   s2_less,k2_less,epsilon_inf, hadm_less,t,precision,sym('0'));
         disp('fin jr')
    end
    if calculate_E == 1
    %     % champ électrique loi d'ohm
        [esignal_top,esignal_bottom]=calcul_signal_boundaries(E_less,s1_less,k1_less,...
                   s2_less,k2_less,epsilon_inf, hadm_less,t,precision,sym('0'));
           disp('fin Er')
    end
end

%% tracé des signaux aux limites
if calculate_boundary_signal == 1
    if calculate_U == 1
        Cfig= figure('name','boundary_signal_u_\theta');
        plot(t,asignal_top)
        hold on
        plot(t,asignal_bottom)
        xlabel('t')
        ylabel('$u_\theta (r,z,t)$')
        legend('$z=1$','$z=0$')
        title(sprintf(['Boundary signals of $u_\\theta$ at $Ha =%4.2f \\times 10^%i$,\n',...
            '$R_\\eta = %4.2f$, $R_\\nu =%4.2f \\times 10^%i$, and $r =%4.2f \\times 10^{%i}$'],mantisse_Ha,exp_Ha,...
            double(Reta),mantisse_Rnu,exp_Rnu,mantisse_r_invest,exp_r_invest))
        ax= gca;
        Cfig.Children(2).TickLabelInterpreter= "latex";
        Cfig.Children(2).FontSize= 12;
        Cfig.Children(1).FontSize= 12;
        Cfig.Children(1).Interpreter= "latex";
        ax.XAxis.FontSize= 16;
        ax.YAxis.FontSize= 16;
        ax.XAxis.TickLabelInterpreter="latex";
        ax.YAxis.TickLabelInterpreter="latex";
    end
    if calculate_B == 1
        Cfig= figure('name','boundary_signal_b_theta');
        plot(t,bsignal_top)
        hold on
        plot(t,bsignal_bottom)
        xlabel('t')
        ylabel('$b_\theta (r,z,t)$')
        legend('$z=1$','$z=0$')
        title(sprintf(['Boundary signals of $b_\\theta$ at $Ha =%4.2f \\times 10^%i$,\n',...
            '$R_\\eta = %4.2f$, $R_\\nu =%4.2f \\times 10^%i$, and $r =%4.2f \\times 10^{%i}$'],mantisse_Ha,exp_Ha,...
            double(Reta),mantisse_Rnu,exp_Rnu,mantisse_r_invest,exp_r_invest))
        ax= gca;
        Cfig.Children(2).TickLabelInterpreter= "latex";
        Cfig.Children(2).FontSize= 12;
        Cfig.Children(1).FontSize= 12;
        Cfig.Children(1).Interpreter= "latex";
        ax.XAxis.FontSize= 16;
        ax.YAxis.FontSize= 16;
        ax.XAxis.TickLabelInterpreter="latex";
        ax.YAxis.TickLabelInterpreter="latex";
    end
    if calculate_V == 1
        Cfig= figure('name','boundary_signal_grad_r phi');
        plot(t,vsignal_top)
        hold on
        plot(t,vsignal_bottom)
        xlabel('t')
        ylabel('$\partial_r \phi(r,z,t)$')
        legend('$z=1$','$z=0$')
        title(sprintf(['Boundary signals of $\\partial_r \\phi$ at $Ha =%4.2f \\times 10^%i$,\n',...
            '$R_\\eta = %4.2f$, $R_\\nu =%4.2f \\times 10^%i$, and $r =%4.2f \\times 10^{%i}$'],mantisse_Ha,exp_Ha,...
            double(Reta),mantisse_Rnu,exp_Rnu,mantisse_r_invest,exp_r_invest))
        ax= gca;
        Cfig.Children(2).TickLabelInterpreter= "latex";
        Cfig.Children(2).FontSize= 12;
        Cfig.Children(1).FontSize= 12;
        Cfig.Children(1).Interpreter= "latex";
        ax.XAxis.FontSize= 16;
        ax.YAxis.FontSize= 16;
        ax.XAxis.TickLabelInterpreter="latex";
        ax.YAxis.TickLabelInterpreter="latex";
    end
    if calculate_E == 1
    
        Cfig= figure('name','boundary_signal_E_r');
        plot(t,esignal_top)
        hold on
        plot(t,esignal_bottom)
        xlabel('t')
        ylabel('$E_r(r,z,t)$')
        legend('$z=1$','$z=0$')
        title(sprintf(['Boundary signals of $E_r$ at $Ha =%4.2f \\times 10^%i$,\n',...
            '$R_\\eta = %4.2f$, $R_\\nu =%4.2f \\times 10^%i$, and $r =%4.2f \\times 10^{%i}$'],mantisse_Ha,exp_Ha,...
            double(Reta),mantisse_Rnu,exp_Rnu,mantisse_r_invest,exp_r_invest))
        ax= gca;
        Cfig.Children(2).TickLabelInterpreter= "latex";
        Cfig.Children(2).FontSize= 12;
        Cfig.Children(1).FontSize= 12;
        Cfig.Children(1).Interpreter= "latex";
        ax.XAxis.FontSize= 16;
        ax.YAxis.FontSize= 16;
        ax.XAxis.TickLabelInterpreter="latex";
        ax.YAxis.TickLabelInterpreter="latex";
    end
    if calculate_J == 1
        Cfig= figure('name','boundary_signal_j_r');
        plot(t,jsignal_top)
        hold on
        plot(t,jsignal_bottom)
        xlabel('t')
        ylabel('$j(r,z,t)$')
        legend('$z=1$','$z=0$')
        title(sprintf(['Boundary signal of $j_r$ at $Ha =%4.2f \\times 10^%i$,\n',...
            '$R_\\eta = %4.2f$, $R_\\nu =%4.2f \\times 10^%i$, and $r =%4.2f \\times 10^{%i}$'],mantisse_Ha,exp_Ha,...
            double(Reta),mantisse_Rnu,exp_Rnu,mantisse_r_invest,exp_r_invest))
        ax= gca;
        Cfig.Children(2).TickLabelInterpreter= "latex";
        Cfig.Children(2).FontSize= 12;
        Cfig.Children(1).FontSize= 12;
        Cfig.Children(1).Interpreter= "latex";
        ax.XAxis.FontSize= 16;
        ax.YAxis.FontSize= 16;
        ax.XAxis.TickLabelInterpreter="latex";
        ax.YAxis.TickLabelInterpreter="latex";
    end
end
 

%%  calculation time evolution along z
nb_clted_variable= calculate_U + calculate_B + calculate_J + calculate_E + calculate_V;
if calculate_time_evolution_along_z == 1
    time=[0:0.2*2*vpa(pi)/epsilon_inf:2*vpa(pi)/(epsilon_inf)];
    length(time)
    
    i=1
    length(time)
    asignal= [];
    vsignal= [];
    bsignal= [];
    jsignal= [];
    esignal= [];
    phisignal= [];
       
    for t = time
        
        if calculate_U == 1
            % champ vitesse
            asignal_inf= calcul_signal_along_z(A_less,s1_less,k1_less,...
                   s2_less,k2_less,epsilon_inf,z,t,precision,0);
            asignal= [asignal ; asignal_inf ];           
            disp('fin u')
        end
        if calculate_V == 1
            % gradient de potentiel
            vsignal_inf= calcul_signal_along_z(V_less,s1_less,k1_less,...
                    s2_less,k2_less,epsilon_inf,z,t,precision,0);
            vsignal= [vsignal ; vsignal_inf ];    
            disp('fin grad_r phi')
        end
        if calculate_B == 1
            % champ magnétique
            bsignal_inf= calcul_signal_along_z(B_less,s1_less,k1_less,...
                   s2_less,k2_less,epsilon_inf,z,t,precision,0);
            bsignal= [bsignal ; bsignal_inf ];
           disp('fin b')
        end
        if calculate_J == 1
            % courant électrique
            jsignal_inf= calcul_signal_along_z(J_less,s1_less,k1_less,...
                    s2_less,k2_less,epsilon_inf,z,t,precision,0);
            jsignal= [jsignal ; jsignal_inf ];
            disp('fin j')
        end
        if calculate_E == 1
            % champ électrique loi d'ohm
            esignal_inf= calcul_signal_along_z(E_less,s1_less,k1_less,...
                   s2_less,k2_less,epsilon_inf,z,t,precision,0);  
            esignal= [esignal ; esignal_inf ];
           disp('fin Er')
        end
    
   
    i= i+1
    end
end
%% figure which show the time evolution
%
% 
set(0,'defaultTextInterpreter','latex')
opengl software
fig_evolution_suivant_z= figure;
if exist('video_struc')
    clear video_struc
end
for i = 1:length(time)
    k= 1;
    if calculate_U == 1
        s_ax= subplot(1,nb_clted_variable,k);
        plot(aenv,z,'r',-aenv,z,'r')
        hold on
        plot(asignal(i,:),z,'k');
        hold off
        title(sprintf('$u_\\theta(r,z,t)$'))
        xlabel('$u_\theta$')
        ylabel('$z$')
        s_ax.FontSize= 10;
        s_ax.TickLabelInterpreter= "latex";
        s_ax.XLabel.FontSize= 10;
        s_ax.YLabel.FontSize= 10;
        s_ax.Title.FontSize= 12;
        k= k+1;
    end
    if calculate_B == 1
        s_ax= subplot(1,nb_clted_variable,k);
        plot(benv,z,'r',-benv,z,'r')
        hold on
        plot(bsignal(i,:),z,'k');
        hold off
        title(sprintf('$b_\\theta (r,z,t)$'))
        xlabel('$b_\theta$')
        ylabel('$z$')
        s_ax.FontSize= 10;
        s_ax.TickLabelInterpreter= "latex";
        s_ax.XLabel.FontSize= 12;
        s_ax.YLabel.FontSize= 12;
        s_ax.Title.FontSize= 12;
        k= k+1;
    end
    if calculate_E == 1
        s_ax= subplot(1,nb_clted_variable,k);
        plot(eenv,z,'r',-eenv,z,'r')
        hold on
        plot(esignal(i,:),z,'k');
        hold off
        title(sprintf('$E_r(r,z,t)$'))
        xlabel('$E_r$')
        ylabel('$z$')
        s_ax.FontSize= 10;
        s_ax.TickLabelInterpreter= "latex";
        s_ax.XLabel.FontSize= 12;
        s_ax.YLabel.FontSize= 12;
        s_ax.Title.FontSize= 12;
        k= k+1;
    end
    if calculate_V == 1
        s_ax= subplot(1,nb_clted_variable,k);
        plot(venv,z,'r',-venv,z,'r')
        hold on
        plot(vsignal(i,:),z,'k');
        hold off
        title(sprintf('$\\partial_r \\phi(r,z,t)$'))
        xlabel('$\partial_r \phi$')
        ylabel('$z$')
        s_ax.FontSize= 10;
        s_ax.TickLabelInterpreter= "latex";
        s_ax.XLabel.FontSize= 12;
        s_ax.YLabel.FontSize= 12;
        s_ax.Title.FontSize= 12;
        k= k+1;
    end
    if calculate_J == 1      
        s_ax= subplot(1,nb_clted_variable,k);
        plot(jenv,z,'r',-jenv,z,'r')
        hold on
        plot(jsignal(i,:),z,'k');
        hold off
        title(sprintf('$j_r(r,z,t)$'))
        xlabel('$j_r$')
        ylabel('$z$')
        hold off
        s_ax.FontSize= 10;
        s_ax.TickLabelInterpreter= "latex";
        s_ax.XLabel.FontSize= 12;
        s_ax.YLabel.FontSize= 12;
        s_ax.Title.FontSize= 12;
        k= k+1;
    end

sgtitle(sprintf(['Enveloppe of different parameters, at $Ha =%4.2f \\times 10^%i$,\n',...
        '$R_\\eta = %4.2f$, $R_\\nu =%4.2f \\times 10^%i$, and $r = %4.2f \\times 10^{%i}$'],...
        mantisse_Ha,exp_Ha, double(Reta),mantisse_Rnu,exp_Rnu,mantisse_r_invest,exp_r_invest))
pause(0.5)
i
video_struc(i)= getframe(fig_evolution_suivant_z);
hold off
end

%% sauvegarde video
if save_video_time_evol_along_z == 1
    video_struc_bis= struct;
    disp('Choose the sought folder to save figure: ')
    selpath = uigetdir(pwd);
    for i=1:length(video_struc)
        video_struc_bis(i)= video_struc(i);
        i
    end
    v= VideoWriter('film_evolution_signal_suivant_z.avi');
    v.FrameRate= 1;
    open(v)
    writeVideo(v,video_struc_bis)
    close(v)
    clear v

end

