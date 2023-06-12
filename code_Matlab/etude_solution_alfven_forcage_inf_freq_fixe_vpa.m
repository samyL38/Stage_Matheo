clear all
clc
%close all


% Calcul des amplitudes pour le potentiel, vitesse, champ magnétique induit
% au niveau des plaques inférieures et supérieures. Le programme utilise la
% fonction alfven_non_homogene qui permet de résoudre la relation de
% dispersion full MHD pour un cas non_homogène (existance d'un nombre
% d'onde transverse k_t)
% On considère ici que seul le courant injecté produit de la vorticité (la
% plaque inférieure n'est pas oscillante)

% Le forcage est produit parles deux plaques a deux frequences données
% différentes. Néanmoins la distribution spatiale du forcage électrique est
% le même sur la plaque supérieure que sur la plaque inférieure

%L_carac est la hauteur du récipient

% Adimensionnement de :
% - j par I0/L_carac^2
% - b par b0= mu_0*I0/(L_carac)
% - TF_b_theta_kt par L_carac^2*b0
% - u par u0= eta/L_carac* b0/B0
% - tau= L_carac^2/eta;
%%%%

%% toolbox flottant haute précision
%addpath('C:\Users\lalloz\Documents\Multiprecision Computing Toolbox')
precision= 4*32;% % Setup default precision to 40 decimal digits (quadruple).
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
%Paramètres physiques
mu_0= vpa(1.25*sym(1e-6));
sigma_gal= vpa(3.46*sym(1e6));
visco_gal= vpa(0.0024); %Viscosité dynamique du Galstn en Pa.s
eta= 1/(sigma_gal * mu_0); %diffusivité mhd
rho_gal= vpa(6440);% masse volumique Galinstan kg/m3
nu_gal= visco_gal/rho_gal;

% Paramètres géométriques/expérimentaux
B0= vpa(5); %champ magnetique uniforme en Tesla
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
va_num= B0/sqrt(rho_gal*mu_0); %Vitesse d'Alfven 
S_num= va_num*L_carac/eta; %Nombre de Lundquist basé sur la hauteur récipient
Ha_num= B0*L_carac*sqrt(sigma_gal/visco_gal); %Nombre d'Hartmann basé sur la hauteur récipient

Tau_Alfven= L_carac/va_num;
Tau_Ha= h^2/(2*nu_gal*Ha_num); 
Tau_joule= rho_gal/(sigma_gal*B0^2);
Tau_2D= Tau_joule*(h/L_trans)^2;
Freq_Alfven= 1/Tau_Alfven;

%%%% Controle parameter: frequency %%%%
frequence_forcage= Freq_Alfven;% ; %en Hz %0.01*tau^-1
% frequence
if ~exist('frequence_forcage') || isempty(frequence_forcage)
    frequence_forcage= input('Write the frequency (in Hz) to test : \n');
end
omega_inf= vpa(2*sym(pi)*frequence_forcage);%390);
epsilon_inf= tau_eta*omega_inf;
Pm_num= vpa(nu_gal/(eta)); %nb de Prandtl magnétique 

%%% New control parameter
Rnu_fixed= 0;
Reta=  tau_eta*2*vpa(pi)*frequence_forcage;
Rnu= tau_nu*2*vpa(pi)*frequence_forcage;

if Rnu_fixed == 1
    Pm_num= Reta/Rnu;
    S_num= Ha_num*sqrt(Reta/Rnu);
    epsilon_inf= Reta;
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
nb_point_z= 25;
nb_kt= 1200;%2500;%3200;%1000;%2100;%
delta_kt= 5*2.25;%0.5;%
r_investigation= vpa(0.005);
kt_adm= [];


cas_single_mode= input('Study of a single transverse mode case ? Yes (1), No (0) \n');

switch cas_single_mode
    
    case 0
        if ~exist('r_investigation') || isempty(r_investigation)
        r_investigation= input('Write the scaled location for the calculation: \n');
        end
        if ~exist('nb_kt') || isempty(nb_kt)
        nb_kt= input('Write the number of transverse mode to test: \n');
        end
        if ~exist('delta_kt') || isempty(delta_kt)
        delta_kt= input('Value of the discretisation for kt ?\n'); %1 exemple, 355   
        end
        ordre_B= 1;
        % location and transverse mode vectors
        J1roots= vpa(besselzero(1,nb_kt,1)); %determine the roots of J1 in order to find the transverse mode
        R= (J1roots(2)-J1roots(1))/vpa(delta_kt); % taille du domaine de définition du forçage magnétique
        kt_adm= J1roots/R; %mode transversaux
        fprintf('Value of the domain R %4.2f \n',double(R))

        % Calcul des coefficients Bb_ri, projections du forcage magnetique suivant
        % les modes propres transversaux
        TF_b_theta_dot_kt= double((1- 0.5*kt_adm*sigma_r*sqrt(vpa(pi)).*exp(-sigma_r^2*kt_adm.^2/8)...
           .*besseli(0.5,sigma_r^2*kt_adm.^2/8))); %adimensionné par L_carac*b0
       Bb_ri=double(kt_adm)/(pi*(double(J1roots.*besselj(2,J1roots))).^2).*TF_b_theta_dot_kt; 
       Bb_ri= vpa(Bb_ri);
    case 1
        if ~exist('kt_adm') || isempty(kt_adm)
        kt_adm=  vpa(input('write the mode to study: \n'));
        end
        type_forcing= input('Do you want to work with Bessel function/coefficent (1) or with a unit forcing (2) ?\n');
        switch type_forcing
            case 2
                Bb_ri= 1;
            case 1
                r_investigation= input('Write the value of the scaled location: \n');
                J1roots= vpa(besselzero(1,2,1)); %determine the roots of J1 in order to find the transverse mode
                R= (J1roots(2)-J1roots(1))/kt_adm;
                %TF_b_theta_dot_kt= double((1- 0.5*kt_adm*sigma_r*sqrt(vpa('pi')).*exp(-sigma_r^2*kt_adm.^2/8)...
                %    .*besseli(0.5,sigma_r^2*kt_adm.^2/8))); %adimensionné par L_carac*b0
                %Bb_ri=double(kt_adm)/(pi*(double(J1roots(1).*besselj(2,J1roots(1)))).^2).*TF_b_theta_dot_kt; 
                %Bb_ri= vpa(Bb_ri);
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

% solution for the bottom plate
disp('Calculating alfven wave solutions for the bottom plate forcing...')
%DigitsOld= digits(precision);
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
disp('end solution searching for bottom plate forcing')



% diminution précision des paramètres
disp('diminution de la précision des paramètre')
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
%    disp('fin A')
%    B_less= vpa(B.*besselj(ordre_B,kt_adm(1:end-1)'.*r_investigation),precis_b);
    disp('fin B')
    V_less= vpa(V.*besselj(ordre_B,kt_adm(1:end-1)'.*r_investigation),precis_b);
    disp('fin V')
%     J_less= vpa(J.*besselj(ordre_B,kt_adm(1:end-1)'.*r_investigation),precis_b);
%     disp('fin J')
%    E_less= vpa(E.*besselj(ordre_B,kt_adm(1:end-1)'.*r_investigation),precis_b);
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
    Ar_less= vpa(Ar.*besselj(ordre_B,kt_adm.*r_investigation),precis_b);
    disp('fin Ar')
    dtAr_less= vpa(dtAr.*besselj(ordre_B,kt_adm.*r_investigation),precis_b);
    disp('fin dtAr')
end
disp('fin diminution précision')

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
    disp('debut calcul enveloppe')
        if length(A_less(:,1))==1
            [venv,vphase,z2]=amplitudes(V_less,s1_less,k1_less,...
               s2_less,k2_less,hadm_less,nz);
           [benv,bphase,z2]=amplitudes(B_less,s1_less,k1_less,...
               s2_less,k2_less,hadm_less,nz);
           [aenv,aphase,z2]=amplitudes(A_less,s1_less,k1_less,...
               s2_less,k2_less,hadm_less,nz);
           [dtarenv,dtarphase,z2]=amplitudes(dtAr_less,s1_less,k1_less,...
               s2_less,k2_less,hadm_less,nz);
           [jenv,jphase,z2]=amplitudes(J_less,s1_less,k1_less,...
               s2_less,k2_less,hadm_less,nz);
           [eenv,ephase,z2]=amplitudes(E_less,s1_less,k1_less,...
               s2_less,k2_less,hadm_less,nz);
           [arenv,arphase,z2]=amplitudes(Ar_less,s1_less,k1_less,...
               s2_less,k2_less,hadm_less,nz);

           disp('fin')
        else
            tic
            disp('debut calcul pour vitesse')
%             [aenv,aphase,z2]=amplitudes_somme_onde(A_less,s1_less,k1_less,...
%                 s2_less,k2_less,hadm_less,nz,precision,epsilon_inf,phase_inf);
            disp('debut calcul pour gradient potentiel')
            [venv,vphase,z2]=amplitudes_somme_onde(V_less,s1_less,k1_less,...
                s2_less,k2_less,hadm_less,nz,precision,epsilon_inf,phase_inf);
            disp('debut calcul pour champ magnetic')
%            [benv,bphase,z2]=amplitudes_somme_onde(B_less,s1_less,k1_less,...
%                s2_less,k2_less,hadm_less,nz,precision,epsilon_inf,phase_inf);
%            disp('debut calcul pour courant')
%             [jenv,jphase,z2]=amplitudes_somme_onde(J_less,s1_less,k1_less,...
%                 s2_less,k2_less,hadm_less,nz,precision,epsilon_inf,phase_inf);
%              disp('debut calcul champ électrique')
%             [eenv,ephase,z2]=amplitudes_somme_onde(E_less,s1_less,k1_less,...
%                 s2_less,k2_less,hadm_less,nz,precision,epsilon_inf,phase_inf);
%             disp('debut calcul potentiel electrique')
%             [phienv,phiphase,z2]=amplitudes_somme_onde(Phi_less,s1_less,k1_less,...
%                 s2_less,k2_less,hadm_less,nz,precision,epsilon_inf,phase_inf);

            toc
        end
    disp('fin calcul enveloppe')
  
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

figure
sgtitle(sprintf(['Enveloppe of different parameters with a simple bottom forcing, \n'...
'non zero Diriclet conditions (freq= %4.1f Hz, S= %2.1f, $r= %4.1f$)'],frequence_forcage,S_num,r_investigation))
% subplot(1,3,1)
% plot(-benv(1:end),z(1:end),benv(1:end),z(1:end))
% title('')
% xlabel('magnetic field')
% ylabel('z axis')

% subplot(1,3,2)
% plot(-aenv(1:end),z(1:end),aenv(1:end),z(1:end))
% title('')
% xlabel('azimuthal velocity')
% ylabel('z axis')

subplot(1,3,2)
plot(-venv(1:end),z(1:end),venv(1:end),z(1:end))
title('')
xlabel('potential gradient ')
ylabel('z axis')

subplot(1,3,3)
% plot(-eenv(1:end),z(1:end),eenv(1:end),z(1:end))
% title('')
% xlabel('electric field')
% ylabel('z axis')

% subplot(1,5,4)
% plot(-jenv(1:end),z(1:end),jenv(1:end),z(1:end))
% title('')
% xlabel('radial current')
% ylabel('z axis')
% 

figure
sgtitle(sprintf(['Phase of different parameters with a simple bottom forcing, \n',...
'Neumann top condition (freq= %4.1f Hz, S= %2.1f, $r= %4.1f$)'],frequence_forcage,S_num,r_investigation))
% subplot(1,3,1)
% plot(bphase(1:end),z(1:end))
% title('')
% xlabel('magnetic field phase')
% ylabel('z axis')

subplot(1,3,2)
% plot(aphase(1:end),z(1:end))
% title('')
% xlabel('azimuthal velocity phase')
% ylabel('z axis')

subplot(1,3,3)
plot(vphase(1:end),z(1:end))
title('')
xlabel('potential gradient phase ')
ylabel('z axis')

% subplot(1,3,3)
% plot(ephase(1:end),z(1:end))
% title('')
% xlabel('electric field')
% ylabel('z axis')

% subplot(1,5,4)
% plot(-jenv(1:end),z(1:end),jenv(1:end),z(1:end))
% title('')
% xlabel('radial current')
% ylabel('z axis')
% 



%% calcul des signaux aux limites (top and bottom)
t=[0:0.01*2*pi/epsilon_inf:2*pi/(epsilon_inf)];
lim= 0.5*length(A(:,1));
tic
    % champ vitesse
%     [asignal_top,asignal_bottom]=calcul_signal_boundaries(A_less,s1_less,k1_less,...
%            s2_less,k2_less,epsilon_inf, hadm_less,t,precision,sym('0'));
%        disp('fin a')
    % gradient de potentiel
    [vsignal_top,vsignal_bottom]=calcul_signal_boundaries(V_less,s1_less,k1_less,...
               s2_less,k2_less,epsilon_inf, hadm_less,t,precision,sym('0'));
       disp('fin v')
    % champ magnétique
%     [bsignal_top,bsignal_bottom]=calcul_signal_boundaries(B_less,s1_less,k1_less,...
%                s2_less,k2_less,epsilon_inf, hadm_less,t,precision,sym('0'));
%        disp('fin b')
%     % courant électrique
%     [jsignal_top,jsignal_bottom]=calcul_signal_boundaries(J_less,s1_less,k1_less,...
%                s2_less,k2_less,epsilon_inf, hadm_less,t,precision,sym('0'));
%     % champ électrique loi d'ohm
    [esignal_top,esignal_bottom]=calcul_signal_boundaries(E_less,s1_less,k1_less,...
               s2_less,k2_less,epsilon_inf, hadm_less,t,precision,sym('0'));
       disp('fin e')
% toc

%% tracé des signaux aux limites

% figure 
% plot(t,jsignal_top)
% hold on
% plot(t,jsignal_bottom)
% xlabel('scaled time')
% ylabel('current')
% legend('signal Jtop','signal Jbottom')

figure 
plot(t,asignal_top)
hold on
plot(t,asignal_bottom)
xlabel('scaled time')
ylabel('current')
legend('signal velocity top','signal velocity down')

figure 
plot(t,vsignal_top)
hold on
plot(t,vsignal_bottom)
xlabel('scaled time')
ylabel('potential gradient')
legend('signal grad(phi) top','signal grad(phi) bottom')

figure 
plot(t,esignal_top)
hold on
plot(t,esignal_bottom)
xlabel('scaled time')
ylabel('electric field')
legend('signal Er top','signal Er bottom')

figure 
plot(t,bsignal_top)
hold on
plot(t,bsignal_bottom)
xlabel('scaled time')
ylabel('mangetic field')
legend('signal b top','signal b bottom')
%%
figure 
plot(t,phisignal_top)
hold on
plot(t,phisignal_bottom)
xlabel('non-scaled time')
ylabel('electric potential')
legend('signal phi top','signal phi bottom')
%%  
time=[0:0.2*2*vpa(pi)/epsilon_inf:2*vpa(pi)/(epsilon_inf)];
length(time)
%%
i=1
asignal= [];
vsignal= [];
bsignal= [];
jsignal= [];
esignal= [];
phisignal= [];
for t = time
asignal_inf= calcul_signal_along_z(A_less,s1_less,k1_less,...
               s2_less,k2_less,epsilon_inf,z,t,precision,0);
% vsignal_inf= calcul_signal_along_z(V_less,s1_less,k1_less,...
%                s2_less,k2_less,epsilon_inf,z,t,precision,0);
            
bsignal_inf= calcul_signal_along_z(B_less,s1_less,k1_less,...
               s2_less,k2_less,epsilon_inf,z,t,precision,0);

% jsignal_inf= calcul_signal_along_z(J_less,s1_less,k1_less,...
%                s2_less,k2_less,epsilon_inf,z,t,precision,0);
           
esignal_inf= calcul_signal_along_z(E_less,s1_less,k1_less,...
               s2_less,k2_less,epsilon_inf,z,t,precision,0);  
           
%phisignal_inf= calcul_signal_along_z(Phi_less,s1_less,k1_less,...
%               s2_less,k2_less,epsilon_inf,z,t,precision,0);                
% j_ohm_signal_inf= calcul_signal_along_z(J_ohm_less,s1_less,k1_less,...
%                s2_less,k2_less,epsilon_inf,z,t,precision,0);

asignal= [asignal ; asignal_inf ];           
%vsignal= [vsignal ; vsignal_inf ];
bsignal= [bsignal ; bsignal_inf ];
% jsignal= [jsignal ; jsignal_inf ];
esignal= [esignal ; esignal_inf ];
%phisignal= [phisignal ; phisignal_inf ];

i= i+1
end

%% save
selpath = uigetdir(['C:\Users\lalloz\Documents\these\onde alfven\etude theorique numerique',...
    '\1_electrode_model_results\one plate forcing case\freq_fixe']);
save([selpath,'\donnees.mat'],'Bb_ri','aenv','aphase','benv','bphase',...
    'venv','vphase','eenv','ephase','asignal','bsignal','esignal','vsignal',...
    'time','Freq_Alfven','Ha_num','R','S_num','frequence_forcage','va_num',...
'kt_adm','kt_max','nb_kt','r_investigation','A_less','B_less','V_less',...
'J_less','E_less','k1_less','k2_less','s1_less','s2_less')  

%%

set(0,'defaultTextInterpreter','latex')
opengl software
fig_evolution_suivant_z= figure;



for i = [1:1:length(time)]
% subplot(1,5,1)
% plot(venv,z,'r',-venv,z,'r')
% hold on
% plot(vsignal(i,:),z,'k');
% hold off
% title('Electric potential')
% xlabel('$V/B_0\eta$')
% ylabel('$z (\omega/\eta)^{1/2}$')

subplot(1,3,2)
plot(aenv,z,'r',-aenv,z,'r')
hold on
plot(asignal(i,:),z,'k');
hold off
title(sprintf('velocity field'))
xlabel('$u/u_0$')
ylabel('$z/L_{carac}$')

subplot(1,3,3)
plot(eenv,z,'r',-eenv,z,'r')
hold on
plot(esignal(i,:),z,'k');
hold off
title('electric field')
xlabel('$E/ E_0$')
ylabel('$z/L_{carac}$')

% subplot(1,2,2)
% plot(venv,z,'r',-venv,z,'r')
% hold on
% plot(vsignal(i,:),z,'k');
% title('potential gradient')
% xlabel('$u/ u_0$')
% ylabel('$z/L_{carac}$')
% hold off

subplot(1,3,1)
plot(benv,z,'r',-benv,z,'r')
hold on
plot(bsignal(i,:),z,'k');
hold off
title('magnetic field')
xlabel('$b/ b_0$')
ylabel('$z/L_{carac}$')
% % 
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
% subplot(1,4,4)
% plot(eenv,z,'r',-eenv,z,'r')
% hold on
% plot(esignal(i,:),z,'k');
% hold off
% title('electric potental')
% xlabel('$E/ E_0$')
% ylabel('$z/L_{carac}$')




sgtitle(sprintf('Enveloppe of different parameters, \n simple bottom forcing'))
pause(0.5)
video_struc(i)= getframe(fig_evolution_suivant_z);
i
hold off



end
%% sauvegarde video

    for i=1:length(video_struc)
        video_struc_bis(i)= video_struc(i);
        i
    end
    v= VideoWriter('film_evolution_signal_suivant_z_B_7T_Finf_5Hz.avi');
    v.FrameRate= 1;
    open(v)
    writeVideo(v,video_struc_bis)
    close(v)
    clear v



          
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




