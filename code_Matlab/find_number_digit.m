%%% Determination number of digits for calculating Full MHD solutions

%% toolbox flottant haute précision
%addpath('C:\Users\lalloz\Documents\Multiprecision Computing Toolbox')


set(0,'defaultTextInterpreter','latex')
opengl software
precis_b= 32;


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


%% Paramètre du problème
syms mu_0 sigma_gal visco_gal rho_gal omega_0 B0 h I0 kt_max epsilon

% Parameters to define
B0_double= []; %champ magnetique uniforme en Tesla
h_double= []; % distance entre la plaque inf et sup en mètre
frequence_forcage_double= [];% ; %en Hz %0.01*tau^-1
nb_kt= [];%1200;%2500;%3200;%1000;%2100;%
r_investigation_double= [];
R_double= [];

if ~exist('B0_double') || isempty(B0_double)
    B0_double= vpa(input('Write the magnetic field intensity in Tesla : \n'));
end
if ~exist('h_double') || isempty(h_double)
    h_double= vpa(input('Write the height of the vessel in meter : \n'));
end
if ~exist('frequence_forcage_double') || isempty(frequence_forcage_double)
    frequence_forcage_double= input('Write the frequency (in Hz) to test : \n');
end
if ~exist('r_investigation_double') || isempty(r_investigation_double)
r_investigation_double= input('Write the scaled location for the calculation: \n');
end
if ~exist('nb_kt') || isempty(nb_kt)
nb_kt= input('Write the number of transverse mode to test: \n');
end
if ~exist('R_double') || isempty(R_double)
    R_double= input('Write the radial size of the vessel (scaled with h) ?\n'); % exemple, 355   
end
ordre_B= 1;

digit_multiple= 2;
precision= digit_multiple*32;% % Setup default precision to 40 decimal digits (quadruple).
error_born= 10^8;
digits(precision);  
find_brk=0;

while find_brk == 0
    precision= digit_multiple*32
    digits(precision);
    %Paramètres physiques
    mu_0= vpa(1.25*sym(1e-6));
    sigma_gal= vpa(3.46*sym(1e6));
    visco_gal= vpa(0.0024); %Viscosité dynamique du Galstn en Pa.s
    eta= 1/(sigma_gal * mu_0); %diffusivité mhd
    rho_gal= vpa(6440);% masse volumique Galinstan kg/m3
    nu_gal= visco_gal/rho_gal;

    % Paramètres géométriques/expérimentaux
    B0= vpa(B0_double); %champ magnetique uniforme en Tesla
    diam_elec= 1e-3; % en mètre
    rayon_elec= diam_elec/2;
    dist_inter_elec= 2*10^-3; % en mètre
    h= vpa(h_double); % distance entre la plaque inf et sup en mètre
    phase_inf= sym(0);

    % Grandeurs caractéristiques
    L_carac=h;%(eta/omega_0)^(1/2);%longueur caractéristique ==> hauteur de la cuve
    tau_eta= L_carac^2/eta;
    tau_nu= L_carac^2/nu_gal;

    % Paramètres adimensionnées
    sigma_r= rayon_elec/L_carac;% diam_adm de l'électrode
    hadm=h/L_carac;
    L_trans= rayon_elec;

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
    frequence_forcage= vpa(frequence_forcage_double);% ; %en Hz %0.01*tau^-1

    omega_inf= vpa(2*sym(pi)*frequence_forcage);%390);
    epsilon_inf= tau_eta*omega_inf;
    Pm_num= vpa(nu_gal/(eta)); %nb de Prandtl magnétique 

    %%% New control parameter
    Reta=  tau_eta*2*vpa(pi)*frequence_forcage;
    Rnu= tau_nu*2*vpa(pi)*frequence_forcage;



    % Définition of the boundary magnetic perturbation
    %(avec solution de Hartmann et Alfven)

    R= vpa(R_double);
    r_investigation= vpa(r_investigation_double);

    % location and transverse mode vectors
    J1roots= vpa(besselzero(1,nb_kt,1)); %determine the roots of J1 in order to find the transverse mode
    kt_adm= J1roots/R; %mode transversaux
    % Calcul des coefficients Bb_ri, projections du forcage magnetique suivant
    % les modes propres transversaux
    TF_b_theta_dot_kt= double((1- 0.5*kt_adm*sigma_r*sqrt(vpa(pi)).*exp(-sigma_r^2*kt_adm.^2/8)...
       .*besseli(0.5,sigma_r^2*kt_adm.^2/8))); %adimensionné par L_carac*b0
    Bb_ri=double(kt_adm)/(pi*(double(J1roots.*besselj(2,J1roots))).^2).*TF_b_theta_dot_kt; 
    Bb_ri= vpa(Bb_ri);

    ordre_B= 1;



    % Calculation of the Full MHD solutions (Hartmann and Alfven solutions)

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
        i;

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

    V_less= vpa(V.*besselj(ordre_B,kt_adm(1:end-1)'.*r_investigation),precis_b);
    disp('fin V')

    Produit1= abs(V_less(:,5).*exp(s2_less));
    Produit2= abs(V_less(:,2).*exp(s1_less));

    nb_mode_err_1= sum( double(Produit1) > error_born);
    nb_mode_err_2= sum( double(Produit1) > error_born);
    nb_mode_err= max(nb_mode_err_1, nb_mode_err_2);

    if nb_mode_err >= 5
        digit_multiple= digit_multiple + ceil(nb_mode_err/5);
        fprintf('number of wrong modes: %i \n',nb_mode_err);
        
    elseif nb_mode_err >=1
        digit_multiple= digit_multiple +1;
        fprintf('number of wrong modes: %i \n',nb_mode_err);
    else 
        find_brk= 1;
        digit_final= digit_multiple*32;
        fprintf('The number of digits needed for given condition is %i \n',digit_final);
    end

end
