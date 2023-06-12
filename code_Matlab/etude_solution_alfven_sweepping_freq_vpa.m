clear all
clc
%close all


% Calcul des amplitudes pour le potentiel, vitesse, champ magnétique induit
% au niveau des plaques inférieures et supérieures. Le programme utilise la
% fonction alfven_non_homogene_vpa qui permet de résoudre la relation de
% dispersion full MHD pour un cas non_homogène (existance d'un nombre
% d'onde transverse k_t)
% On considère ici que seul le courant injecté produit de la vorticité (la
% plaque inférieure n'est pas oscillante)

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
precision=12*32;% Setup default precision to 40 decimal digits (quadruple).
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

%For r_investigation < 0.2
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
diam_elec= 1e-3; % en mètre
rayon_elec= diam_elec/2;
dist_inter_elec= 20*10^-3; % en mètre
h= vpa(0.1); % distance entre la plaque inf et sup en mètre

% Grandeurs caractéristiques
L_carac=h;%(eta/omega_0)^(1/2);%longueur caractéristique ==> hauteur de la cuve
tau_eta= L_carac^2/eta;
tau_nu= L_carac^2/nu_gal;
I0=2; % Injection courant en Ampère
j0= I0/(L_carac)^2; %%% forme de jz = j0/(pi*sigma_r^2)*exp(-(r/sigma_r)^2)
b0= mu_0*I0/(L_carac); %%% forme de b_theta= b0*(1-exp(-((r)/sigma_r).^2))./(r);
E0= eta/L_carac * b0;
joule_time= rho_gal/(sigma_gal*B0^2);
alfven_time=h*sqrt(rho_gal*mu_0)/B0;
% 2D time
TwoD_time= joule_time*10.^2;
TwoD_time_norm= TwoD_time*0.5;

% 2D time
Lpara= h;
Ltrans= 0.1*h;
TwoD_time= joule_time*(Lpara/Ltrans)^2;
TwoD_velocity= h/TwoD_time;


% Paramètres adimensionnées
sigma_r= rayon_elec/L_carac;% diam_adm de l'électrode
hadm=h/L_carac;



% Paramètre tracé
Rmax= hadm;%r_exp_adm;
nr= 5*10^2;%floor(Rmax);
nz= 100;




% Nombre adm
Pm_num= vpa(10^3);%nu_gal/(eta); %nb de Prandtl magnétique 
va_num= B0/sqrt(rho_gal*mu_0); %Vitesse d'Alfven 
S_num=va_num*L_carac/eta; %Nombre de Lundquist basé sur la hauteur récipient
Ha_num= B0*L_carac*sqrt(sigma_gal/visco_gal); %Nombre d'Hartmann basé sur la hauteur récipient

Freq_Alfven= B0/sqrt(rho_gal*mu_0)/L_carac;

%% Champ à calculer
% Entrer ici les champs à calculer
calculate_U= 0;
calculate_V= [];
calculate_B= 0;
calculate_E= [];
calculate_J= [];
save_data= [];
if ~exist('calculate_U') || isempty(calculate_U)
    calculate_U=input('Do you want to calculate the velcoity field ? Yes (1), No (0) \n');
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
    calculate_J=input('Do you want to calculate the current density ? Yes (1), No (0) \n');
end
if ~exist('save_data') || isempty(save_data)
    save_data=input(...
        'Do you want to save data ? Yes (1), No (0) \n');
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

% parameters by default
nb_point_z= 25;
nb_kt= 100;%1200;%2500;%3200;%1000;%2100;%
delta_kt= 2.25;%0.5;%
r_investigation= vpa(0.005);
kt_adm= [];
freq_min= [];
freq_min= [];
nb_omega= [];

if ~exist('freq_min') || isempty(freq_min)
        freq_min= input('Write the minimum frequency (in Hz) to test : \n');
end
if ~exist('freq_max') || isempty(freq_max)
        freq_max= input('Write the maximal frequency (in Hz) to test : \n');
end
if ~exist('nb_omega') || isempty(nb_omega)
        nb_omega= input('Write the number of frequencies to test: \n');
end


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
        R= (J1roots(2)-J1roots(1))/delta_kt; % taille du domaine de définition du forçage magnétique
        kt_adm= J1roots/R; %mode transversaux
        fprintf('Value of the domain R %4.2f \n',double(R))

        % Calcul des coefficients Bb_ri, projections du forcage magnetique suivant
        % les modes propres transversaux
        TF_b_theta_dot_kt= ((1- 0.5*kt_adm*sigma_r*sqrt(vpa(pi)).*exp(-sigma_r^2*kt_adm.^2/8)...
           .*besseli(0.5,sigma_r^2*kt_adm.^2/8))); %adimensionné par L_carac*b0
       Bb_ri=(kt_adm)./(pi*(J1roots.*besselj(2,J1roots)).^2).*TF_b_theta_dot_kt; 
       %Bb_ri= 2./(R*besselj(2,J1roots)).^2 /(2*pi).*TF_b_theta_dot_kt;

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
                ordre_B= 1;1;
        end
end


% plage frequence
omega_min= vpa(2*vpa(pi)*freq_min);
omega_max= vpa(2*vpa(pi)*freq_max);
omega_vect= linspace(omega_min,omega_max,nb_omega);
freq_vect= 1/(2*vpa(pi))* omega_vect;
epsilon_vect= tau*omega_vect;

%%% New control parameter
Rnu_fixed= 0;
Reta_vect=  tau_eta*2*vpa(pi)*freq_vect;
Rnu_vect= tau_nu*2*vpa(pi)*freq_vect;




%% Calcul des solutions de la relation de dispersion full MHD 
%(avec solution de Hartmann et Alfven)

%initialisation matrice
ratio_ampl_mat= vpa(zeros(nb_omega,1));
phase_up_mat= vpa(zeros(nb_omega,1));
phase_bot_mat= vpa(zeros(nb_omega,1));
dephasage_mat= vpa(zeros(nb_omega,1));

% frequency sweeping ampl/phase shift calculation

% initialisation vecteur signal max top/bottom
time= [];
gradV_env= [];
gradV_phase=[];
E_env= [];
E_phase=[];
J_env= [];
J_phase= [];



for j= 1:nb_omega

    epsilon= epsilon_vect(j); 
    disp('Calculating alfven wave solutions for the bottom plate forcing...')
    DigitsOld= digits(precision);
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

if isempty(r_investigation)

    % diminution précision des paramètres
    disp('diminution de la précision des paramètre')
    s1_less= vpa(s1,precis_b);
    k1_less= vpa(k1,precis_b);
    s2_less= vpa(s2,precis_b);
    k2_less= vpa(k2,precis_b);
    hadm_less= vpa(hadm,precis_b);
    A_less= vpa(A,precis_b);
    B_less= vpa(B,precis_b);
    V_less= vpa(V,precis_b);
    J_less= vpa(J,precis_b);    
    E_less= vpa(E,precis_b);    
    
elseif   ~ isempty( r_investigation) 
    



    % Ajout des fonctions propres aux coefficient et diminution de la
    % précisison
        s1_less= vpa(s1,precis_b);
        k1_less= vpa(k1,precis_b);
        s2_less= vpa(s2,precis_b);
        k2_less= vpa(k2,precis_b);
        hadm_less= vpa(hadm,precis_b);
    %     A_less= vpa(A.*besselj(ordre_B,kt_adm(1:end-1)'*r_investigation),precis_b);
    %     B_less= vpa(B.*besselj(ordre_B,kt_adm(1:end-1)'*r_investigation),precis_b);
    if calculate_V == 1
         V_less= vpa(V.*besselj(ordre_B,kt_adm(1:end-1)'.*r_investigation),precis_b);
    end
    if calculate_J == 1
         J_less= vpa(J.*besselj(ordre_B,kt_adm(1:end-1)'*r_investigation),precis_b);
    end
    if calculate_E == 1
         E_less= vpa(E.*besselj(ordre_B,kt_adm(1:end-1)'.*r_investigation),precis_b);
    end
    %     Fr_phi= vpa((1- besselj(ordre_B-1,kt_adm(1:end-1)'.*r_investigation))./kt_adm(1:end-1)');
    %     Fr_phi(1)= vpa(sym(0));
    %     Phi_less= vpa(V.*Fr_phi,precis_b);
    %  
end
   disp('fin diminution précision')    
   j  
   fprintf('progressing %4.2f%% ..\n',j*100/nb_omega)
 
nz= 2;
if calculate_V == 1
    disp('debut calcul enveloppe gradient phi')
    [gradV_env_j,gradV_phase_j,z2]= amplitudes_somme_onde(V_less,s1_less,k1_less,...
             s2_less,k2_less,hadm_less,nz,precision,epsilon,sym(0)); 
    gradV_env= [gradV_env ; gradV_env_j];
    gradV_phase= [gradV_phase ; gradV_phase_j];

    disp('fin calcul enveloppe grad potentiel')
end
if calculate_E
    disp('debut calcul enveloppe electric field')
    [E_env_j,E_phase_j,z2]= amplitudes_somme_onde(E_less,s1_less,k1_less,...
            s2_less,k2_less,hadm_less,nz,precision,epsilon,sym(0));    
    E_env= [E_env ; E_env_j];
    E_phase= [E_phase ; E_phase_j];
    disp('fin calcul enveloppe electric field')
end
if calculate_J
    disp('debut calcul enveloppe radial cuurent')
    [J_env_j,J_phase_j,z2]= amplitudes_somme_onde(J_less,s1_less,k1_less,...
            s2_less,k2_less,hadm_less,nz,precision,epsilon,sym(0));   
    J_env= [J_env ; J_env_j];
    J_phase= [J_phase ; J_phase_j];
    disp('fin calcul enveloppe electric field')
end

disp('fin calcul enveloppe')

                    
end

%% frequence (in the case if omega and epsilon have not been saved

mantisse_Ha= Ha_num/10^(floor(log10(Ha_num)));
exp_Ha= floor(log10(Ha_num));
choice_freq= 'freq_adm_Fav'; % write either 'Reta' or 'freq', freq_adm_Fav
switch choice_freq
    case 'Reta'
        xvect= double(epsilon_vect)/(2*pi);
        label_x= '$R_{\eta}/2\pi$';
        label_tle= '$R_{\eta}/2\pi$';
    case 'freq'
        xvect= freq_vect;
        label_x= '$f_{inj}$ (Hz)';
        label_tle= '$f_{inf}$';

    case 'freq_adm_Fav'
    xvect= freq_vect/double(Freq_Alfven);
    label_x= '$\sqrt{R_\nu R_\eta}/Ha$ (Hz)';
    label_tle= '$f_{inf}/Fav$';
end

label_gradV_ampl= '$\\sup\\left\\{\\nabla_r\\phi(x= %4.2f ,y= %4.2f, z)\\right\\}/E_0$';
label_gradV_phase= '$\\varphi_{\\nabla_r\\phi}(z)$';
label_gradV_tle_grad= '$\partial_x\phi$';


label_E_ampl= '$\\sup\\left\\{E_r(x= %4.2f ,y= %4.2f, z)\\right\\}/E_0$';
label_E_phase= '$\\varphi_{E_r}(z)$';
label_E_tle_grad= '$E_r$';

label_j_ampl= '$\\sup\\left\\{j_r(x= %4.2f ,y= %4.2f, z)\\right\\}/E_0$';
label_j_phase= '$\\varphi_{j_r}(z)$';
label_j_tle_grad= '$j_r$';


%% calculation attenuation coefficient and phase shift at the boundaries
i_lim= length(freq_vect);
if calculate_E
    E_max_top= E_env(:,end);  
    E_max_bot= E_env(:,1); 
    E_att_coef= log(E_max_top/E_max_bot); 
    E_phase_top= unwrap(E_phase(1:i_lim,end));  %unwrap
    E_phase_bot= (E_phase(1:i_lim,1)); %unwrap
    E_phase_top=mk_same_congruence_phase_shift(E_phase_top);   
    E_dephasage= E_phase_top - E_phase_bot;
end

if calculate_V
    gradV_max_top= gradV_env(:,end); 
    gradV_max_bot= gradV_env(:,1);
    gradV_att_coef= log(E_max_top/gradV_max_bot); 
    gradV_phase_top= unwrap(gradV_phase(:,end)); %unwrap
    gradV_phase_bot= gradV_phase(:,1); 
    gradV_phase_top=mk_same_congruence_phase_shift(gradV_phase_top);
    gradV_dephasage= gradV_phase_top - gradV_phase_bot;
end

if calculate_J
    J_max_top= J_env(:,end); 
    J_max_bot= J_env(:,1);
    J_att_coef= log(J_max_top/J_max_bot); 
    J_phase_top= unwrap(J_phase(:,end)); %unwrap
    J_phase_bot= J_phase(:,1); 
    J_phase_top=mk_same_congruence_phase_shift(J_phase_top);
    J_dephasage= J_phase_top - J_phase_bot;
end



% %% transit time and velocity
% gradV_time_delay= -double(gradV_dephasage./(2*vpa(pi)*double(freq_vect)'));
% E_time_delay= -double(E_dephasage./(2*vpa(pi)*double(freq_vect)'));
% gradV_velocity= double(h./gradV_time_delay);
% E_velocity= double(h./E_time_delay);
% 
% gradV_velocity_adm= double(gradV_velocity/va_num);
% E_velocity_adm= double(E_velocity/va_num);



%% save
if save_data == 1
    disp('choose the folder to save data')
    selpath = uigetdir(['C:\Users\lalloz\Documents\these\onde alfven\etude theorique numerique',...
        '\1_electrode_model_results\one plate forcing case\freq_sweeping']);
    save([selpath,'\donnees.mat'],'Bb_ri','gradV_env','gradV_phase','Freq_Alfven',...
        'Ha_num','R','S_num','freq_vect','freq_max','freq_min','va_num',...
    'E_env','E_phase','kt_adm','kt_max','delta_kt','nb_kt','r_investigation','sigma_r')  

end
    
%% Figure bot & top gradV amplitude vs frequency
if calculate_V == 1
    figure
    plot(xvect,gradV_max_bot)
    hold on
    plot(xvect,gradV_max_top)
    ylabel(sprintf('$|\\partial\\phi(r,z)|$ '))
    xlabel(label_x)
    ax= gca;
    ax.TickLabelInterpreter="latex";
    grid on 
    l_fig= legend('$z = 0$','$z = 1$','Interpreter','Latex');
    title(sprintf(['Potentiel gradient amplitude versus %s \n'...
        ' ($Ha= %4.1f \\times 10^%i, r= %4.2f$)'],label_tle,mantisse_Ha,exp_Ha,double(r_investigation))) 
    set(gca, 'YScale', 'log')
    f_ampl_ax= gca;
    f_ampl_ax.TickLabelInterpreter= 'latex';
    f_ampl_ax.FontSize= 10;
    f_ampl_ax.YLabel.FontSize = 14;
    f_ampl_ax.XLabel.FontSize = 14;
    l_fig.FontSize= 12;
    
    % Figure spatial attenuation coefficient
    figure('name','attenuation coef GradV versus freq')
    %yyaxis left
    plot(xvect,gradV_att_coef)
    hold on
    xlabel(label_x)
    ylabel(sprintf('$\\alpha\\left(r = %4.1f, z= 1\\right)$',r_investigation))
    l_fig= legend(sprintf('$Ha= %4.2e$',Ha_num),'interpreter','latex');
    title(sprintf(['Spatial attenuation coefficient of %s versus %s\n'...
    ' ($Ha= %4.1f \\times 10^%i$, $r= %4.2f$)'],label_gradV_tle_grad,label_tle,...
    mantisse_Ha,exp_Ha,double(r_investigation)))
    %txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
    %t= text(700,0.4,txt);
    %t.Units= 'normalized';
    f_att_coef_ax= gca;
    f_att_coef_ax.TickLabelInterpreter= 'latex';
    f_att_coef_ax.FontSize= 10;
    f_att_coef_ax.YLabel.FontSize = 14;
    f_att_coef_ax.XLabel.FontSize = 14;
    l_fig.FontSize= 12;
    
    grid on
    
    % Figure phase shift versus Freq
    f_phase_shift_ax= figure('name','phase shift versus freq');
    %yyaxis left
    plot(xvect,gradV_dephasage)
    hold on
    xlabel(label_x)
    ylabel(sprintf('$\\varphi_{\Delta}\\left(r = %4.1f, z= 1\\right)$',r_investigation))
    l_fig= legend(sprintf('$Ha= %4.2e$',Ha_num),'interpreter','latex');
    title(sprintf(['Phase shift of %s versus %s\n'...
    ' ($Ha= %4.1f \\times 10^%i, r= %4.2f$)'],label_gradV_tle_grad,label_tle,...
    mantisse_Ha,exp_Ha,double(r_investigation))) 
    f_phase_shift_ax.TickLabelInterpreter= 'latex';
    f_phase_shift_ax.FontSize= 12;
    f_phase_shift_ax.YLabel.FontSize = 14;
    f_phase_shift_ax.XLabel.FontSize = 14;
    l_fig.FontSize= 12;
    grid on

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% A finir pour les autres paramètres %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Figure transit time scaled with the alfven time coefficient
figure
%yyaxis left
plot(freq_vect*TwoD_time,gradV_time_delay*Freq_Alfven)%/Freq_Alfven
hold on
xlabel('$\mathit{\Omega}/\left(2\pi \right)$')
ylabel('$\tau_{\nabla\phi}/\tau_A$')
title(sprintf(['scaled transit time for the potential gradient\n'...
' ($S= %4.1f, r= %4.2f$)'],double(S_num),double(r_investigation)))
%txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
%t= text(700,0.4,txt);
%t.Units= 'normalized';
f_transit_time_ax= gca;
f_transit_time_ax.TickLabelInterpreter= 'latex';
f_transit_time_ax.FontSize= 10;
f_transit_time_ax.YLabel.FontSize = 15;
f_transit_time_ax.XLabel.FontSize = 12;
grid on


