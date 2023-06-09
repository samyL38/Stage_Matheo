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
    path_com_function= dir('common_functions');
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
tau= L_carac^2/eta;
I0=2; % Injection courant en Ampère
j0= I0/(L_carac)^2; %%% forme de jz = j0/(pi*sigma_r^2)*exp(-(r/sigma_r)^2)
b0= mu_0*I0/(L_carac); %%% forme de b_theta= b0*(1-exp(-((r)/sigma_r).^2))./(r);
E0= eta/L_carac * b0;
joule_time= rho_gal/(sigma_gal*B0^2);
alfven_time=h*sqrt(rho_gal*mu_0)/B0;

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

%% Calcul des solutions de la relation de dispersion full MHD 
%(avec solution de Hartmann et Alfven)

% entrée des paramètres d'intérêt : 
% r_min : distance minimale en mètre depuis le centre de l'électrode
% r_max : distance maximale en mètre depuis le centre de l'électrode
% nb_point_r= nombre de points de calcul
% nb_kt= nombre de mode propre transversaux à utiliser pour l'approximation
%       du forcage magnetique (nb_kt est adimensionné par L_carac= 0.1m)
% delta_kt= pas entre deux modes transversaux
r_investigation= [];

freq_min= input('Write the minimum frequency (in Hz) to test : \n');
freq_max= input('Write the maximal frequency (in Hz) to test : \n');
nb_omega= input('Write the number of frequencies to test: \n');
cas_single_mode= input('Study of a single transverse mode case ? Yes (1), No (0) \n');

switch cas_single_mode
    
    case 0
        r_investigation= input('Write the scaled location for the frequency sweeping: \n');
        nb_kt= input('Write the number of transverse mode to test: \n');
        delta_kt= input('Value of the discretisation for kt ?\n'); %1 exemple, 355   
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
        kt_adm=  vpa(input('write the mode to study: \n'));
        
        type_forcing= input('Do you want to work with Bessel function/coefficent (1) or with a unit forcing (2) ?\n');
        switch type_forcing
            case 2
                Bb_ri= 1;
            case 1
                r_investigation= input('Write the value of the scaled location: \n');
                J1roots= vpa(besselzero(1,3,1)); %determine the roots of J1 in order to find the transverse mode
                R= (J1roots(2)-J1roots(1))/0.5;
                TF_b_theta_dot_kt= (1- 0.5*sigma_r*kt_adm*sqrt(vpa(pi)).*exp(-sigma_r^2*kt_adm.^2/8)...
                    .*besseli(0.5,sigma_r^2*kt_adm.^2/8)); %adimensionné par L_carac*b0
                Bb_ri= kt_adm./(vpa(pi)*(J1roots.*besselj(2,J1roots)).^2).*TF_b_theta_dot_kt; 
                %Bb_ri= 2./(R*besselj(2,kt_adm*R)).^2 /(2*pi).*TF_b_theta_dot_kt;
                ordre_B= 1;
        end
        R= Inf;
        nb_kt= 1;

end

% plage frequence
omega_min= vpa(2*vpa(pi)*freq_min);
omega_max= vpa(2*vpa(pi)*freq_max);
omega_mat= linspace(omega_min,omega_max,nb_omega);
freq_mat= 1/(2*vpa(pi))* omega_mat;
epsilon_mat= tau*omega_mat;





%% Calcul des solutions de la relation de dispersion full MHD 
%(avec solution de Hartmann et Alfven)

%initialisation matrice
ratio_ampl_mat= vpa(zeros(nb_omega,1));
phase_up_mat= vpa(zeros(nb_omega,1));
phase_bot_mat= vpa(zeros(nb_omega,1));
dephasage_mat= vpa(zeros(nb_omega,1));

% frequency sweeping ampl/phase shift calculation

% initialisation vecteur signal max top/bottom
E_top= [];
E_bot= [];
gradV_top= [];
gradV_bot= [];
time= [];
B_top= [];
B_bot= [];

gradV_env= [];
gradV_phase=[];
E_env= [];
E_phase=[];



for j= 1:nb_omega


    epsilon= epsilon_mat(j); 
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
         V_less= vpa(V.*besselj(ordre_B,kt_adm(1:end-1)'.*r_investigation),precis_b);
    %     J_less= vpa(J.*besselj(ordre_B,kt_adm(1:end-1)'*r_investigation),precis_b);
         E_less= vpa(E.*besselj(ordre_B,kt_adm(1:end-1)'.*r_investigation),precis_b);
    %     Fr_phi= vpa((1- besselj(ordre_B-1,kt_adm(1:end-1)'.*r_investigation))./kt_adm(1:end-1)');
    %     Fr_phi(1)= vpa(sym(0));
    %     Phi_less= vpa(V.*Fr_phi,precis_b);
    %  
end
   disp('fin diminution précision')    
   j  
   fprintf('progressing %4.2f%% ..\n',j*100/nb_omega)
  % % Calculation of the wave enveloppe and the phases for each field
%     z=0:hadm_less/(nz-1):hadm_less;
%     disp('debut calcul enveloppe gradient phi')
%         if length(A_less(:,1))==1
%             [benv,bphase,z2]=amplitudes(B_less,s1_less,k1_less,...
%                s2_less,k2_less,hadm_less,nz,precision);
%         else
%             [benv,bphase,z2]=amplitudes_somme_onde(B_less,s1_less,k1_less,...
%                 s2_less,k2_less,hadm_less,nz,precision);
%             disp('fin enveloppe grad potentiel')
%         end
%     disp('fin calcul enveloppe')
%     
% % Calculation of the wave enveloppe and the phases for each field
% z=0:hadm_less/(nz-1):hadm_less;
%     disp('debut calcul enveloppe gradient phi')
%     if length(A_less(:,1))==1
%         [venv_j,vphase_j,z2]=amplitudes(V_less,s1_less,k1_less,...
%            s2_less,k2_less,hadm_less,2,precision);
%     else
%     [venv_j,vphase_j,z2]=amplitudes_somme_onde(V_less,s1_less,k1_less,...
%         s2_less,k2_less,hadm_less,2,precision,epsilon,vpa(0));
%     end
%     venv= [venv ; venv_j];
%     vphase= [vphase ; vphase_j];
%     disp('fin enveloppe grad potentiel')
% %     end
% disp('fin calcul enveloppe')

%%% wave enveloppe

% % Calculation of the wave enveloppe and the phases for each field
%     z=0:hadm_less/(nz-1):hadm_less;
nz= 2;
disp('debut calcul enveloppe gradient phi')
[gradV_env_j,gradV_phase_j,z2]= amplitudes_somme_onde(V_less,s1_less,k1_less,...
         s2_less,k2_less,hadm_less,nz,precision,epsilon,sym(0));     

disp('fin calcul enveloppe grad potentiel')

disp('debut calcul enveloppe electric field')
[E_env_j,E_phase_j,z2]= amplitudes_somme_onde(E_less,s1_less,k1_less,...
        s2_less,k2_less,hadm_less,nz,precision,epsilon,sym(0));      
disp('fin calcul enveloppe electric field')

disp('fin calcul enveloppe')

gradV_env= [gradV_env ; gradV_env_j];
gradV_phase= [gradV_phase ; gradV_phase_j];
E_env= [E_env ; E_env_j];
E_phase= [E_phase ; E_phase_j];
%     


% calcul des signaux aux limites (top and bottom)
% ti=[0:0.01*2*pi/epsilon:2*pi/(epsilon)];
% lim= 0.5*length(A(:,1));
% %     % champ vitesse
% %     [asignal_top,asignal_bottom]=calcul_signal_boundaries(A_less,s1_less,k1_less,...
% %            s2_less,k2_less,epsilon_inf, hadm_less,t,precision);
%     % gradient de potentiel
%     [vsignal_top_j,vsignal_bottom_j]=calcul_signal_boundaries(V_less,s1_less,k1_less,...
%                s2_less,k2_less,epsilon, hadm_less,ti,precision);
% %         % champ magnétique
% %         [bsignal_top_j,bsignal_bottom_j]=calcul_signal_boundaries(B_less,s1_less,k1_less,...
% %                    s2_less,k2_less,epsilon, hadm_less,ti,precision);
% %     % courant électrique
% %     [jsignal_top,jsignal_bottom]=calcul_signal_boundaries(J_less,s1_less,k1_less,...
% %                s2_less,k2_less,epsilon_inf, hadm_less,t,precision);
%     % champ électrique loi d'ohm
%     [esignal_top_j,esignal_bottom_j]=calcul_signal_boundaries(E_less,s1_less,k1_less,...
%                s2_less,k2_less,epsilon, hadm_less,ti,precision);
% % 
% disp('fin calcul des signaux aux limites')
%     gradV_top_j= vsignal_top_j;
%     gradV_bot_j= vsignal_bottom_j;
%     E_top_j= esignal_top_j;
%     E_bot_j= esignal_bottom_j;
%     B_top_j= bsignal_top_j;
%     B_bot_j= bsignal_bottom_j;           
%     
%     gradV_top= [gradV_top ; gradV_top_j];
%     gradV_bot= [gradV_bot ; gradV_bot_j];
%     time= [time ; ti];
%     E_top= [E_top ; E_top_j];
%     E_bot= [E_bot ; E_bot_j];                          
%     B_top= [B_top ; B_top_j];
%     B_bot= [B_bot ; B_bot_j];           
%                    
end

%% frequence (in the case if omega and epsilon have not been saved
omega_mat= freq_mat*2*vpa(pi);
epsilon_mat=tau*omega_mat;


mantisse_Ha= Ha_num/10^(floor(log10(Ha_num)));
exp_Ha= floor(log10(Ha_num));
choice_freq= 'freq_adm_Fav'; % write either 'Reta' or 'freq', freq_adm_Fav
switch choice_freq
    case 'Reta'
        xvect= double(epsilon_mat)/(2*pi);
        label_x= '$R_{\eta}/2\pi$';
        label_tle= '$R_{\eta}/2\pi$';
    case 'freq'
        xvect= freq_mat;
        label_x= '$f_{inj}$ (Hz)';
        label_tle= '$f_{inf}$';

    case 'freq_adm_Fav'
    xvect= freq_mat/double(Freq_Alfven);
    label_x= '$\sqrt{R_\nu R_\eta}/Ha$ (Hz)';
    label_tle= '$f_{inf}/Fav$';
end

label_gradV_ampl= '$\\sup\\left\\{\\nabla_r\\phi(x= %4.2f ,y= %4.2f, z)\\right\\}/E_0$';
label_gradV_phase= '$\\varphi_{\\nabla_r\\phi}(z)$';
label_gradV_tle_grad= '$\nabla_x\phi$';


label_E_ampl= '$\\sup\\left\\{E_r(x= %4.2f ,y= %4.2f, z)\\right\\}/E_0$';
label_E_phase= '$\\varphi_{E_r}(z)$';
label_E_tle_grad= '$E_r$';


%% amplitude ratio
i_lim= length(freq_mat);
gradV_max_top= gradV_env(1:i_lim,end); 
gradV_max_bot= gradV_env(1:i_lim,1);
E_max_top= E_env(1:i_lim,end);  
E_max_bot= E_env(1:i_lim,1); 

gradV_phase_top= unwrap(gradV_phase(1:i_lim,end)); %unwrap
gradV_phase_bot= gradV_phase(1:i_lim,1); 
E_phase_top= unwrap(E_phase(1:i_lim,end));  %unwrap
E_phase_bot= (E_phase(1:i_lim,1)); %unwrap

gradV_phase_top=mk_same_congruence_phase_shift(gradV_phase_top);
E_phase_top=mk_same_congruence_phase_shift(E_phase_top);

gradV_ampl_ratio= double(gradV_max_top./gradV_max_bot); 
E_ampl_ratio= double(E_max_top./E_max_bot);

gradV_att_coef= log(gradV_ampl_ratio); 
E_att_coef= log(E_ampl_ratio);

E_dephasage= E_phase_top - E_phase_bot;
gradV_dephasage= gradV_phase_top - gradV_phase_bot;

%% transit time and velocity
gradV_time_delay= -double(gradV_dephasage./(2*vpa(pi)*double(freq_mat)'));
E_time_delay= -double(E_dephasage./(2*vpa(pi)*double(freq_mat)'));
gradV_velocity= double(h./gradV_time_delay);
E_velocity= double(h./E_time_delay);

gradV_velocity_adm= double(gradV_velocity/va_num);
E_velocity_adm= double(E_velocity/va_num);

% 2D time
TwoD_time= joule_time*10.^2;
TwoD_time_norm= TwoD_time*0.5;

%% Phase
% r_ampl= zeros(nb_omega,1);
% t_delay= zeros(nb_omega,1);
% for i = 1:nb_omega
% [i_lag_bot_i, ampl_bot_i] =find(double(gradV_bot(i,:)) == max(double(gradV_bot(i,:))));
% [i_lag_top_i, ampl_top_i] =find(double(gradV_top(i,:)) == max(double(gradV_top(i,:))));
% r_ampl(i)= ampl_top_i/ampl_bot_i;
% t_delay= double(time(i,i_lag_top_i)) - double(time(i,i_lag_bot_i));
% t_delay_adm= t_delay*double(Freq_Alfven);
% i
% end
% 

%% save
selpath = uigetdir(['C:\Users\lalloz\Documents\these\onde alfven\etude theorique numerique',...
    '\1_electrode_model_results\one plate forcing case\freq_sweeping']);
save([selpath,'\donnees.mat'],'Bb_ri','gradV_env','gradV_phase','Freq_Alfven',...
    'Ha_num','R','S_num','freq_mat','freq_max','freq_min','va_num',...
'E_env','E_phase','kt_adm','kt_max','delta_kt','nb_kt','r_investigation','sigma_r')  

%% Figure bot & top gradV amplitude vs frequency

figure

subplot(2,1,1)  %yyaxis left
plot(xvect,gradV_ampl_ratio)
hold on 

xlabel(label_x)
ylabel(sprintf('$|\\nabla\\phi (z = 1)|/|\\nabla\\phi (z = 0)|$'))
ax= gca;
ax.TickLabelInterpreter="latex";
title(sprintf(['Potentiel gradient amplitude ratio in function of %s \n'...
    ' ($Ha= %4.1f \\times 10^%i, r= %4.2f$)'],label_tle,mantisse_Ha,exp_Ha,double(r_investigation))) 
txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
%t= text(700,0.4,txt);
 %t.Units= 'normalized';
grid on
 
subplot(2,1,2) 
%hold on
plot(xvect,gradV_max_bot)
hold on
plot(xvect,gradV_max_top)
ylabel(sprintf('$|\\nabla\\phi(z)|$ '))
xlabel(label_x)
ax= gca;
ax.TickLabelInterpreter="latex";
grid on

legend('$z = 0$','$z = 1$','Interpreter','Latex')
title(sprintf(['Potentiel gradient amplitude in function of %s \n'...
    ' ($Ha= %4.1f \\times 10^%i, r= %4.2f$)'],label_tle,mantisse_Ha,exp_Ha,double(r_investigation))) 
set(gca, 'YScale', 'log')

%% Figure spatial attenuation coefficient
figure('name','attenuation coef GradV versus Reta')
%yyaxis left
plot(xvect,gradV_att_coef)
hold on
xlabel(label_x)
ylabel(sprintf('$\\alpha\\left(r = %4.1f, z= 1\\right)$',r_investigation))
l_fig= legend(sprintf('$Ha= %4.2e$',Ha_num),'interpreter','latex');
title(sprintf(['Spatial attenuation coefficient for the potential gradient\n'...
' ($S= %4.1f, r= %4.2f$)'],double(S_num),double(r_investigation)))
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

%% Figure phase shift versus Freq

figure('name','phase shift versus  versus F/Falfven')
%yyaxis left
plot(xvect,gradV_dephasage)
hold on
xlabel(label_x)
ylabel(sprintf('$\\varphi_{\Delta}\\left(r = %4.1f, z= 1\\right)$',r_investigation))
l_fig= legend(sprintf('$Ha= %4.2e$',Ha_num),'interpreter','latex');
title(sprintf(['Spatial attenuation coefficient for the potential gradient\n'...
' ($Ha= %4.1f \\times 10^%i, r= %4.2f$)'],label_tle,mantisse_Ha,exp_Ha,double(r_investigation))) 
%txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
%t= text(700,0.4,txt);
%t.Units= 'normalized';
f_att_coef_ax= gca;
f_att_coef_ax.TickLabelInterpreter= 'latex';
f_att_coef_ax.FontSize= 12;
f_att_coef_ax.YLabel.FontSize = 14;
f_att_coef_ax.XLabel.FontSize = 14;
l_fig.FontSize= 12;
grid on


%%
figure('name','phase shift versus Reta')
%yyaxis left
plot(xvect,gradV_dephasage)
hold on
xlabel(label_x)
ylabel(sprintf('$\\mathcal{D}\\varphi\\left(r = %4.1f, z= 1\\right)$',r_investigation))
l_fig= legend(sprintf('$Ha= %4.2e$',Ha_num),'interpreter','latex');
title(sprintf(['Spatial attenuation coefficient versus %s\n'...
' ($Ha= %4.1f \\times 10^%i, r= %4.2f$)'],label_tle,mantisse_Ha,exp_Ha,double(r_investigation))) 
%txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
%t= text(700,0.4,txt);
%t.Units= 'normalized';
f_att_coef_ax= gca;
f_att_coef_ax.TickLabelInterpreter= 'latex';
f_att_coef_ax.FontSize= 12;
f_att_coef_ax.YLabel.FontSize = 14;
f_att_coef_ax.XLabel.FontSize = 14;
l_fig.FontSize= 12;
grid on


%% Figure transit time scaled with the alfven time coefficient
figure
%yyaxis left
plot(freq_mat*TwoD_time,gradV_time_delay*Freq_Alfven)%/Freq_Alfven
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


%% Figure bot & top gradV phases vs frequency

figure
subplot(2,1,1)  %yyaxis left
plot(freq_mat,gradV_dephasage)
xlabel('$\tilde{f}_{inj}$ (Hz)')
ylabel(sprintf('$\\Delta\\,\\varphi \\left( \\nabla\\phi \\right)\\;$ (rad)'))

title(sprintf(['Potential gradient phase shift in function of the signal frequency'...
    ' ($S= %4.1f$, $r= %4.2f$)'],double(S_num),double(r_investigation)))
% txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
% t= text(700,0.04,txt);
% t.Units= 'normalized';

subplot(2,1,2) 
plot(freq_mat,gradV_phase_bot)
hold on
plot(freq_mat,gradV_phase_top)
ylabel(sprintf('$\\varphi\\left(\\nabla\\phi (z)\\right)\\;$ (rad) '))
xlabel('$\tilde{f}_{inj}$ (Hz)')

legend('$(z=0)$','$(z=1)$','Interpreter','Latex')
title(sprintf(['Potential gradient phase in function of the signal frequency'...
    ' (S= %4.1f, $r= %4.2f$)'],double(S_num),double(r_investigation)))

%% Figure Phase difference

figure
plot(freq_mat,gradV_dephasage)
hold on 
xlabel('$\tilde{f}_{inj}$ (Hz)')
ylabel(sprintf('$\\mathcal{D}\\varphi$'))
title(sprintf(['Phase difference of the signal in function of the signal frequency'...
    '\n($S= %4.1f$, $r= %4.2f$)'],double(S_num),double(r_investigation))) 


%% Figure transit time

% physic transit time
figure
plot(freq_mat,gradV_time_delay)
hold on 
xlabel('$\tilde{f}_{inj}$ (Hz)')
ylabel(sprintf('$\\tau_{transit}$ (s)'))
title(sprintf(['transit time of the signal in function of the signal frequency'...
    '\n($S= %4.1f$, $r= %4.2f$)'],double(S_num),double(r_investigation))) 
%txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);

% scaled with 2Dtime transit time
figure
plot(freq_mat,gradV_time_delay/TwoD_time_norm)
hold on 
xlabel('$\tilde{f}_{inj}$ (Hz)')
ylabel(sprintf('$\\tau_{transit}$/$\\tau_{2D}$ '))
title(sprintf(['Scaled transit time of the signal in function of the signal frequency'...
    '\n($S= %4.1f$, $r= %4.2f$)'],double(S_num),double(r_investigation))) 
%txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
set(gca,'YScale', 'log')

% scaled with Alfven_time transit time
figure
plot(freq_mat,gradV_time_delay/alfven_time)
hold on 
xlabel('$\tilde{f}_{inj}$ (Hz)')
ylabel(sprintf('$\\tau_{transit}$/$\\tau_{Alfven}$ '))
title(sprintf(['Scaled transit time of the signal in function of the signal frequency'...
    '\n($S= %4.1f$, $r= %4.2f$)'],double(S_num),double(r_investigation))) 
%txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
set(gca,'YScale', 'log')
%% Figure gradV velocity

% physic velocity
figure
plot(freq_mat,-gradV_velocity)
hold on 
xlabel('$\tilde{f}_{inj}$ (Hz)')
ylabel(sprintf('$v_\\varphi$ (m/s)'))
title(sprintf(['Velocity of the wave in function of the signal frequency'...
    '\n($S= %4.1f$, $r= %4.2f$)'],double(S_num),double(r_investigation))) 
%txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);

% scaled velocity adm
figure
plot(freq_mat,-gradV_velocity/va_num)
hold on 
xlabel('$\tilde{f}_{inj}$ (Hz)')
ylabel(sprintf('$v_\\varphi/v_{Alfven}$'))
title(sprintf(['Velocity of the wave non scaled with the Alfven velocity'...
    '\nin function of the signal frequency'...
    '\n($S= %4.1f$, $r= %4.2f$)'],double(S_num),double(r_investigation))) 
%txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

%% Figure bot & top E amplitude vs frequency

figure
subplot(2,1,1)  %yyaxis left
plot(xvect,E_ampl_ratio)
xlabel(label_x)
ylabel(sprintf('$|E(z = 1)|/|E (z = 0)|$'))
grid on
title(sprintf(['Electric field amplitude ratio in function of the signal frequency'...
    ' ($S= %4.1f$, $r= %4.2f$)'],double(S_num),double(r_investigation))) 
% txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
% t= text(700,0.04,txt);
% t.Units= 'normalized';

subplot(2,1,2) 
plot(xvect,E_max_bot)
hold on
plot(xvect,E_max_top)
ylabel(sprintf('$|E(z)|$'))
xlabel(label_x)
grid on
legend('$z=0$','$z=1$','Interpreter','Latex')
title(sprintf(['Electric field amplitude in function of the signal frequency'...
    ' (S= %4.1f, $r= %4.2f$)'],double(S_num),double(r_investigation))) 
set(gca, 'YScale', 'log')

%% Figure bot & top E phases vs frequency

figure
subplot(2,1,1)  %yyaxis left
plot(xvect,E_dephasage)
xlabel('$\tilde{f}_{inj}$ (Hz)')
ylabel(sprintf('$\\Delta\\,\\varphi \\left( E \\right)\\;$ (rad)'))


title(sprintf(['Electric field phase shift in function of the signal frequency'...
    ' ($S= %4.1f$, $r= %4.2f$)'],double(S_num),double(r_investigation))) 
% txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
% t= text(700,0.04,txt);
% t.Units= 'normalized';

subplot(2,1,2) 
plot(xvect,E_phase_bot)
hold on
plot(xvect,E_phase_top)
ylabel(sprintf('$\\varphi\\left(E(z)\\right)\\;$ (rad)'))
xlabel('$\tilde{f}_{inj}$ (Hz)')

legend('$z=0$','$z=1$','Interpreter','Latex')
title(sprintf(['Electric field phase in function of the signal frequency'...
    ' (S= %4.1f, $r= %4.2f$)'],double(S_num),double(r_investigation))) 

%% Figure transit time E

% physic transit time
figure
plot(freq_mat,E_time_delay)
hold on 
xlabel('$\tilde{f}_{inj}$ (Hz)')
ylabel(sprintf('$\\tau_{transit}$ (s)'))
title(sprintf(['transit time of E in function of the signal frequency'...
    '\n($S= %4.1f$, $r= %4.2f$)'],double(S_num),double(r_investigation))) 
%txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);

% scaled with 2Dtime transit time
figure
plot(freq_mat,E_time_delay/TwoD_time_norm)
hold on 
xlabel('$\tilde{f}_{inj}$ (Hz)')
ylabel(sprintf('$\\tau_{transit}$/$\\tau_{2D}$ '))
title(sprintf(['Scaled transit time of E in function of the signal frequency'...
    '\n($S= %4.1f$, $r= %4.2f$)'],double(S_num),double(r_investigation))) 
%txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
set(gca,'YScale', 'log')

% scaled with Alfven_time transit time
figure
plot(freq_mat,E_time_delay/alfven_time)
hold on 
xlabel('$\tilde{f}_{inj}$ (Hz)')
ylabel(sprintf('$\\tau_{transit}$/$\\tau_{Alfven}$ '))
title(sprintf(['Scaled transit time of the signal in function of the signal frequency'...
    '\n($S= %4.1f$, $r= %4.2f$)'],double(S_num),double(r_investigation))) 
%txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
set(gca,'YScale', 'log')

%% Figure E velocity

figure
plot(freq_mat,E_velocity)
hold on 
xlabel('$\tilde{f}_{inj}$ (Hz)')
ylabel(sprintf('$v_\\varphi$ (m/s)'))
title(sprintf(['Velocity of the wave in function of the signal frequency'...
    '\n($S= %4.1f$, $r= %4.2f$)'],double(S_num),double(r_investigation))) 
%txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);

%% Figure E velocity adm

figure
plot(freq_mat,E_velocity/va_num)
hold on 
xlabel('$\tilde{f}_{inj}$ (Hz)')
ylabel(sprintf('$v_\\varphi / v_{Alfven}$'))
title(sprintf(['Phase velocity scaled with the Alfven velocity'...
    '\nin function of the signal frequency'...
    '\n($S= %4.1f$, $r= %4.2f$)'],double(S_num),double(r_investigation))) 
%txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);

%% Figure amplitude ratio vs frequency
figure

plot(freq_mat/Freq_Alfven,E_ampl_ratio,'-*')
ylabel('amplitude ratio')
xlabel('$f_{inj}$')

legend('electric field ratio')
title(sprintf(['Amplitude ratio between the upper and the bottom signal \n',... at %4.2f cm \n from the centre of the electrode'...
    ' (S= %4.1f , Pm= %4.2e, Ha= %4.2e)'],double(S_num),double(Pm_num),double(Ha_num)));%double(r_investigation)*10,double(S_num),double(Pm_num),double(Ha_num)))
% txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
% texte= text(1000,0.78,txt);
% texte.Units= 'normalized';

%% amplitude ratio vs frequency, sump up of different S number
Freq_Alfven_S_14p5= 3/sqrt(rho_gal*mu_0)/L_carac;
%gradV_ampl_ratio_S_53= gradV_ampl_ratio; % à faire pour chaque S
%Freq_S_53= freq_mat;
%% 
%adm freq sweeping figure

figure('Name','Fig_gradV_ratio_adm_freq_sweeping_diff_S')
plot(Freq_S_14p5/Freq_Alfven_S_14p5,gradV_ampl_ratio_S_14p5,'-*')
ylabel('amplitude ratio')
xlabel('signal frequency (Hz)')
title(sprintf(['Amplitude ratio of the potential gradient between the upper and the bottom signal at %4.2f cm \n from the centre of the electrode'...
' ($Pm= %4.2e, Ha= %4.2e$)'],double(r_investigation)*10,double(Pm_num),double(Ha_num)))
txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
texte= text(1000,0.78,txt);
hold on
plot(Freq_S_24/Freq_Alfven_S_24,gradV_ampl_ratio_S_24,'-*')
hold on
plot(Freq_S_33/Freq_Alfven_S_33,gradV_ampl_ratio_S_33,'-*')
hold on
plot(Freq_S_53/Freq_Alfven_S_53,gradV_ampl_ratio_S_53,'-*')
legend('S = 14,5','S = 24','S = 33','S = 53')

%% freq sweeping figure

figure('Name','Fig_gradV_ratio_freq_sweeping_diff_S')
plot(Freq_S_14p5,gradV_ampl_ratio_S_14p5,'-*')
ylabel('amplitude ratio')
xlabel('signal frequency (Hz)')
title(sprintf(['Amplitude ratio of the potential gradient between the upper and the bottom signal at %4.2f cm \n from the centre of the electrode'...
' ($Pm= %4.2e, Ha= %4.2e$)'],double(r_investigation)*10,double(Pm_num),double(Ha_num)))
txt= sprintf('$ \\nabla\\phi_0 = %4.2e \\,V/m $',E0);
texte= text(1000,0.78,txt);
hold on
plot(Freq_S_24,gradV_ampl_ratio_S_24,'-*')
hold on
plot(Freq_S_33,gradV_ampl_ratio_S_33,'-*')
hold on
plot(Freq_S_53,gradV_ampl_ratio_S_53,'-*')
legend('S = 14,5','S = 24','S = 33','S = 53')

%%
figure
plot(freq_mat,gradV_max_bot)
ylabel('amplitude of the potential gradient')
xlabel('signal frequency (Hz)')

%% figure enveloppe et phase

    if show_fig == 1
    Fig_enveloppe_phase_onde= figure('name','fig_enveloppe_onde');
    subplot(1,2,1)
    %axis([-Bmax Bmax 0 hadm])
    title('Electric potential')
    xlabel('$V/ \left(b_0\,L_{carac}\,\omega\right)$')
    ylabel('$z/L_{carac}$')
    hold on
    plot(benv,z,'r',-benv,z,'r')

    subplot(1,2,2)
    %axis([-Vmax Vmax 0 hadm])
    plot(bphase,z,'k')
    ylabel('$z/L_{carac}$')
    xlabel('Phase to forcing (rad)')
    legend('Electric potential')
    xticks(-pi/2:pi/4:pi/2) 
    xticklabels({'-\pi/2','-\pi/4','0','\pi/4','\pi/2'})
    end 

    if show_fig == 1
    Fig_enveloppe_phase_onde= figure('name','fig_enveloppe_onde');
    subplot(1,2,1)
    %axis([-Bmax Bmax 0 hadm])
    title('Electric potential')
    xlabel('$V/ \left(b_0\,L_{carac}\,\omega\right)$')
    ylabel('$z/L_{carac}$')
    hold on
    plot(venv,z,'r',-venv,z,'r')

    subplot(1,2,2)
    %axis([-Vmax Vmax 0 hadm])
    plot(vphase,z,'k')
    ylabel('$z/L_{carac}$')
    xlabel('Phase to forcing (rad)')
    legend('Electric potential')
    xticks(-pi/2:pi/4:pi/2) 
    xticklabels({'-\pi/2','-\pi/4','0','\pi/4','\pi/2'})
    end 
    
    
% extraction ratio ampl / phase
    disp('debut calcul ration ampl et dephasage')
    ratio_ampl_mat(j)= venv(end)/venv(1);
    phase_bot_mat(j)= vphase(1);
    phase_up_mat(j)= vphase(end);
    dephasage_mat(j)= vphase(end)-vphase(1);
    disp('fin calcul ration ampl et dephasage')


%%
dephasage_bis_mat= dephasage_mat;
for i= 1:19
    if phase_up_mat(i+1)*phase_up_mat(i)<0 && phase_up_mat(i+1)>0 && phase_bot_mat(i+1)<0
        dephasage_bis_mat(i+1)= dephasage_mat(i+1)- sym(pi);
    end
end
    
%%
set(0,'defaultTextInterpreter','latex')
opengl software
Fig_ampl_ratio_phase_sft= figure('name','fig_ratio_phase_shift');
subplot(2,1,1)
%axis([-Bmax Bmax 0 hadm])
title('Amplitude ratio between the upper and the bottom signal')
xlabel('signal frequency (Hz)')
ylabel('amplitude ratio')
hold on
plot((1/(2*pi*tau))*epsilon_mat,ratio_ampl_mat)

subplot(2,1,2)
%axis([-Vmax Vmax 0 hadm])
plot((1/(2*pi*tau))*epsilon_mat,dephasage_bis_mat)
title('Phase shift between the upper and the bottom signal')
ylabel('phase shift')
xlabel('signal frequency (Hz)')
yticks(-pi/2:pi/4:pi/2) 
yticklabels({'-\pi/2','-\pi/4','0','\pi/4','\pi/2'})

