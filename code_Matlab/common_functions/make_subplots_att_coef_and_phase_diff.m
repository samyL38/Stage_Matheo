%%%
% Make plots attenuation coef and phase diff
%%%
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

%% cumulative parameter
Ha_vect=[];
Reta_vect= [];


%% Import data

load("donnees.mat",'B_env', 'B_phase', 'Ha_num', 'Reta', 'Rnu', 'epsilon',...
    'frequence_forcage', 'gradV_env', 'gradV_phase','hadm','nb_point_r',...
    'nb_point_z', 'r_investigation', 'r_max', 'r_min')


%% Parameters

%Paramètres physiques
mu_0= vpa(1.25*sym(1e-6));
sigma_gal= vpa(3.46*sym(1e6));
visco_gal= vpa(0.0024); %Viscosité dynamique du Galstn en Pa.s
eta= 1/(sigma_gal * mu_0); %diffusivité mhd
rho_gal= vpa(6440);% masse volumique Galinstan kg/m3
nu_gal= visco_gal/rho_gal;

% geometric parameter
h= vpa(0.1); % distance entre la plaque inf et sup en mètre
hadm= h/h;

z_vect= double(0:hadm/(nb_point_z-1):hadm);
[R_mat , Z_mat]= meshgrid(double(r_investigation),z_vect);
nb_line_env= 40;
nb_line_phase= 20;
r_inv_lin= linspace(r_investigation(1),r_investigation(end),20);

attenuation_coef= log(transpose(gradV_env./gradV_env(1:end,1)));
phase_diff= (gradV_phase-gradV_phase(:,1))';

linecontour_gradV_env= logspace(log10(1e-2),...
    log10(5e5),nb_line_env);
linecontour_gradV_env_norm = logspace(log10(1e-8),...
    log10(1),nb_line_env);
linecontour_gradV_phase= linspace(min(double(gradV_phase),[],'all',"omitnan"),...
    max(double(gradV_phase),[],'all'),nb_line_phase);
linecontour_attenuation_coef_env= 20;%linspace(-9,1,20);%real([linspace(-12,...
%-0.1,15), linspace(0.1,1,10)]);%4
linecontour_phase_diff= 20;
%linecontour_A_env= logspace(-2,...
%    2,nb_line_env);

% limits for the colorbar
gradV_ampl_inf= 1;%min(abs(double(gradV_env(:,1:end-1))),[],'all');
gradV_ampl_sup= 10^4;%max(abs(double(gradV_env(:,1:end-1))),[],'all');

depth_penetration_eta= 1/sqrt(vpa(pi)*frequence_forcage/eta);
depth_penetration_eta_adm= depth_penetration_eta/h;

%% add cumulative param
Ha_vect= [Ha_vect, Ha_num];
Reta_vect= [Reta_vect, Reta];

%% contours attenaution coef alone single Ha, Reta

figure('Name','attenuation coef vs r and z');
s=contourf(R_mat(:,1:end),Z_mat(:,1:end),attenuation_coef,linecontour_attenuation_coef_env);%linecontour_attenuation_coef_env);
hold on
plot(1/sqrt(Ha_num)*ones(nb_point_z,1),z_vect,'--w','linewidth',1.5)
hold on
plot(double(r_inv_lin),double(depth_penetration_eta_adm*ones(length(r_inv_lin),1)),'.k','linewidth',1.5)
hold on
contour(R_mat,Z_mat,attenuation_coef,[0 0],'k','linewidth',2);
colormap(BWR2(100))
c= colorbar;
caxis([-14,14]);
c.Limits=([min(double(attenuation_coef),[],'all'),...
    max(double(attenuation_coef),[],'all')]); %[linecontour_attenuation_coef_env(1),linecontour_attenuation_coef_env(end)] ;%
c.Label.String = '$\alpha \left( r, z\right)$';
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
title(sprintf(['Mapping of the radial potential gradient amplitude for \n $Ha= %4.2e,\\, R_{\\eta}/2\\pi = %4.2f,',...
           ' \\,R_\\nu/2\\pi = %4.2e$'],Ha_num,epsilon/(2*vpa(pi)),Rnu/(2*vpa(pi))))


%% figure for attenuation coefficient
figure
tiledlayout(3,4)
Cfig= gcf;
%%
nexttile()
contourf(R_mat(:,1:end),Z_mat(:,1:end),attenuation_coef,linecontour_attenuation_coef_env);%linecontour_attenuation_coef_env);
hold on
plot(1/sqrt(Ha_num)*ones(nb_point_z,1),z_vect,'--w','linewidth',1.5)
hold on
plot(double(r_inv_lin),double(depth_penetration_eta_adm*ones(length(r_inv_lin),1)),'.k','linewidth',1.5)
hold on
contour(R_mat,Z_mat,attenuation_coef,[0 0],'k','linewidth',2);
colormap(BWR2(100))
c= colorbar;
caxis([-14,14]);
c.Label.Interpreter= 'latex';
c.Label.FontSize= 14;
c.TickLabelInterpreter='latex';
f= gca;
f.ColorScale= 'linear';
f.TickLabelInterpreter= 'latex';
%title(sprintf(['Mapping of the radial potential gradient amplitude for \n $Ha= %4.2e,\\, R_{\\eta}/2\\pi = %4.2f,',...
%           ' \\,R_\\nu/2\\pi = %4.2e$'],Ha_num,epsilon/(2*vpa(pi)),Rnu/(2*vpa(pi))))

colorbar('off')
%%
cb= colorbar;
cb.Layout.Tile = 'north';
Cfig.Children.TileSpacing= "tight";
Cfig.Children.Padding= "loose";

%% Ha and Reta labels

% addings of label and content characteristics for attenuation coef
for i= 2:13
    Cfig.Children.Children(i).FontSize= 10;
    Cfig.Children.Children(i).TickLabelInterpreter="latex";
end

Cfig.Children.Children(1).TickLabelInterpreter="latex";
Cfig.Children.Children(1).FontSize=10;
Cfig.Children.Children(1).Label.String= '$\alpha \left( r,z\right)$';
Cfig.Children.Children(1).Label.Interpreter="latex";
Cfig.Children.Children(1).Label.FontSize= 14;

for i= 2:5
    Cfig.Children.Children(i).XLabel.String= '$r$';
    Cfig.Children.Children(i).XLabel.Interpreter="latex";
    Cfig.Children.Children(i).XLabel.FontSize= 14;
end

for i= [5,9,13]
    Cfig.Children.Children(i).YLabel.String= '$z$';
    Cfig.Children.Children(i).YLabel.Interpreter="latex";
    Cfig.Children.Children(i).YLabel.FontSize= 14;
    Cfig.Children.Children(i).YLabel.Rotation= 0;
    Cfig.Children.Children(i).YLabel.HorizontalAlignment="right";
end

for i=[2,6,10]
    Cfig.Children.Children(i).Subtitle.String=sprintf('$Ha = %4.0f$',Ha_vect(end-i+2)) ;
    Cfig.Children.Children(i).Subtitle.Position=[1.2,0.5,0];
    Cfig.Children.Children(i).Subtitle.Units="normalized";
    Cfig.Children.Children(i).Subtitle.Interpreter="latex";
    Cfig.Children.Children(i).Subtitle.FontSize= 14;
    Cfig.Children.Children(i).Subtitle.Rotation = 90;
end

for i= 2:5
    Cfig.Children.Children(i).Title.String=sprintf('$R_\\eta = %4.0f$',Reta_vect(end-i+2)) ;
    Cfig.Children.Children(i).Title.Units="normalized";
    Cfig.Children.Children(i).Title.Position=[0.5,-.3,0];    
    Cfig.Children.Children(i).Title.Interpreter="latex";
    Cfig.Children.Children(i).Title.FontSize= 14;
    Cfig.Children.Children(i).Title.Rotation = 0;
end

%% export file
Cfig= gcf;
Cfig.Units="centimeters";
Cfig.Position=[4,-11,21,29];
%%
exportgraphics(Cfig,'Att_coef_alpha_summary.eps','ContentType','vector')

%%
%%%%%%%%
% phase difference
%%%%%%%%

% contours phase diff alone single Ha, Reta

figure('Name','gradV phase diff vs r and z')
contourf(R_mat,Z_mat,phase_diff,linecontour_phase_diff);%linecontour_gradV_phase);
hold on
contour(R_mat,Z_mat,phase_diff,[0 0],'k','linewidth',2);
colormap(BWR2(100))
c= colorbar;

caxis([-max(abs(phase_diff),[],'all'), max(abs(phase_diff),[],'all')]);
c.Limits= ([min(double(phase_diff),[],'all'),max(double(phase_diff),[],'all')]);
%c.Limits= ([min(double(transit_time),[],'all'),0.01])
c.Label.String= ['$\Delta\varphi\left[\partial_r\phi\left(r,z\right)\right]$'];
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


%% figure for phase diff
figure
tiledlayout(3,4)
Cfig= gcf;
%%
nexttile()
contourf(R_mat,Z_mat,phase_diff,linecontour_phase_diff);%linecontour_gradV_phase);
hold on
contour(R_mat,Z_mat,phase_diff,[0 0],'k','linewidth',2);
colormap(BWR2(100))
c= colorbar;
caxis([-max(abs(phase_diff),[],'all'), max(abs(phase_diff),[],'all')]);
c.Limits= ([min(double(phase_diff),[],'all'),max(double(phase_diff),[],'all')]);
c.Label.Interpreter= 'latex';
c.Label.FontSize= 14;
c.TickLabelInterpreter='latex';
f= gca;
f.ColorScale= 'linear';
f.TickLabelInterpreter= 'latex';
%title(sprintf(['Mapping of the radial potential gradient amplitude for \n $Ha= %4.2e,\\, R_{\\eta}/2\\pi = %4.2f,',...
%           ' \\,R_\\nu/2\\pi = %4.2e$'],Ha_num,epsilon/(2*vpa(pi)),Rnu/(2*vpa(pi))))

colorbar('off')
%%
cb= colorbar;
cb.Layout.Tile = 'north';
Cfig.Children.TileSpacing= "tight";
Cfig.Children.Padding= "loose";

%% Ha and Reta labels

% addings of label and content characteristics for attenuation coef
for i= 2:13
    Cfig.Children.Children(i).FontSize= 10;
    Cfig.Children.Children(i).TickLabelInterpreter="latex";
end

Cfig.Children.Children(1).TickLabelInterpreter="latex";
Cfig.Children.Children(1).Ticks= [-6*pi  -5*pi -4*pi -3*pi -5*pi/2 -2*pi ...
    -7*pi/4 -3*pi/2 -5*pi/4 -pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi];
Cfig.Children.Children(1).TickLabels= {'$-6\pi$', '$-5\pi$','$-4\pi$',...
    '$-3\pi$', '$-5\pi/2$','$-2\pi$', '$-7\pi/4$', '$-3\pi/2$','$-5\pi/4$',...
    '$-\pi$', '$-3\pi/4$','$-\pi/2$', '$-\pi/4$', '$0$', '$\pi/4$',...
    '$\pi/2$', '$3\pi/4$', '$\pi$'};
Cfig.Children.Children(1).FontSize=10;
Cfig.Children.Children(1).Label.String= '$\Delta\varphi\left[\partial_r\phi\left(r,z\right)\right]$';
Cfig.Children.Children(1).Label.Interpreter="latex";
Cfig.Children.Children(1).Label.FontSize= 14;

for i= 2:5
    Cfig.Children.Children(i).XLabel.String= '$r$';
    Cfig.Children.Children(i).XLabel.Interpreter="latex";
    Cfig.Children.Children(i).XLabel.FontSize= 14;
end

for i= [5,9,13]
    Cfig.Children.Children(i).YLabel.String= '$z$';
    Cfig.Children.Children(i).YLabel.Interpreter="latex";
    Cfig.Children.Children(i).YLabel.FontSize= 14;
    Cfig.Children.Children(i).YLabel.Rotation= 0;
    Cfig.Children.Children(i).YLabel.HorizontalAlignment="right";
end

for i=[2,6,10]
    Cfig.Children.Children(i).Subtitle.String=sprintf('$Ha = %4.0f$',Ha_vect(end-i+2)) ;
    Cfig.Children.Children(i).Subtitle.Position=[1.2,0.5,0];
    Cfig.Children.Children(i).Subtitle.Units="normalized";
    Cfig.Children.Children(i).Subtitle.Interpreter="latex";
    Cfig.Children.Children(i).Subtitle.FontSize= 14;
    Cfig.Children.Children(i).Subtitle.Rotation = 90;
end

for i= 2:5
    Cfig.Children.Children(i).Title.String=sprintf('$R_\\eta = %4.0f$',Reta_vect(end-i+2)) ;
    Cfig.Children.Children(i).Title.Units="normalized";
    Cfig.Children.Children(i).Title.Position=[0.5,-.3,0];    
    Cfig.Children.Children(i).Title.Interpreter="latex";
    Cfig.Children.Children(i).Title.FontSize= 14;
    Cfig.Children.Children(i).Title.Rotation = 0;
end

%% export file
Cfig= gcf;
Cfig.Units="centimeters";
Cfig.Position=[4,-11,21,29];
%%
exportgraphics(Cfig,'phase_diff_summary.eps','ContentType','vector')

%% for given Ha
Cfig= gcf;
Cfig.Units="centimeters";

%%
Cfig_phase;
exportgraphics(Cfig_phase,'phase_diff_Ha_38e4.pdf','ContentType','vector')

%%
Cfig_alpha;
exportgraphics(Cfig_alpha,'att_coef_Ha_38e4.pdf','ContentType','vector')
