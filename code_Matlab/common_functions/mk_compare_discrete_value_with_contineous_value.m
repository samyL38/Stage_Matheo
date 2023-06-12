%%%
% Comparaison des mesures attendues des gradients de potentiel en fonction 
% de points de mesures dicrets en fréquence avec l'évolution continue en 
% fréquences de ces gradients de potentiel, obtenue a partir du modèle
%%%%



Cfig= gcf  % doit être une figure qui montre le coefficient d'att ainsi que
            % les gradients de potentiel en z=0 et z=1;

% frequnces discrete
Reta_exp_global_mat= [57:40:1.2*double(Freq_Alfven)]*double(tau)*2*pi;
Reta_exp_zoom_mat= [double(Freq_Alfven)/2 - 50:10:...
    double(Freq_Alfven)/2 + 50]*double(tau)*2*pi;
Reta_exp_low_mat= [3:5:63]*double(tau)*2*pi;

nb_global_pt= 40;
nb_zoom_pt= 15;
nb_low_pt= 20;

Reta_exp_global_mat= linspace(0.08*Freq_Alfven,1.2*Freq_Alfven,nb_global_pt)*double(tau)*2*pi;

Reta_exp_zoom_mat= [linspace(Freq_Alfven*(1/2 -0.07),Freq_Alfven*(1/2+0.07),nb_zoom_pt),...
    linspace(Freq_Alfven*(1 -0.06),Freq_Alfven*(1+0.06),nb_zoom_pt)] *double(tau)*2*pi;

Reta_exp_low_mat= linspace(2,0.08*Freq_Alfven,nb_low_pt)*double(tau)*2*pi;


% attenuation coefficient
att_exp_global_mat= interp1(double(Cfig.Children(3).Children.XData),...
    Cfig.Children(3).Children.YData,[Reta_exp_low_mat,Reta_exp_global_mat],"cubic");
att_exp_zoom_mat= interp1(double(Cfig.Children(3).Children.XData),...
    Cfig.Children(3).Children.YData,Reta_exp_zoom_mat,"cubic");
att_exp_low_mat= interp1(double(Cfig.Children(3).Children.XData),...
    Cfig.Children(3).Children.YData,[Reta_exp_low_mat,Reta_exp_global_mat],"cubic");
%%
% gradient de potentiel
phi_z_1_exp_global_mat= interp1(double(Cfig.Children(2).Children(1).XData),...
    Cfig.Children(2).Children(1).YData,[Reta_exp_low_mat,Reta_exp_global_mat],"cubic");
phi_z_1_exp_zoom_mat= interp1(double(Cfig.Children(2).Children(1).XData),...
    Cfig.Children(2).Children(1).YData,Reta_exp_zoom_mat,"cubic");
phi_z_1_exp_low_mat= interp1(double(Cfig.Children(2).Children(1).XData),...
    Cfig.Children(2).Children(1).YData,[Reta_exp_low_mat,Reta_exp_global_mat],"cubic");

phi_z_0_exp_global_mat= interp1(double(Cfig.Children(2).Children(2).XData),...
    Cfig.Children(2).Children(2).YData,[Reta_exp_low_mat,Reta_exp_global_mat],"cubic");
phi_z_0_exp_zoom_mat= interp1(double(Cfig.Children(2).Children(2).XData),...
    Cfig.Children(2).Children(2).YData,Reta_exp_zoom_mat,"cubic");
phi_z_0_exp_low_mat= interp1(double(Cfig.Children(2).Children(2).XData),...
    Cfig.Children(2).Children(2).YData,[Reta_exp_low_mat,Reta_exp_global_mat],"cubic");

%% figure medium freq
figure
plot([Reta_exp_low_mat,Reta_exp_global_mat,Reta_exp_zoom_mat],...
    [phi_z_0_exp_global_mat, phi_z_0_exp_zoom_mat],'+',"MarkerFaceColor" ...
 ,"#0072BD",'LineWidth',1)
hold on
plot([Reta_exp_low_mat,Reta_exp_global_mat,Reta_exp_zoom_mat],...
    [phi_z_1_exp_global_mat, phi_z_1_exp_zoom_mat],'+',"MarkerFaceColor" ...
 ,"#D95319",'LineWidth',1)

figure
plot([Reta_exp_low_mat,Reta_exp_global_mat,Reta_exp_zoom_mat],...
    [att_exp_global_mat, att_exp_zoom_mat],'k+','LineWidth',1)

%% low frequency
figure
plot([Reta_exp_low_mat,Reta_exp_global_mat],phi_z_0_exp_low_mat,'o',"MarkerEdgeColor" ...
 ,"#0072BD",'LineWidth',1)
hold on
plot([Reta_exp_low_mat,Reta_exp_global_mat],phi_z_1_exp_low_mat,'o',"MarkerEdgeColor" ...
 ,"#D95319",'LineWidth',1)

figure
plot([Reta_exp_low_mat,Reta_exp_global_mat],att_exp_low_mat,'ko','LineWidth',1)

