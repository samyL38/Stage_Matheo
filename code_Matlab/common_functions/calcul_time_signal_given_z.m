function [signal]=calcul_time_signal_given_z(A_num,s1_num,...
    k1_num,s2_num,k2_num,f,z,t,precision,phase_num)

% calcul_signal_boundaries calcul l'évolution temporelle du signal A au
% niveau des plaques inférieure et supérieure. A est une quantité décrivant des ondes
% d'Alfven de la forme :
%A(t,z) = exp(s1*z)*(A(1)*cos(t+k1*z)-A(2)*sin(t+k1*z))+
%	+exp(-s1*z)*(A(3)*cos(t-k1*z)-A(4)*sin(t-k1*z))+
%	+exp(s2*z)*(A(5)*cos(t+k2*z)-A(6)*sin(t+k2*z))+
%	+exp(-s2*z)*(A(7)*cos(t-k2*z)-A(8)*sin(t-k2*z))+
% A est notamment une vitesse, un champ magnetique induit ou un potentiel
% électrique
% Samy Lalloz 21/01/2021 


% Parameters
%npts=100;
% Definition des variables symboliques
syms A s1 k1 s2 k2 dephasage

if ~exist('precision','var') || isempty(precision)
    precision= 16;
end
digits(precision)

A = vpa(A_num);
s1 = vpa(s1_num);
k1 = vpa(k1_num);
s2 = vpa(s2_num);
k2 = vpa(k2_num);


% envelope and phase
if ~exist('phase_num','var')
    phase= sym(0);
elseif exist('phase_num','var')
    phase = vpa(phase_num);
end

signal_i=  exp(z*s1).*(A(:,1).*cos(f*t+phase+z*k1)-A(:,2).*sin(f*t+phase+z*k1))+...
        exp(-z*s1).*(A(:,3).*cos(f*t+phase-z*k1)-A(:,4).*sin(f*t+phase-z*k1))+...
        exp(z*s2).*(A(:,5).*cos(f*t+phase+z*k2)-A(:,6).*sin(f*t+phase+z*k2))+...
        exp(-z*s2).*(A(:,7).*cos(f*t+phase-z*k2)-A(:,8).*sin(f*t+phase-z*k2));
if length(k1) > 1
    signal=sum(signal_i);
else
    signal= signal_i;
end
end



