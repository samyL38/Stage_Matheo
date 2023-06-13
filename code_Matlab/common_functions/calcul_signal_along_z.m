function [signal]=calcul_signal_along_z(A_num,s1_num,k1_num,s2_num,k2_num,f,h,t,precision,phase_num)

%calcul_signal_along_z calcule l'amplitude de l'oscillation d'une quantité
% A le long de h à n instant t. A décrit des ondes Alfven de la forme
%A(t,z) = exp(s1*z)*(A(1)*cos(t+k1*z)-A(2)*sin(t+k1*z))+
%	+exp(-s1*z)*(A(3)*cos(t-k1*z)-A(4)*sin(t-k1*z))+
%	+exp(s2*z)*(A(5)*cos(t+k2*z)-A(6)*sin(t+k2*z))+
%	+exp(-s2*z)*(A(7)*cos(t-k2*z)-A(8)*sin(t-k2*z))+


% Typcally A est une vitesse, un champ magnetique induit ou un potentiel 
% électrique 
% Samy Lalloz 21/01/21

% Parameters
%npts=100;
% Definition des variables symboliques

    

syms A s1 k1 s2 k2 dephasage
if ~exist('precision','var')
    precision= 16;
end
digits(precision)

A = vpa(A_num);
s1 = vpa(s1_num);
k1 = vpa(k1_num);
s2 = vpa(s2_num);
k2 = vpa(k2_num);

if ~exist('dephasage_num','var')
    phase= sym(0);
else
    phase= phase_num;
end

s=exp(h.*s1).*(A(:,1).*cos(f*t+phase+h.*k1)-A(:,2).*sin(f*t+phase+h.*k1))+...
    exp(-h.*s1).*(A(:,3).*cos(f*t+phase-h.*k1)-A(:,4).*sin(f*t+phase-h.*k1))+...
    exp(h.*s2).*(A(:,5).*cos(f*t+phase+h.*k2)-A(:,6).*sin(f*t+phase+h.*k2))+...
    exp(-h.*s2).*(A(:,7).*cos(f*t+phase-h.*k2)-A(:,8).*sin(f*t+phase-h.*k2));
if length(k1) > 1
    signal=sum(s);
else
    signal= s;
end
