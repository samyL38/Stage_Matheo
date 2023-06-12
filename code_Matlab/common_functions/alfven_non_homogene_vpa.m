function [s1, k1, s2, k2, U, B, Vr, Jr, Er,Ar, dtAr ]=alfven_non_homogene_vpa(Pm_num,S_num,h_num,kti_num,...
    Bbi,epsilon_num,precision,plate)

% Calcule la vitesse, le champ magnetique et le potentiel electrique pour
% les ondes d'Alfven, ces dernières étant existé par un forçage magnétique
% inhomogène au niveau des plaques supérieure ou inférieure.
% Si aucune plaque n'est spécifiée, la forçage est automatiquement
% localisé au niveau de la plaque inférieure
% Samy Lalloz 21/01/2021


%SYNPOSIS
%	function [U, B, V]=alfven_2w(Pm,S,h)
%
% INPUT
% Pm 	= nu/eta, Magnetic Prandtl number
% S	= V_A/(eta*omega)^(1/2) Lunquist number  
% h	=h_dimensional/L_carac, Distance between top and bottom 
%	wall normalised the height of the vessel
% epsilon = omega*tau, angular frequency normalized by the caracteristic
% magnetic diffusion time along the height of the vessel

1j; % nombre imaginaire

% Definition des variables symboliques
syms Ub Ut Bb Bt
syms Pm S kti epsilon h

% precision
if ~exist('precision','var')
    precision= 64;
end
digits(precision)

% Application valeur numérique aux nombres adm
Pm= vpa(Pm_num);
S= vpa(S_num);
kti= vpa(kti_num);
epsilon= vpa(epsilon_num);
h= vpa(h_num);

% DEFINITION OF  BOUNDARY CONDITIONS

% Amplitude of velocity oscillation at bottom and top walls

Ub=vpa(0);
Ut=vpa(0);
 
% Amplitude of magnetic field oscillations at bottom and top walls

if exist('plate','var')
    switch plate
        case {'upper','top','superior','0'}
            Bt=vpa(Bbi);%vpa(1/(2*sym(pi)) * TF_b_theta_kti.*besselj(0,kti*ri)*delta_kti);
            Bb=vpa(0);
        case {'bottom','lower','inferior','1'}
            Bb=vpa(Bbi);%vpa(1/(2*sym(pi)) * TF_b_theta_kti.*besselj(0,kti*ri)*delta_kti);
            Bt=vpa(0);
        otherwise
            error("case doesn't exist")
    end
else
    Bb=vpa(Bbi);%vpa(1/(2*sym(pi)) * TF_b_theta_kti.*besselj(0,kti*ri)*delta_kti);
    Bt=vpa(0);
end
%Bt= vpa(1/20.5);
%Diss=[Pm 0 (1+Pm)*1j-S.^2 0 -1];
CoefRD1= [Pm 0 S^2+2*Pm*kti^2-1j*epsilon*(1+Pm) 0 Pm*kti^4-1j*epsilon*(1+Pm)*kti^2-epsilon^2];
CoefRD2= [Pm 0 S^2+2*Pm*kti^2+1j*epsilon*(1+Pm) 0 Pm*kti^4+1j*epsilon*(1+Pm)*kti^2-epsilon^2];

ro=roots(CoefRD1);

[srro,isort]=sort(real(ro));
sro=ro(isort);
posro=find(real(sro)<=10^-15);
s1=abs(imag(sro(posro(1))));
s2=abs(imag(sro(posro(2))));
k1=abs(real(sro(posro(1))));
k2=abs(real(sro(posro(2))));

% BOUNDARY CONDITIONS MATRICES
% No-slip at top and bottom wall

sk1=vpa(sin(k1*h));
ck1=vpa(cos(k1*h));
ep1=vpa(exp(s1*h));
em1=vpa(exp(-s1*h));


N1= [ 0 1 0 1;...
    1 0 1 0;...
    -sk1*ep1 -ck1*ep1 sk1*em1 -ck1*em1;...
    ck1*ep1 -sk1*ep1 ck1*em1 sk1*em1];

sk2=vpa(sin(k2*h));
ck2=vpa(cos(k2*h));
ep2=vpa(exp(s2*h));
em2=vpa(exp(-s2*h));

N2= [ 0 1 0 1;...
    1 0 1 0;...
    -sk2*ep2 -ck2*ep2 sk2*em2 -ck2*em2;...
    ck2*ep2 -sk2*ep2 ck2*em2 sk2*em2];



Uw=[0  Ub 0 Ut].';
% Matrice identité alternant 1 et -1
Id4= [[1 0 0 0]; [0 -1 0 0]; [0 0 1 0]; [0 0 0 -1]];
% ==> Remarque : cette matrice est utilisée (et non la matrice identité)
% car les coffecients devant les termes sin(w*t_k*z) sont définis avec un 
% signe moins

% Matrice dérivée temporelle
Dt= [[0 -epsilon 0 0]; [epsilon 0 0 0]; [0 0 0 -epsilon]; [0 0 epsilon 0]];

% Dirichlet boundary conditions for the magnetic field at top and bottom walls
Z00=zeros(2);
Dz10=[[s1, -k1];[-k1 -s1]];
Dz1=[[Dz10 Z00];[Z00 -Dz10]];
M10=[[Pm*(k1^2-s1^2) 2*k1*s1*Pm-epsilon]; [2*k1*s1*Pm-epsilon -Pm*(k1^2-s1^2)]];
M1=[[M10 Z00];[Z00 M10]];
BU1=S^(-2)*(Dz1\(M1+Pm*kti^2*Id4));
Q1=N1*BU1;

Dz20=[[s2, -k2];[-k2 -s2]];
Dz2=[[Dz20, Z00];[Z00 -Dz20]];
M20=[[Pm*(k2^2-s2^2) 2*k2*s2*Pm-epsilon]; [2*k2*s2*Pm-epsilon -Pm*(k2^2-s2^2)]];
M2=[[M20 Z00];[Z00 M20]];
BU2=S^(-2)*(Dz2\(M2+Pm*kti^2*Id4));
Q2=(N2*BU2);
Bw=[0  Bb 0 Bt].';



%Matrice dérivée double suivant z
Dzz10=[[(s1^2-k1^2) -2*k1*s1]; [-2*k1*s1 -(s1^2-k1^2)]];
Dzz1=[[Dzz10 Z00];[Z00 Dzz10]];
Dzz20=[[(s2^2-k2^2) -2*k2*s2]; [-2*k2*s2 -(s2^2-k2^2)]];
Dzz2=[[Dzz20 Z00];[Z00 Dzz20]];


% Mat= [ N1 N2 ; Q1 Q2];
% eig_mat= eig(Mat);
% eig_mat= sort(real(eig_mat));
% nb_eig_mat= find(eig_mat > 0);
% ratio= eig_mat(nb_eig_mat(end))/eig_mat(nb_eig_mat(1));
% req_precision= log10(ratio);
%%
% SOLUTIONS FOR U, B and J

% The general form of the signal is :
%A(t,z) = exp(s1*z)*(A(1)*cos(t+k1*z)-A(2)*sin(t+k1*z))+
%	+exp(-s1*z)*(A(3)*cos(t-k1*z)-A(4)*sin(t-k1*z))+
%	+exp(s2*z)*(A(5)*cos(t+k2*z)-A(6)*sin(t+k2*z))+
%	+exp(-s2*z)*(A(7)*cos(t-k2*z)-A(8)*sin(t-k2*z))

% Be carefull to the minus before each coefficient of the sines. This minus
% should not absorb by the matrix writing !!

% Coefficients in velocity field expression from system of all 8 boundary conditions

A= [[ N1 N2 ];[Q1 Q2]]\[ Uw ; Bw];
U=A.'; %(change variable name)

% Coefficients in magnetic field expression from system of all 8 boundary conditions
B=vpa(zeros(1,8));
B(1:4)=BU1*A(1:4);
B(5:8)=BU2*A(5:8);
 
% matrice rotationnel
Rot_z10=[[s1, -k1];[k1 s1]];
Rot_z1=[[Rot_z10 Z00];[Z00 -Rot_z10]];
Rot_z20=[[s2, -k2];[k2 s2]];
Rot_z2=[[Rot_z20, Z00];[Z00 -Rot_z20]];


% B_grad1= Rot_z1*B(1:4)';
% B_grad2= Rot_z2*B(5:8)';
% B_grad= [B_grad1 ; B_grad2]';

% % Coefficient for electric potential field
% Dt=vpa(zeros(4));
% Dt(1,2)=-epsilon;
% Dt(2,1)=-epsilon;
% Dt(3,4)=-epsilon;
% Dt(4,3)=-epsilon;
% VU1=(Dz1\Dt-Dt/(Dz1))*BU1;
% VU2=(Dz2\Dt-Dt/(Dz2))*BU2;
% V=vpa(zeros(1,8));
% V(1:4)=VU1*A(1:4);
% V(5:8)=VU2*A(5:8);

% % Coefficient for electric field
% E=vpa(zeros(1,8));
% EU1=(Dz1\Dt)*BU1;
% EU2=(Dz2\Dt)*BU2;
% E(1:4)=(-Dz1\Dt)*B(1:4)';%EU1*A(1:4);
% E(5:8)=(-Dz2\Dt)*B(5:8)';%EU2*A(5:8);


% Coefficient for electric current
Jr=vpa(zeros(1,8));
Jr(1:4)= -Rot_z1*B(1:4)';
Jr(5:8)= -Rot_z2*B(5:8)';


% Coefficient for the vector potential
Ar=vpa(zeros(1,8));
Ar(1:4)=((Dzz1-kti^2*Id4)\Dz1) *B(1:4)'; %((Dzz1-kti^2*Id4)\Dz1) *B(1:4)';
Ar(5:8)= ((Dzz2-kti^2*Id4)\Dz2) *B(5:8)'; %((Dzz2-kti^2*Id4)\Dz2) *B(5:8)';

% Coefficient gradient du potentiel suivant er
 Vr=vpa(zeros(1,8));
Vr(1:4)= Jr(1:4)' + Dt*Ar(1:4)' - U(1:4)';
Vr(5:8)= Jr(5:8)' + Dt*Ar(5:8)' - U(5:8)';

 dtAr=vpa(zeros(1,8));
dtAr(1:4)= Dt*Ar(1:4)' ;
dtAr(5:8)= Dt*Ar(5:8)' ;



% Coefficient for the electric potential
Er=vpa(zeros(1,8));
Er(1:4)= Jr(1:4) - U(1:4);
Er(5:8)= Jr(5:8) - U(5:8);
% 

% % Coefficient for the potential
% Phi=vpa(zeros(1,8));
% phi(1:4)= Vr(1:4);
% phi(5:8)= Vr(5:8);


end
