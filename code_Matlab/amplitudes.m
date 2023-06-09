function [env,phase,z,Abot,Atop,Amin,zmin,Amax,zmax]=...
    amplitudes(A_num,s1_num,k1_num,s2_num,k2_num,h,npts,precision)

%A Calculates, bottom (z=0), top(z=h), max and  min amplitudes of oscillations for quantity A describing Alfven wave of the form
%A(t,z) = exp(s1*z)*(A(1)*cos(t+k1*z)-A(2)*sin(t+k1*z))+
%	+exp(-s1*z)*(A(3)*cos(t-k1*z)-A(4)*sin(t-k1*z))+
%	+exp(s2*z)*(A(5)*cos(t+k2*z)-A(6)*sin(t+k2*z))+
%	+exp(-s2*z)*(A(7)*cos(t-k2*z)-A(8)*sin(t-k2*z))+
% Typcally A is velocity, induced magnetic field or electric potential
% Alban Potherat 28/05/2018

% Parameters
%npts=100;

syms A s1 k1 s2 k2
if ~exist('precision','var')
    precision= 16;
end
digits(precision)

A = vpa(A_num);
s1 = vpa(s1_num);
k1 = vpa(k1_num);
s2 = vpa(s2_num);
k2 = vpa(k2_num);


% envelope and phase
z=0:h/(npts-1):h;

f=cos(k1*z).*(A(1)*exp(s1*z)+A(3)*exp(-s1*z))...
+sin(k1*z).*(-A(2)*exp(s1*z)+A(4)*exp(-s1*z))...
+cos(k2*z).*(A(5)*exp(s2*z)+A(7)*exp(-s2*z))...
+sin(k2*z).*(-A(6)*exp(s2*z)+A(8)*exp(-s2*z));

g=cos(k1*z).*(-A(2)*exp(s1*z)-A(4)*exp(-s1*z))...
+sin(k1*z).*(-A(1)*exp(s1*z)+A(3)*exp(-s1*z))...
+cos(k2*z).*(-A(6)*exp(s2*z)-A(8)*exp(-s2*z))...
+sin(k2*z).*(-A(5)*exp(s2*z)+A(7)*exp(-s2*z));


env=(f.^2+g.^2).^(1/2);
phase=atan(-g./f);


% Max and Min amplitudes
[Amin, imin]=min(env);
zmin=z(imin);
[Amax, imax]=max(env);
zmax=z(imax);

% Amplitudes at top and bottom
Abot=env(1);
Atop=env(end);
