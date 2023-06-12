function Hank_T_forcing= Hankel_T_type_forcing(type_forcing,kt_vect,r_ext,r0,precision)

% Function which return the Hankel transform based of J1 for a given
% axisymetric forcing type
% Two type are code:
% - linear forcing: forcing starting from 0 at r=0 to 2 at r= R (in order
% have a unit integral
% - forcing derived from a gaussien distribution: [(1-exp(-r^2))/r]  
% For this forcing, the caracteristic width of the gaussian needs to be
% specified (i.e. its standard deviation)

arguments
    type_forcing  char {mustBeMember(type_forcing,{'linear','Gaussian_deriv'})}
    kt_vect (1,:) {double, sym}     
    r_ext = ''
    r0= ''
    precision(1,1) {double}= 0
end
narginchk(4,5)
if precision ==0
    precision= 10*32; 
else
end% Setup default precision to 40 decimal digits (quadruple).
digits(precision);  

switch type_forcing
    case 'linear'
        if ~exist('r_ext')
            error('external radius not specified')
        else
        Hank_T_forcing= r_ext^2./kt_vect.*besselj(2,kt_vect*r_ext);
        end
        
    case 'Gaussian_deriv'
        if ~exist('r0')
            error('Caracteristic length of the gaussian not specified') 
        else
        Hank_T_forcing= (2*vpa(pi))^-1*(1./kt_vect- 0.5*r0*sqrt(vpa(pi))...
            .*exp(-r0^2*kt_vect.^2/8).*besseli(0.5,r0^2*kt_vect.^2/8));
        end
end
