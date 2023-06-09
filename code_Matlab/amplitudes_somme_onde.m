function [env,phase,z]=amplitudes_somme_onde(A_num,s1_num,...
    k1_num,s2_num,k2_num,h,npts,precision,pulsation_inf,phase_inf,...
    pulsation_sup,phase_sup)

%A Calculates, bottom (z=0), top(z=h), max and  min amplitudes of oscillations for quantity A describing Alfven wave of the form
%A(t,z) = exp(s1*z)*(A(1)*cos(t+k1*z)-A(2)*sin(t+k1*z))+
%	+exp(-s1*z)*(A(3)*cos(t-k1*z)-A(4)*sin(t-k1*z))+
%	+exp(s2*z)*(A(5)*cos(t+k2*z)-A(6)*sin(t+k2*z))+
%	+exp(-s2*z)*(A(7)*cos(t-k2*z)-A(8)*sin(t-k2*z))+
%
% Enveloppe calculée suivant la forme A(t,z)= f(z)*cos(wt) + g(z)*sin(wt)
% dephasage_num est le déphasage du forçage magnétique supérieur sur le
% forçage magnétique inférieur

% Le paramètre "type_freq" renseigne si, lors d'un forçage double, la
% fréquence supérieure est égale ou non à la fréquence inférieure. Lorsque
% les fréquences sont les mêmes renseigner : 0, 'egales'
% Lorsque les frequences sont différentes rensiegner : 1, "differentes"

% Les n/2 premiers coefficients de A sont considérés comme les solutions du
% forçage inférieur et les n/2 derniers coeffcients comme les solutions du 
% forçage supérieur. Nénmoins la fonction 'amplitudes_somme_onde.m' est aussi
% concue pour calculer l'enveloppe du signal de l'onde d'Alfven pour un
% nombre de modes transversaux supérieur à 1 et ce même lors d'un forçage
% simple (plaque inférieure ou plaque supérieure). Ce cas de forçage simple
% est implicitement traité lorsque 'dephasage_num' et 'type_freq' ne sont
% pas spécifiés. Attention !! ce programme ne traite pas le cas d'un
% forçage simple avec un mode transversal. Pour ce dernier, se référer à la
% fonction 'amplitudes.m'

% Typcally A is velocity, induced magnetic field or electric potential
% Samy Lalloz 21/01/2021

% Parameters
%npts=100;
% Definition des variables symboliques
    syms A s1 k1 s2 k2 phi
    digits(precision)

    A = vpa(A_num);
    s1 = vpa(s1_num);
    k1 = vpa(k1_num);
    s2 = vpa(s2_num);
    k2 = vpa(k2_num);
    z=0:h/(npts-1):h;
    
    if exist('pulsation_inf') && exist('pulsation_sup')
        maxi= max(pulsation_inf,pulsation_sup);
        mini= min(pulsation_inf,pulsation_sup);
    elseif exist('pulsation_sup') 
        maxi= pulsation_sup;
        mini= pulsation_sup;   
    elseif exist('pulsation_inf') 
        maxi= pulsation_inf;
        mini= pulsation_inf;
    else
        maxi= 1000;
        mini= 1000;
    end
        
        
    time= [-4*vpa(pi)/mini:0.01*2*vpa(pi)/maxi:4*vpa(pi)/mini]';
    L_time= length(time);
    nb_mode= length(A(:,1)); % et non 0.5 * nb_transverse mode
    
    if mod(nb_mode,2) ~= 0 || nb_mode== 2
    
        f= cos(k1.*z).*(A(:,1).*exp(s1.*z)+A(:,3).*exp(-s1.*z))...
       +sin(k1.*z).*(-A(:,2).*exp(s1.*z)+A(:,4).*exp(-s1.*z))...
       +cos(k2.*z).*(A(:,5).*exp(s2.*z)+A(:,7).*exp(-s2.*z))...
       +sin(k2.*z).*(-A(:,6).*exp(s2.*z)+A(:,8).*exp(-s2.*z));

        g= cos(k1.*z).*(-A(:,2).*exp(s1.*z)-A(:,4).*exp(-s1.*z))...
       +sin(k1.*z).*(-A(:,1).*exp(s1.*z)+A(:,3).*exp(-s1.*z))...
       +cos(k2.*z).*(-A(:,6).*exp(s2.*z)-A(:,8).*exp(-s2.*z))...
       +sin(k2.*z).*(-A(:,5).*exp(s2.*z)+A(:,7).*exp(-s2.*z));
    
        if nb_mode>1
            sf= sum(f);
            sg= sum(g);
        else
            sf= f;
            sg= g;
        end
        if  exist('phase_inf','var')         
            phi= vpa(phase_inf);
            sf_phi=sf;
            sg_phi=sg;            
            sf= sf_phi*cos(phi) + sg_phi*sin(phi);
            sg= -sf_phi*sin(phi) + sg_phi*cos(phi);
            clear sf_phi sg_phi
        end
        signal= sf.*cos(pulsation_inf*time)+sg.*sin(pulsation_inf*time);
        signal_ana= hilbert(double(signal));
        env= abs(signal_ana(floor(0.5*end)+1,:));      %
        phase= unwrap(angle(signal_ana(floor(0.5*end)+1,:)));   
        
%         env=(sf.^2+sg.^2).^(1/2);
%         phase=atan(-sg./sf);
        
    else 
        L= nb_mode/2;
        % envelope and phase
        % finf et ginf sont les coefficients pour les solutions de
        % l'oscillation de la plaque inférieure.
        finf= cos(k1(1:L).*z).*(A(1:L,1).*exp(s1(1:L).*z)+A(1:L,3).*exp(-s1(1:L).*z))...
       +sin(k1(1:L).*z).*(-A(1:L,2).*exp(s1(1:L).*z)+A(1:L,4).*exp(-s1(1:L).*z))...
       +cos(k2(1:L).*z).*(A(1:L,5).*exp(s2(1:L).*z)+A(1:L,7).*exp(-s2(1:L).*z))...
       +sin(k2(1:L).*z).*(-A(1:L,6).*exp(s2(1:L).*z)+A(1:L,8).*exp(-s2(1:L).*z));

        ginf= cos(k1(1:L).*z).*(-A(1:L,2).*exp(s1(1:L).*z)-A(1:L,4).*exp(-s1(1:L).*z))...
       +sin(k1(1:L).*z).*(-A(1:L,1).*exp(s1(1:L).*z)+A(1:L,3).*exp(-s1(1:L).*z))...
       +cos(k2(1:L).*z).*(-A(1:L,6).*exp(s2(1:L).*z)-A(1:L,8).*exp(-s2(1:L).*z))...
       +sin(k2(1:L).*z).*(-A(1:L,5).*exp(s2(1:L).*z)+A(1:L,7).*exp(-s2(1:L).*z));

        % fsup et gsup sont les coefficients des solutions de l'oscillation
        % avec la plaque supérieure
        fsup= cos(k1(L+1:end).*z).*(A(L+1:end,1).*exp(s1(L+1:end).*z)+A(L+1:end,3).*exp(-s1(L+1:end).*z))...
       +sin(k1(L+1:end).*z).*(-A(L+1:end,2).*exp(s1(L+1:end).*z)+A(L+1:end,4).*exp(-s1(L+1:end).*z))...
       +cos(k2(L+1:end).*z).*(A(L+1:end,5).*exp(s2(L+1:end).*z)+A(L+1:end,7).*exp(-s2(L+1:end).*z))...
       +sin(k2(L+1:end).*z).*(-A(L+1:end,6).*exp(s2(L+1:end).*z)+A(L+1:end,8).*exp(-s2(L+1:end).*z));

        gsup= cos(k1(L+1:end).*z).*(-A(L+1:end,2).*exp(s1(L+1:end).*z)-A(L+1:end,4).*exp(-s1(L+1:end).*z))...
       +sin(k1(L+1:end).*z).*(-A(L+1:end,1).*exp(s1(L+1:end).*z)+A(L+1:end,3).*exp(-s1(L+1:end).*z))...
       +cos(k2(L+1:end).*z).*(-A(L+1:end,6).*exp(s2(L+1:end).*z)-A(L+1:end,8).*exp(-s2(L+1:end).*z))...
       +sin(k2(L+1:end).*z).*(-A(L+1:end,5).*exp(s2(L+1:end).*z)+A(L+1:end,7).*exp(-s2(L+1:end).*z));  

        if  exist('phase_sup','var')         
            phi= vpa(phase_sup);
            fsup_phi=fsup;
            gsup_phi=gsup;            
            fsup= fsup_phi*cos(phi) + gsup_phi*sin(phi);
            gsup= -fsup_phi*sin(phi) + gsup_phi*cos(phi);
            clear fsup_phi gsup_phi
        end
        if L>1
            sfinf= sum(finf);
            sfsup= sum(fsup);
            sginf= sum(ginf);
            sgsup= sum(gsup);
        else
            sfinf= finf;
            sfsup= fsup;
            sginf= ginf;
            sgsup= gsup;
        end
        if exist('pulsation_inf','var') && exist('pulsation_sup','var')
            if  exist('phase_inf','var')         
                phi= vpa(phase_inf);
                sfinf_phi=sfinf;
                sginf_phi=sginf;            
                sfinf= sfinf_phi*cos(phi) + sginf_phi*sin(phi);
                sginf= - sfinf_phi*sin(phi) + sginf_phi*cos(phi);
                clear sfinf_phi sginf_phi
            end
            
            switch pulsation_sup
                case pulsation_inf
                    signal= (sfinf+sfsup).*cos(pulsation_inf*time)+(sginf+sgsup).*sin(pulsation_inf*time);
                otherwise
                    signal= sfinf.*cos(pulsation_inf*time)+sginf.*sin(pulsation_inf*time)+...
                        sfsup.*cos(pulsation_sup*time)+sgsup.*sin(pulsation_sup*time);
            end
        elseif exist('pulsation_inf','var')
            sf= sfinf + sfsup;
            sg= sginf + sgsup;
            if  exist('phase_inf','var')         
                phi= vpa(phase_inf);
                sf_phi=sfinf+sfsup;
                sg_phi=sginf+sgsup;            
                sf= sf_phi*cos(phi) + sg_phi*sin(phi);
                sg= - sf_phi*sin(phi) + sg_phi*cos(phi);
                clear sf_phi sg_phi
            end
            
            signal= sf.*cos(pulsation_inf*time)+sg.*sin(pulsation_inf*time);
        else
            sf= sfinf + sfsup;
            sg= sginf + qgsup;
            if  exist('phase_inf','var')         
                phi= vpa(phase_inf);
                sf_phi=sfinf+sfsup;
                sg_phi=sginf+sgsup;            
                sf= sf_phi*cos(phi) + sg_phi*sin(phi);
                sg= - sf_phi*sin(phi) + sg_phi*cos(phi);
                clear sf_phi sg_phi
            end
            signal= sf.*cos(2*pi*time)+sg.*sin(2*pi*time);
        end
        signal_ana= hilbert(double(signal));
        env= abs(signal_ana(floor(0.5*end)+1,:));   
        phase= unwrap(angle(signal_ana(floor(0.5*end)+1,:)));   
    end    
%     % détermination des coefficient f et g avec un déphasage de phi
%     % entre les oscillations des forçages inf et sup
%     % on a : 
%     % Asup_phi(z,t)  = fsup(z)*cos(wt+phi) + gsup(z)*sin(wt+phi)
%     %                = (fsup*cos(phi) +gsup*sin(phi))*cos(wt) +...
%     %                = (-fsup*sin(phi) +gsup*cos(phi))*sin(wt)
%     %                = fsup_phi*cos(wt) + gsup_phi*sin(wt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%                 
%     if ~exist('dephasage_num','var')
% 
%         % envelope and phase
%         
%         f=cos(k1.*z).*(A(:,1).*exp(s1.*z)+A(:,3).*exp(-s1.*z))...
%         +sin(k1.*z).*(-A(:,2).*exp(s1.*z)+A(:,4).*exp(-s1.*z))...
%         +cos(k2.*z).*(A(:,5).*exp(s2.*z)+A(:,7).*exp(-s2.*z))...
%         +sin(k2.*z).*(-A(:,6).*exp(s2.*z)+A(:,8).*exp(-s2.*z));
% 
%         g=cos(k1.*z).*(-A(:,2).*exp(s1.*z)-A(:,4).*exp(-s1.*z))...
%         +sin(k1.*z).*(-A(:,1).*exp(s1.*z)+A(:,3).*exp(-s1.*z))...
%         +cos(k2.*z).*(-A(:,6).*exp(s2.*z)-A(:,8).*exp(-s2.*z))...
%         +sin(k2.*z).*(-A(:,5).*exp(s2.*z)+A(:,7).*exp(-s2.*z));
% 
% 
%         sf= sum(f);
%         sg= sum(g);
% 
%         env=((sf.^2)+(sg.^2)).^0.5;
%         phase=-atan((sg)./(sf));
%         
%     elseif  exist('dephasage_num','var')  
%         
%         phi= vpa(dephasage_num);
%         L= 0.5*length(A(:,1));
% 
%         % envelope and phase
%         % finf et ginf sont les coefficients pour les solutions de
%         % l'oscillation de la plaque inférieure. Phase = 0 car la phase des
%         % solutions de l'oscillation inférieure est la phase de ref
%         
%         finf=  cos(k1(1:L).*z).*(A(1:L,1).*exp(s1(1:L).*z)+A(1:L,3).*exp(-s1(1:L).*z))...
%            +sin(k1(1:L).*z).*(-A(1:L,2).*exp(s1(1:L).*z)+A(1:L,4).*exp(-s1(1:L).*z))...
%            +cos(k2(1:L).*z).*(A(1:L,5).*exp(s2(1:L).*z)+A(1:L,7).*exp(-s2(1:L).*z))...
%            +sin(k2(1:L).*z).*(-A(1:L,6).*exp(s2(1:L).*z)+A(1:L,8).*exp(-s2(1:L).*z));
% 
%         ginf=  cos(k1(1:L).*z).*(-A(1:L,2).*exp(s1(1:L).*z)-A(1:L,4).*exp(-s1(1:L).*z))...
%            +sin(k1(1:L).*z).*(-A(1:L,1).*exp(s1(1:L).*z)+A(1:L,3).*exp(-s1(1:L).*z))...
%            +cos(k2(1:L).*z).*(-A(1:L,6).*exp(s2(1:L).*z)-A(1:L,8).*exp(-s2(1:L).*z))...
%            +sin(k2(1:L).*z).*(-A(1:L,5).*exp(s2(1:L).*z)+A(1:L,7).*exp(-s2(1:L).*z));
% 
%        
%        % fsup et gsup sont les coefficients des solutions de l'oscillation
%        % avec la plaque supérieure sans déphasage par rapport à
%        % l'oscillation de la plaque inférieure cad
%        % Asup(z,t) = fsup*cos(wt) + gsup*sin(wt)
%         fsup=  cos(k1(L+1:end).*z).*(A(L+1:end,1).*exp(s1(L+1:end).*z)+A(L+1:end,3).*exp(-s1(L+1:end).*z))...
%            +sin(k1(L+1:end).*z).*(-A(L+1:end,2).*exp(s1(L+1:end).*z)+A(L+1:end,4).*exp(-s1(L+1:end).*z))...
%            +cos(k2(L+1:end).*z).*(A(L+1:end,5).*exp(s2(L+1:end).*z)+A(L+1:end,7).*exp(-s2(L+1:end).*z))...
%            +sin(k2(L+1:end).*z).*(-A(L+1:end,6).*exp(s2(L+1:end).*z)+A(L+1:end,8).*exp(-s2(L+1:end).*z));
% 
%         gsup=  cos(k1(L+1:end).*z).*(-A(L+1:end,2).*exp(s1(L+1:end).*z)-A(L+1:end,4).*exp(-s1(L+1:end).*z))...
%            +sin(k1(L+1:end).*z).*(-A(L+1:end,1).*exp(s1(L+1:end).*z)+A(L+1:end,3).*exp(-s1(L+1:end).*z))...
%            +cos(k2(L+1:end).*z).*(-A(L+1:end,6).*exp(s2(L+1:end).*z)-A(L+1:end,8).*exp(-s2(L+1:end).*z))...
%            +sin(k2(L+1:end).*z).*(-A(L+1:end,5).*exp(s2(L+1:end).*z)+A(L+1:end,7).*exp(-s2(L+1:end).*z));
%        
%        % détermination des coefficient f et g avec un déphasage de phi
%        % entre les oscillations des forçages inf et sup
%        % on a : 
%        % Asup_phi(z,t)  = fsup(z)*cos(wt+phi) + gsup(z)*sin(wt+phi)
%        %                = (fsup*cos(phi) +gsup*sin(phi))*cos(wt) +...
%        %                = (-fsup*sin(phi) +gsup*cos(phi))*sin(wt)
%        %                = fsup_phi*cos(wt) + gsup_phi*sin(wt)
%        
%        fsup_phi= fsup*cos(phi) + gsup*sin(phi);
%        gsup_phi= -fsup*sin(phi) + gsup*cos(phi);
%        
%        sf= sum([finf ; fsup_phi]);
%        sg= sum([ginf ; gsup_phi]);
% 
%         env=((sf.^2)+(sg.^2)).^0.5;
%         phase=-atan((sg)./(sf));
%     end  
% elseif exist('type_freq','var')
%     switch type_freq
%         case (
% end     
% % Max and Min amplitudes
% %[Amin, imin]=min(env);
% %zmin=z(imin);
% %[Amax, imax]=max(env);
% %zmax=z(imax);
% 
% % Amplitudes at top and bottom
% %Abot=env(1);
% %Atop=env(end);
