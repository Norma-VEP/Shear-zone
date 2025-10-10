% Calculating plastic strain
strainm(ind_plast) = strainm(ind_plast) + E2ndm(ind_plast)*dt;
% Calculating viscous strain
if Diffusion_creep == 0
    strainvm(~ind_plast) = strainvm(~ind_plast) + (E2ndm(~ind_plast).*dt);  % strain dislocation creep
elseif Diffusion_creep == 1 && timestep>1
    strainvm(~ind_plast) = strainvm(~ind_plast) + (E2ndm(~ind_plast).*dt.*(1-def_modem(~ind_plast)));  % strain dislocation creep
    straindm(~ind_plast) = straindm(~ind_plast) + (E2ndm(~ind_plast).*dt.*def_modem(~ind_plast));       % strain diffusion creep

    % Grain size evolution
    
%=======================================================================
    % for Quartz
    ind=find(ind_plast==0 & Im>2 & Im<=10);
    ff = 5521e6.*exp(-(31.28e3-2.009e-5.*Pm)./(RK.*Tm)); % Fugacity after Shinevar et al. (2015)
    p = 3; % from Tokle
    Qg = 128e3; % from Tokle
    K_g = (0.19.*(ff).^1.33).*10.^(-p*6); % from Tokle
    lam = 0.015; % from Tokle
    GBE=270e-3; % Grain Boundary Energy after Hiraga et al. (2002)
    
    % equilibrium after Austin and Evans (2007)
%     if timestep<20
        dm(ind)  = (   (   (K_g(ind).*exp(-(Qg./RK./Tm(ind)))*(p.^-1).* 3.1416.*GBE) ./ (lam.*2.*T2ndm(ind).*(E2ndm(ind).*(1-def_modem(ind))))     ).^(1./(1+p)));
%     else
        %     % grain size over time after Austin and Evans (2007)
%         dm(ind) = dm(ind) + ((-((lam.*2.*T2ndm(ind).*(E2ndm(ind).*(1-def_modem(ind)))).*dm(ind).^2) ./ (3.1416.*GBE)) + K_g(ind) .* exp(-(Qg./RK./Tm(ind))) .* (p.^-1).* dm(ind).^(1-p))*dt;    
%     end
%=======================================================================
    % for Anorthite
    ind=find(ind_plast==0 & Im>=11 & Im<=12);
    p = 2.6; % from Dresen et al. (1996)
    Qg = 365e3; % from Dresen et al. (1996)
    K_g = 2.59e-4; % from Dresen et al. (1996)
    lam = 0.1; % from Austin % Evans (2007)
    GBE=1; % from Austin % Evans (2007)
    
    % equilibrium after Austin and Evans (2007)
%     if timestep<20
        dm(ind)  = (   (   (K_g.*exp(-(Qg./RK./Tm(ind)))*(p.^-1).* 3.1416.*GBE) ./ (lam.*2.*T2ndm(ind).*(E2ndm(ind).*(1-def_modem(ind))))     ).^(1./(1+p)));
%     else
        %     % grain size over time after Austin and Evans (2007)
%         dm(ind) = dm(ind) + ((-((lam.*2.*T2ndm(ind).*(E2ndm(ind).*(1-def_modem(ind)))).*dm(ind).^2) ./ (3.1416.*GBE)) + K_g .* exp(-(Qg./RK./Tm(ind))) .* (p.^-1).* dm(ind).^(1-p))*dt;
%     end
%=======================================================================
    % for Olivine
    ind=find(ind_plast==0 & Im>=15 & Im<=17);
%     p = 2; % from Austin & Evans (2007)
    p = 3.2; % from Speciale
%     Qg = 520e3; % from Austin & Evans (2007)
    Qg = 620e3+Pm.*5e-6; % from Speciale
    K_g = 1800; % from Speciale
%     K_g = 70000; % from VanDerWal et al. (1990)
    lam = 0.1; % Austin & Evans (2007)
    GBE=1.4; % Grain Boundary Energy after Duyster and St?ckhert (2001)
    
    % equilibrium after Austin and Evans (2007)
%     if timestep<20
        dm(ind)  = (   (   (K_g.*exp(-(Qg(ind)./RK./Tm(ind)))*(p.^-1).* 3.1416.*GBE) ./ (lam.*2.*T2ndm(ind).*(E2ndm(ind).*(1-def_modem(ind))))     ).^(1./(1+p)));
%     else
            % grain size over time after Austin and Evans (2007)
        dm(ind) = dm(ind) + ((-((lam.*2.*T2ndm(ind).*(E2ndm(ind).*(1-def_modem(ind)))).*dm(ind).^2) ./ (3.1416.*GBE)) + K_g .* exp(-(Qg(ind)./RK./Tm(ind))) .* (p.^-1).* dm(ind).^(1-p))*dt;
    end
    
%     ind=Im==6; % Biotite
%     dm(ind)=0;

%     ind=find(ind_plast==0 & Im>2 & dm>10);
%     dm(ind)=10;
%     ind=find(ind_plast==0 & Im>2 & dm<1e-8);
%     dm(ind)=1e-8;
end




% Strain weakeing
ind = strainm > w_1;
phim(ind)    = phi(Im(ind))-(strainm(ind)-w_1)/(w_2-w_1).*(phi(Im(ind))-phi_w(Im(ind)));
Cm(ind)      = C(Im(ind))-(strainm(ind)-w_1)/(w_2-w_1).*(C(Im(ind))-C_w(Im(ind)));

ind = strainm > w_2;
phim(ind)    = phi_w(Im(ind));
Cm(ind)      = C_w(Im(ind));