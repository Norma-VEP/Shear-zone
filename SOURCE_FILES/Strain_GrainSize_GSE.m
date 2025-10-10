% if timestep <= 30
%     ind=find(Im>14 & ym<160e3);
%     dm(ind)=5*10.^(-3.+(ym(ind)-43e3)./117e3.*1.301);
%     ind=find(Im>14 & ym>=160e3);
%     dm(ind)=1e-1;
% end

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
    if length(find(ind_plast==1 & Im>1 & Im<=3)) < 1e7
        ind=find(ind_plast==0 & Im>1 & Im<=3);
        ff = 5521e6.*exp(-(31.28e3-2.009e-5.*Pm)./(RK.*Tm)); % Fugacity after Shinevar et al. (2015)
        p = 3; % from Tokle
        Qg = 134e3; % from Tokle
        K_g = (0.261.*(ff./1e6).^1.38).*10.^(-p*6); % from Tokle
        lam = 1e-4 * 10.^(-sind(((strainvm+straindm)+3/2)*180/3)+1)*1.5;
        kind = find(strainvm+straindm>=3);
        lam(kind)=0.015;

        GBE=1; % Grain Boundary Energy after Tokle and Hirth (2021)
        
        % equilibrium after Austin and Evans (2007)
        %     if timestep==12
        %         dm(ind)  = (   (   (K_g(ind).*exp(-(Qg./RK./Tm))*(p.^-1).* 3.1416.*GBE) ./ (lam.*T2ndm(ind).*(E2ndm(ind).*(1-def_modem(ind))))     ).^(1./(1+p)));
        %     else
        %     % grain size over time after Austin and Evans (2007)
        ddm(ind) = ((-((lam(ind).*T2ndm(ind).*(E2ndm(ind).*(1-def_modem(ind)))).*dm(ind).^2) ./ (3.1416.*GBE)) + K_g(ind) .* exp(-(Qg./RK./Tm)) .* (p.^-1).* dm(ind).^(1-p));
        
        
        %     end
        %     ind=find(ind_plast==1 & Im>1 & Im<=3);
        %     dm(ind)  = (   (   (K_g(ind).*exp(-(Qg./RK./Tm))*(p.^-1).* 3.1416.*GBE) ./ (lam.*T2ndm(ind).*(E2ndm(ind).*(1-0)))     ).^(1./(1+p)));
    end
    %=======================================================================
    % for Anorthite
    if length(find(ind_plast==1 & Im>=4 & Im<=5)) <= 4.0208e+04
        ind=find(ind_plast==0 & Im>=4 & Im<=5);
        p = 2.6; % from Dresen et al. (1996)
        Qg = 365e3; % from Dresen et al. (1996)
        K_g = 2.59e-4; %.*ff./0.1e6; % from Dresen et al. (1996)
        lam = 1e-4 * 10.^(-sind(((strainvm+straindm)+3/2)*180/3)+1)*1.0;
        kind = find(strainvm+straindm>=3);
        lam(kind)=0.010;

        GBE=1; % from Austin % Evans (2007)
        
        % equilibrium after Austin and Evans (2007)
        %     if timestep==12
        %         dm(ind)  = (   (   (K_g.*exp(-(Qg./RK./Tm))*(p.^-1).* 3.1416.*GBE) ./ (lam.*T2ndm(ind).*(E2ndm(ind).*(1-def_modem(ind))))     ).^(1./(1+p)));
        %     else
        %     % grain size over time after Austin and Evans (2007)
        ddm(ind) = ((-((lam(ind).*T2ndm(ind).*(E2ndm(ind).*(1-def_modem(ind)))).*dm(ind).^2) ./ (3.1416.*GBE)) + K_g .* exp(-(Qg./RK./Tm)) .* (p.^-1).* dm(ind).^(1-p));
        
        %         dm(ind) = dm(ind) + ddm(ind)*dt;
        
        %    end
        %     ind=find(ind_plast==1 & Im>=4 & Im<=5) ;
        %     dm(ind)  = (   (   (K_g.*exp(-(Qg./RK./Tm))*(p.^-1).* 3.1416.*GBE) ./ (lam.*T2ndm(ind).*(E2ndm(ind).*(1-0)))     ).^(1./(1+p)));
        
        
    end
%=======================================================================
%     % for Olivine
%     ind=find(ind_plast==0 & Im>=15 & Im<=17);
% %     p = 2; % from Austin & Evans (2007)
%     p = 3.2; % from Speciale
% %     Qg = 520e3; % from Austin & Evans (2007)
%     Qg = 620e3+Pm.*5e-6; % from Speciale
%     K_g = 3e-6*OH_const; % from Speciale
% %     K_g = 70000; % from VanDerWal et al. (1990)
%     lam = 0.01; % Austin & Evans (2007)
%     GBE=1.4; % Grain Boundary Energy after Duyster and St?ckhert (2001)
%     
%     % equilibrium after Austin and Evans (2007)
% %     if timestep<20
% %         dm(ind)  = (   (   (K_g.*exp(-(Qg(ind)./RK./Tm(ind)))*(p.^-1).* 3.1416.*GBE) ./ (lam.*2.*T2ndm(ind).*(E2ndm(ind).*(1-def_modem(ind))))     ).^(1./(1+p)));
% %     else
%         %     % grain size over time after Austin and Evans (2007)
%         ddm(ind) = ((-((lam.*T2ndm(ind).*(E2ndm(ind).*(1-def_modem(ind)))).*dm(ind).^2) ./ (3.1416.*GBE)) + K_g .* exp(-(Qg(ind)./RK./Tm(ind))) .* (p.^-1).* dm(ind).^(1-p));
%         dm(ind) = dm(ind) + ddm(ind)*dt;
% %     end
%     ind=find((ind_plast==1 | dm>1 | dm<1e-6) & Im>=15 & Im<=17);
%     dm(ind)  = (   (   (K_g.*exp(-(Qg(ind)./RK./Tm(ind)))*(p.^-1).* 3.1416.*GBE) ./ (lam.*T2ndm(ind).*(E2ndm(ind).*(1-0)))     ).^(1./(1+p)));
% 
    
%     ind=Im==6; % Biotite
%     dm(ind)=0;

    ind=find(ind_plast~=0);
    ddm(ind)=0;

    ind=find(ddm~=0 & ddm>1e-13);
    ddm(ind)=1e-13;
    ind=find(ddm~=0 & ddm<-1e-12);
    ddm(ind)=-1e-12;

    dm = dm + ddm*dt;

    ind=find(ind_plast==0 & Im>1 & dm>10);
    dm(ind)=10;
    ind=find(ind_plast==0 & Im>1 & dm<1e-8);
    dm(ind)=1e-8;
end




% Strain weakeing
ind = strainm > w_1;
phim(ind)    = phi(Im(ind))-(strainm(ind)-w_1)/(w_2-w_1).*(phi(Im(ind))-phi_w(Im(ind)));
Cm(ind)      = C(Im(ind))-(strainm(ind)-w_1)/(w_2-w_1).*(C(Im(ind))-C_w(Im(ind)));

ind = strainm > w_2;
phim(ind)    = phi_w(Im(ind));
Cm(ind)      = C_w(Im(ind));