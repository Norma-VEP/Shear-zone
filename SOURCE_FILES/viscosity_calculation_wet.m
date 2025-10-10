
%==============================================================
% Initial effective viscosity
eta_cur         =   etam;
eta_diffm       =   etam;

% Power-law viscosity
if Powerlaw==1
    co=(Qm./(RK.*Tm));
    ind=Im==5;
    co(ind)=((Qm(ind)+Pm(ind).*3.8e-5)./(RK.*Tm));
    ind=co>150;
    co(ind)=150;
    
    ff = 5521e6.*exp(-(31.28e3-2.009e-5.*Pm)./(8.314.*Tm));
    ind=Im==3;
    etam(ind)=1./(1.75*10^-12*ff(ind)./1e6*10^(-6*4));
    ind=Im==2;
    etam(ind)=1.9953e21.*50e6./ff(ind);
    ind=Im==5;
    etam(ind)=1./(10^0.2*ff(ind)./1e6*10^(-6*3));

%     ind=Im>=15;
%     etam(ind)=1./(90*10^(-6*3.5)*(OH_const./1e6).^1.2);    

    %         eta_ddd     =   ((etam.*exp(co).*E2ndm).^(1./nm))./2./E2ndm;
%     eta_cur     =   2.^((1-nm)./nm).*3.^(-(1+nm)./(2.*nm)).*etam.^(1./nm).*E2ndm.^((1-nm)./nm).*exp(co./nm);    
    eta_cur     =   0.5.*etam.^(1./nm).*E2ndm.^((1-nm)./nm).*exp(co./nm);
    
    if Diffusion_creep == 1 && timestep>1
        co=(Qdm./(RK.*Tm));
        ind=Im==5;
        co(ind)=((Qdm(ind)+Pm(ind).*3.8e-5)./(RK.*Tm));
        ind=co>150;
        co(ind)=150;

%         ind=Im>=15;
%         etadm(ind)=1./(1e6*(OH_const./1e6).^1.0*10^(-6)*10^(-3*6));

        % diffusion creep of matrix
        ind=Im==2|Im==3;
        etadm(ind)=1.5849e18.*50e6./ff(ind);
        ind=Im==5;
        etadm(ind)=1./(10^-0.7*ff(ind)./1e6*10^(-6*3)*10^(-6));

        
        ind=etadm~=0;
        def_modem(1,:) = 0;

        eta_diffm(ind) =   0.5.*etadm(ind).*dm(ind).^-mm(ind).*exp(co(ind));
        def_modem(ind) = 1./eta_diffm(ind) ./ (1./eta_cur(ind) + 1./eta_diffm(ind)); % 0 = dislocation creep | 1 = diffusion creep
        eta_cur(ind) = 1./(1./eta_cur(ind)+1./eta_diffm(ind));
        
    end
    
    ind=nm==1;
    eta_cur(ind)=etam(ind);
end

ind = eta_cur > eta_max;
eta_cur(ind) = eta_max;

% Brittle viscosity
Txxm_new        =   Txxm.*(1-(mum*Ddt)./(mum*Ddt+eta_cur))+2*eta_cur.*Exxm.*(mum*Ddt)./(mum*Ddt+eta_cur);
Txym_new        =   Txym.*(1-(mum*Ddt)./(mum*Ddt+eta_cur))+2*eta_cur.*Exym.*(mum*Ddt)./(mum*Ddt+eta_cur);

T2ndm_new       =   (Txxm_new.^2 + Txym_new.^2).^0.5;
T2ndm           =   (Txxm.^2 + Txym.^2).^0.5;

Yield       =   sind(phim).*(1 - lambdam) .* Pm + Cm .* cosd(phim);
ind = Yield < 0;
Yield(ind) = 0;

ind_plast = Yield < T2ndm_new & Im>0;
eta_cur(ind_plast)    =   0.5*Yield(ind_plast)./E2ndm(ind_plast);
Txxm(ind_plast)       =   Txxm(ind_plast).*Yield(ind_plast)./(T2ndm(ind_plast));
Txym(ind_plast)       =   Txym(ind_plast).*Yield(ind_plast)./(T2ndm(ind_plast));

if Diffusion_creep == 1
    def_modem(ind_plast) = -1;
end

ind = eta_cur > eta_max;
eta_cur(ind) = eta_max;
ind = eta_cur < eta_min;
eta_cur(ind) = eta_min;

eta_effm = eta_cur;
%==========================================================
