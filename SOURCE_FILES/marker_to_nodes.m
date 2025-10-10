% UPDATING MARKER DENSITY !!!
% ind=(Im==9 | Im==10) & Pm>1e9 & Tm>500+273;
% rhom0(ind)=3400;

if Temperature == 1
    if TPdep_dens == 1
        rhom    = rhom0     .* (1 + tbm.*(Pm-1e5)) .* (1 - tam.*(Tm-298));
        rhocpm  = rhocpm0   .* (1 + tbm.*(Pm-1e5)) .* (1 - tam.*(Tm-298));
    else
        rhom    = rhom0;
        rhocpm  = rhocpm0;
    end
else
    rhom = rhom0;
end

% define indexes
j = fix(xm/dx+0.5)+1;
i = fix(ym/dy+0.5)+1;

ind= j<1;
j(ind)=1;
ind= j>nx;
j(ind)=nx;
ind= i<1;
i(ind)=1;
ind= i>ny;
i(ind)=ny;

% updating eta_s (i,j)
ind_move = (j-1).*(ny+1)+i;
wtijm              =   (1-abs(xm-x_Vx(j))/dx).*(1-abs(ym-y_Vy(i))/dy);
ind=Im==-1;
wtijm(ind)=0;
for fil=1:length(wtijm)
    eta_s(ind_move(fil))    =   eta_s(ind_move(fil))+(eta_effm(fil).*wtijm(fil));
    Toxy(ind_move(fil))     =   Toxy(ind_move(fil))+(Txym(fil).*wtijm(fil));
    mu_s(ind_move(fil))     =   mu_s(ind_move(fil))+(wtijm(fil)./mum(fil));
    rho_s(ind_move(fil))    =   rho_s(ind_move(fil))+(rhom(fil).*wtijm(fil));
    w(ind_move(fil))        =   w(ind_move(fil))+(wm(fil).*wtijm(fil));
    strain(ind_move(fil))   =   strain(ind_move(fil))+(strainm(fil).*wtijm(fil));
    strainv(ind_move(fil))   =   strainv(ind_move(fil))+(strainvm(fil).*wtijm(fil));
    work(ind_move(fil))     =   work(ind_move(fil))+(workm(fil).*wtijm(fil));
    lambda_s(ind_move(fil)) =   lambda_s(ind_move(fil))+(lambdam(fil).*wtijm(fil));
    if Diffusion_creep == 1
        straind(ind_move(fil))   =   straind(ind_move(fil))+(straindm(fil).*wtijm(fil));
        def_mode(ind_move(fil))  =   def_mode(ind_move(fil))+(def_modem(fil).*wtijm(fil));
        D_grain(ind_move(fil))   =   D_grain(ind_move(fil))+(dm(fil).*wtijm(fil));
        Dd_grain(ind_move(fil))  =   Dd_grain(ind_move(fil))+(ddm(fil).*wtijm(fil));
    end
    wts(ind_move(fil))      =   wts(ind_move(fil))+(wtijm(fil));
end

% define indexes
j = fix(xm/dx)+1;
i = fix(ym/dy)+1;

ind= j<1;
j(ind)=1;
ind= j>nx;
j(ind)=nx;
ind= i<1;
i(ind)=1;
ind= i>ny;
i(ind)=ny;

% updating eta_p(i+1,j+1) and rho_p(i+1,j+1)
ind_move = (j)*(ny+1)+i+1;
wtijm           =   (1-abs((xm-x_P(j+1))/dx)).*(1-abs((ym-y_P(i+1))/dy));
ind=Im==-1;
wtijm(ind)=0;
for fil=1:length(wtijm)
    eta_p(ind_move(fil))  =   eta_p(ind_move(fil))+eta_effm(fil)*wtijm(fil);
    Toxx(ind_move(fil))   =   Toxx(ind_move(fil))+Txxm(fil)*wtijm(fil);
    mu_p(ind_move(fil))   =   mu_p(ind_move(fil))+wtijm(fil)/mum(fil);
    rho_p(ind_move(fil))  =   rho_p(ind_move(fil))+rhom(fil)*wtijm(fil);
    wtp(ind_move(fil))    =   wtp(ind_move(fil)) + wtijm(fil);
end

if Temperature == 1
    if ra_heat == 1
        for fil=1:length(wtijm)
            Hr(ind_move(fil))     =   Hr(ind_move(fil)) + hrm(fil)*wtijm(fil);
        end
    end
    if adiab_heat == 1 || TPdep_dens == 1
        for fil=1:length(wtijm)
            Ta(ind_move(fil))     =   Ta(ind_move(fil)) + tam(fil)*wtijm(fil);
        end
    end
    for fil=1:length(wtijm)
        Temp(ind_move(fil))   =   Temp(ind_move(fil))+Tm(fil)*rhocpm(fil)*wtijm(fil);
        rhocpwt(ind_move(fil))=   rhocpwt(ind_move(fil)) + rhocpm(fil)*wtijm(fil);
        cp_p(ind_move(fil))   =   cp_p(ind_move(fil))+cpm(fil)*wtijm(fil);
        rhocp_p(ind_move(fil))=   rhocp_p(ind_move(fil))+rhocpm(fil)*wtijm(fil);
    end
end

% define indexes
j = fix(xm/dx)+1;
i = fix(ym/dy+0.5)+1;

ind= j<1;
j(ind)=1;
ind= j>nx;
j(ind)=nx;
ind= i<1;
i(ind)=1;
ind= i>ny;
i(ind)=ny;

% updating rho, k_vy (i,j+1)
ind_move = (j)*(ny+1)+i;
wtijm           =   (1-abs(xm-x_Vy(j+1))/dx).*(1-abs(ym-y_Vy(i))/dy);
ind=Im==-1;
wtijm(ind)=0;
for fil=1:length(wtijm)
    rho_vy(ind_move(fil))    =   rho_vy(ind_move(fil))+rhom(fil)*wtijm(fil);
    wtvy(ind_move(fil))      =   wtvy(ind_move(fil)) + wtijm(fil);
end
if Temperature == 1
    for fil=1:length(wtijm)
        k_vy(ind_move(fil))      =   k_vy(ind_move(fil))+km(fil)*wtijm(fil);
    end
end

% define indexes
j = fix(xm/dx+0.5)+1;
i = fix(ym/dy)+1;

ind= j<1;
j(ind)=1;
ind= j>nx;
j(ind)=nx;
ind= i<1;
i(ind)=1;
ind= i>ny;
i(ind)=ny;

% update k_vx (i+1,j)
ind_move = (j-1)*(ny+1)+i+1;
wtijm           =   (1-abs(ym-y_Vx(i+1))/dy).*(1-abs(xm-x_Vx(j))/dx);
ind=Im==-1;
wtijm(ind)=0;
for fil=1:length(wtijm)
    rho_vx(ind_move(fil))=   rho_vx(ind_move(fil))+rhom(fil)*wtijm(fil);
    wtvx(ind_move(fil))  =   wtvx(ind_move(fil)) + wtijm(fil);
end
if Temperature == 1
    for fil=1:length(wtijm)
        k_vx(ind_move(fil))  =   k_vx(ind_move(fil))+km(fil)*wtijm(fil);
    end
end
clear wtijm

% Divide by weights (i,j)
ind=wts>0;
eta_s(ind) = eta_s(ind)./wts(ind);
Toxy(ind)  = Toxy(ind)./wts(ind);
mu_s(ind)  = 1./(mu_s(ind)./wts(ind));
w(ind)     = w(ind)./wts(ind);
strain(ind)= strain(ind)./wts(ind);
strainv(ind)= strainv(ind)./wts(ind);
work(ind)  = work(ind)./wts(ind);
lambda_s(ind)= lambda_s(ind)./wts(ind);
rho_s(ind) = rho_s(ind)./wts(ind);
if Diffusion_creep == 1
    straind(ind)= straind(ind)./wts(ind);
    def_mode(ind)= def_mode(ind)./wts(ind);
    D_grain(ind)= D_grain(ind)./wts(ind);
    Dd_grain(ind)= Dd_grain(ind)./wts(ind);
end
    
ind=wts==0;
eta_s(ind) = eta_s_old(ind);
Toxy(ind)  = Toxy_old(ind);
mu_s(ind)  = mu_s_old(ind);
lambda_s(ind)= lambda_s_old(ind);
work(ind)  = work_old(ind);
w(ind)     = w_old(ind);
rho_s(ind) = rho_s_old(ind);
strain(ind)= strain_old(ind);
strainv(ind)= strainv_old(ind);
if Diffusion_creep == 1
    straind(ind)= straind_old(ind);
    def_mode(ind)= def_mode_old(ind);
    D_grain(ind)= D_grain_old(ind);
    Dd_grain(ind)= Dd_grain_old(ind);
end

ind=wtp>0;
eta_p(ind)  = eta_p(ind)./wtp(ind);
Toxx(ind)   = Toxx(ind)./wtp(ind);
mu_p(ind)   = 1./(mu_p(ind)./wtp(ind));
rho_p(ind)  = rho_p(ind)./wtp(ind);
if Temperature == 1
    if ra_heat == 1
        Hr(ind)     = Hr(ind)./wtp(ind);
    end
    if adiab_heat == 1 || TPdep_dens == 1
        Ta(ind)     = Ta(ind)./wtp(ind);
    end
    Temp(ind)   = Temp(ind)./rhocpwt(ind);
    cp_p(ind)   = cp_p(ind)./wtp(ind);
    rhocp_p(ind)= rhocp_p(ind)./wtp(ind);
end
ind=wtp==0;
eta_p(ind)  = eta_p_old(ind);
Toxx(ind)   = Toxx_old(ind);
mu_p(ind)   = mu_p_old(ind);
rho_p(ind)  = rho_p_old(ind);
if Temperature == 1
    if ra_heat == 1
        Hr(ind)     = Hr_old(ind);
    end
    if adiab_heat == 1 || TPdep_dens == 1
        Ta(ind)     = Ta_old(ind);
    end
    Temp(ind)   = Temp_old(ind);
    cp_p(ind)   = cp_p_old(ind);
    rhocp_p(ind)= rhocp_p_old(ind);
end

ind=wtvy>0;
rho_vy(ind) = rho_vy(ind)./wtvy(ind);
if Temperature == 1
    k_vy(ind)   = k_vy(ind)./wtvy(ind);
end
ind=wtvy==0;
rho_vy(ind) = rho_vy_old(ind);
if Temperature == 1
    k_vy(ind)   = k_vy_old(ind);
end
ind=wtvx>0;
rho_vx(ind) = rho_vx(ind)./wtvx(ind);
if Temperature == 1
    k_vx(ind)   = k_vx(ind)./wtvx(ind);
end
ind=wtvy==0;
rho_vx(ind) = rho_vx_old(ind);
if Temperature == 1
    k_vx(ind)   = k_vx_old(ind);
end

ddm(:) = 0;