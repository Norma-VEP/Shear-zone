
% Save old values
eta_s_old   = eta_s;
Toxy_old    = Toxy;
mu_s_old    = mu_s;
w_old       = w;
rho_s_old   = rho_s;
eta_p_old   = eta_p;
Toxx_old    = Toxx;
mu_p_old    = mu_p;
rho_p_old   = rho_p;
rho_vy_old  = rho_vy;
rho_vx_old  = rho_vx;
% if timestep == 1
%     rho_vy_old(:)=2500;
%     rho_vx_old(:)=2500;
% end
strain_old  = strain;
strainv_old  = strainv;
work_old  = work;
lambda_s_old= lambda_s;
if Temperature == 1
    if ra_heat ==1
        Hr_old      = Hr;
    end
    if adiab_heat == 1 || TPdep_dens == 1
        Ta_old      = Ta;
    end
    k_vx_old    = k_vx;
    k_vy_old    = k_vy;
%     if timestep == 1
%         k_vy_old(:)=2.5;
%         k_vx_old(:)=2.5;
%     end
    Temp_old    = Temp;
    cp_p_old    = cp_p;
    rhocp_p_old = rhocp_p;
end
if Diffusion_creep == 1
    def_mode_old   =   def_mode;
    straind_old    =   straind;
    D_grain_old    =   D_grain;
    Dd_grain_old    =   Dd_grain;
end


% interpolate rho and viscosity to nodes
Toxx    =   zeros((ny+1)*(nx+1),1);
Toxy    =   zeros((ny+1)*(nx+1),1);
eta_s   =   zeros((ny+1)*(nx+1),1);
eta_p   =   zeros((ny+1)*(nx+1),1);
mu_s    =   zeros((ny+1)*(nx+1),1);
mu_p    =   zeros((ny+1)*(nx+1),1);
rho_vy  =   zeros((ny+1)*(nx+1),1);
rho_vx  =   zeros((ny+1)*(nx+1),1);
rho_p   =   zeros((ny+1)*(nx+1),1);
rho_s   =   zeros((ny+1)*(nx+1),1);
wts     =   zeros((ny+1)*(nx+1),1);
wtp     =   zeros((ny+1)*(nx+1),1);
wtvy    =   zeros((ny+1)*(nx+1),1);
wtvx    =   zeros((ny+1)*(nx+1),1);
Exx     =   zeros((ny+1)*(nx+1),1);
Exy     =   zeros((ny+1)*(nx+1),1);
w       =   zeros((ny+1)*(nx+1),1);
strain  =   zeros((ny+1)*(nx+1),1);
strainv  =   zeros((ny+1)*(nx+1),1);
work    =   zeros((ny+1)*(nx+1),1);
lambda_s=   zeros((ny+1)*(nx+1),1);
if Temperature == 1
    if ra_heat ==1
        Hr      =   zeros((ny+1)*(nx+1),1);
    end
    if adiab_heat == 1 || TPdep_dens == 1
        Ta     =   zeros((ny+1)*(nx+1),1);
    end
    cp_p    =   zeros((ny+1)*(nx+1),1);
    rhocp_p =   zeros((ny+1)*(nx+1),1);
    Temp    =   zeros((ny+1)*(nx+1),1);
    k_vy    =   zeros((ny+1)*(nx+1),1);
    k_vx    =   zeros((ny+1)*(nx+1),1);
    rhocpwt =   zeros((ny+1)*(nx+1),1);
end
if Diffusion_creep == 1
    def_mode   =   zeros((ny+1)*(nx+1),1);
    straind    =   zeros((ny+1)*(nx+1),1);
    D_grain    =   zeros((ny+1)*(nx+1),1);
    Dd_grain    =   zeros((ny+1)*(nx+1),1);
end

