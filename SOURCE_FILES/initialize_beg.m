% length of delta x and delta y
dx  =   Lx/(nx-1);
x   =   0:dx:Lx;
dy  =   Ly/(ny-1);
y   =   0:dy:Ly;

% initiate coordinates
x_Vy    =   -dx/2:dx:Lx+dx/2;
y_Vy    =   0:dy:Ly+dy;

x_Vx    =   0:dx:Lx+dx;
y_Vx    =   -dy/2:dy:Ly+dy/2;

x_P     =   -dx/2:dx:Lx+dx/2;
y_P     =   -dy/2:dy:Ly+dy/2;

% % initiate 2D grids
% % Basic points
% [x2d_b,y2d_b]   =   meshgrid(x_Vx(1:end-1),y_Vy(1:end-1));
% % Pressure points
% [x2d_p,y2d_p]   =   meshgrid(x_P,y_P);
% % Velocity points (without ghost points)
% [x2d_vx,y2d_vx]   =   meshgrid(x_Vx(1:end-1),y_Vx);
% [x2d_vy,y2d_vy]   =   meshgrid(x_Vy,y_Vy(1:end-1));

% Number of Markers
nxm     =   mx*(nx-1);
nym     =   my*(ny-1);
marknum =   nxm*nym;

left_income   = zeros(1,ny-1);
left_income_y = zeros(1,ny);
right_income  = zeros(1,ny-1);
right_income_y= zeros(1,ny);
top_income    = zeros(1,nx-1);
bottom_income = zeros(1,nx-1);

xm      =   zeros(1,marknum);   % x-coordinates of markers
ym      =   zeros(1,marknum);   % x-coordinates of markers
etam    =   zeros(1,marknum);   % viscosity of markers
mum     =   zeros(1,marknum);   % elasticity of markers
wm      =   zeros(1,marknum);   % rotation of markers
rhom0   =   zeros(1,marknum);   % initial density of markers
phim    =   zeros(1,marknum);   % friction angles of markers
Cm      =   zeros(1,marknum);   % cohesion of markers
Txxm    =   zeros(1,marknum);   % pure stress of markers
Txym    =   zeros(1,marknum);   % shear stress of markers
txxm    =   zeros(1,marknum);   % to save
txym    =   zeros(1,marknum);   % to save
dtxxm    =   zeros(1,marknum);   % pure stress change for timne step of markers
dtxym    =   zeros(1,marknum);  % shear stress change for timne step of markers
Exxm    =   zeros(1,marknum);   % plane strain rate on markers
Exym    =   zeros(1,marknum);   % Shear strain rates on markers
Exxn    =   zeros(1,marknum);   % plane strain rate on markers (averaged from nodes)
Exyn    =   zeros(1,marknum);   % Shear strain rates on markers (averaged from nodes)
E2ndm   =   zeros(1,marknum);   % 2nd invariant of strain rate of markers
T2ndm   =   zeros(1,marknum);   % 2nd invariant of stress of markers
E2ndn   =   zeros(1,marknum);   % 2nd invariant of strain rate of markers (averaged from nodes)
Pm      =   zeros(1,marknum);   % Pressure of markers
Im      =   zeros(1,marknum);   % Phase of markers
lambdam =   zeros(1,marknum);   % lambda of markers
strainm =   zeros(1,marknum);   % plastic strain of markers
strainvm =   zeros(1,marknum);   % plastic strain of markers
workm   =   zeros(1,marknum);   % work of markers
alpha   =   zeros(1,marknum);   % rotation of markers
if Temperature == 1
    if ra_heat == 1
        hrm     =   zeros(1,marknum); % radioactive heting of markers
    end
    if adiab_heat == 1 && TPdep_dens == 0
        tam     =   zeros(1,marknum); % adiabatic heat of markers
    end
    if TPdep_dens == 1
        tam     =   zeros(1,marknum);   % expansion factor of markers
        tbm     =   zeros(1,marknum);   % expansion factor of markers
    end 
    Tm      =   zeros(1,marknum);   % Temperature of markers
    km      =   zeros(1,marknum);
    cpm     =   zeros(1,marknum);
    rhocpm0 =   zeros(1,marknum);
end
if Powerlaw == 1
    RK = 8.3145;
    nm      =   zeros(1,marknum);   % power-law component of markers
    Qm      =   zeros(1,marknum);   % activation energy of markers
end
if Diffusion_creep == 1
    dm          =   zeros(1,marknum);   % grain size
    ddm         =   zeros(1,marknum);   % growth rate
    straindm    =   zeros(1,marknum);   % diffusion creep strain
    mm          =   zeros(1,marknum);   % diffusion creep coefficient
    etadm       =   zeros(1,marknum);   % diffusion creep pre
    eta_diffm   =   zeros(1,marknum);   % diffusion creep pre
    def_modem   =   zeros(1,marknum);   % diff vs. disl creep
    Qdm         =   zeros(1,marknum);
end
    
eta_s   =   zeros((ny+1)*(nx+1),1);     % viscosity of nodes corners
Toxy    =   zeros((ny+1)*(nx+1),1);     % old shear stress on nodes
mu_s    =   zeros((ny+1)*(nx+1),1);     % elastic module on node corners
w       =   zeros((ny+1)*(nx+1),1);
rho_s   =   zeros((ny+1)*(nx+1),1);     % density on corner of nodes
eta_p   =   zeros((ny+1)*(nx+1),1);     % viscosity in the center of nodes
Toxx    =   zeros((ny+1)*(nx+1),1);     % old pure stress
mu_p    =   zeros((ny+1)*(nx+1),1);     % elastic module in node center
rho_p   =   zeros((ny+1)*(nx+1),1);     % density on node center
rho_vy  =   zeros((ny+1)*(nx+1),1);     % density on vy location
rho_vx  =   zeros((ny+1)*(nx+1),1);     % density on vx location
strain  =   zeros((ny+1)*(nx+1),1);     % strain on nodes
strainv  =   zeros((ny+1)*(nx+1),1);     % strain on nodes
work    =   zeros((ny+1)*(nx+1),1);     % work on nodes
lambda_s=   zeros((ny+1)*(nx+1),1);     % lambda on node corners
if Temperature == 1
    if ra_heat == 1
        Hr      =   zeros((ny+1)*(nx+1),1);
    end
    if adiab_heat == 1 || TPdep_dens == 1
        Ta      =   zeros((ny+1)*(nx+1),1);
    end 
    k_vx    =   zeros((ny+1)*(nx+1),1);
    k_vy    =   zeros((ny+1)*(nx+1),1);
    Temp    =   zeros((ny+1)*(nx+1),1);
    cp_p    =   zeros((ny+1)*(nx+1),1);
    rhocp_p =   zeros((ny+1)*(nx+1),1);
end
if Diffusion_creep == 1
    def_mode   =   zeros((ny+1)*(nx+1),1);
    straind    =   zeros((ny+1)*(nx+1),1);
    D_grain    =   zeros((ny+1)*(nx+1),1);
    Dd_grain   =   zeros((ny+1)*(nx+1),1);
end
    
% Matrix formation
R   =   zeros((nx+1)*(ny+1)*3,1);
LT  =   spalloc((nx+1)*(ny+1),(nx+1)*(ny+1),((nx+1)*(ny+1))*11);    % number of temperature coefficients
RT  =   zeros((nx+1)*(ny+1),1);
LS  =   spalloc((nx-1)*surface_nodes+1,(nx-1)*surface_nodes+1,((nx-1)*surface_nodes+1)*5);
RS  =   zeros((nx-1)*surface_nodes+1,1);
% Right hand side
S   =   zeros((nx+1)*(ny+1)*3,1);

iterations=zeros(2,1);

