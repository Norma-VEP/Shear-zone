%% 3D Conservative Stokes / Marker in cell / plastic iterations
% VisoElastoPlasic approach
% by Jonas, March/April 2012
clear all; 

addpath('SOURCE_FILES')                        % FUNCTIONS
mkdir OUTPUT

create_output       = 50;
create_breakpoint   = 200;

brutus              = 1;    % 1 if you run it on cluster (without plotting every timestep)
Temperature         = 2;
Powerlaw            = 1;
Diffusion_creep     = 1;
Surface_diffusion   = 0;
Beam_function       = 0;    % 1: beam function of elastic bending at the bottom. Requires material below the bottom and no boundary velocities
delta_prograd       = 0;    % 1: delta progradation on the left side of the model

if brutus == 0
    cd COLORMAPS
        mat = dir('*.mat');
        for i=1:length(mat)
            load(mat(i).name);
        end    
    cd ..
end

% Length of model
Lx   =  0.05;
Ly   =  0.05;

% number of nodes
nx      =   401;
ny      =   401;
p_node  =   201;   % between 3 and nx-2

% markers per node
mx  =   4;
my  =   4;

        % Phase    Visc       EModul    Dens     Phi     PhiW      Coh         CohW      lamb      Hr       Ta      Tb       cp      Temp       kk       n       Q      DifCreep    grains      m      Qd
ROCKS= [    1      1e18        1e11        1       0       1       1e20        1e20        0        0        0       0      3e6       273      200       1       0             0        0       0       0   ;   % AIR
%             2      1e18        1e11      1e3       0       1       1620        1e20        0        0        0       0      3e3       273      200       1       0             0        0       0       0   ;   % WATER
            2   1.9953e21      1e11     2800      35      35       25e6        25e6      0.0     2e-6    2e-5   45e-13      1e3       273      2.5      3.0  220e3     1.5849e18     2e-2      -2   220e3   ;   % Quartz (Brodie and Rutter, 2000)
            3   5.7143e34      1e11     2800      35      35       25e6        25e6      0.00    2e-6    2e-5   45e-13      1e3       273      2.5      4.0  125e3     1.5849e18     2e-2      -2   220e3   ;   % Quartz (Tokle et al. 2019; Brodie and Rutter, 2000)
            4   3.1623e19      1e11     2700      35      35       25e6        25e6      0.0     2e-7    2e-5   45e-13      1e3       273      1.6      3.0  235e3     7.9433e22     2e-2      -3   153e3   ;   % Wet Anorthite60 (Rybacki and Dresen, 2004)
            5   2.5119e15      1e11     2700      35      35       25e6        25e6      0.0     2e-7    2e-5   45e-13      1e3       273      1.6      3.0  345e3     1.9953e22     2e-2      -3   159e3   ;   % Wet Anorthite (Rybacki et al., 2006; siehe visc calc)
            6    4.88e133      1e11     3000       5       5      130e6       130e6      0.0     2e-6    2e-5   45e-13      1e3       273      2.5     18.0   51e3             0        0       0       0  ];   % Biotite (Shea and Kronenberg, 1992)
%             7    7.7e252       1e11     3000       5       5      130e6       130e6      0.0     2e-6    2e-5   45e-13      1e3       273      2.5     31.0   98e3             0        0       0       0  ];   % Biotite (Shea and Kronenberg, 1992)
            
                 
Phase   =   ROCKS(:,1)';
eta     =   ROCKS(:,2)';
mu      =   ROCKS(:,3)';
rho     =   ROCKS(:,4)';
phi     =   ROCKS(:,5)';
phi_w   =   ROCKS(:,6)';
C       =   ROCKS(:,7)';
C_w     =   ROCKS(:,8)';
lambda  =   ROCKS(:,9)';
if Temperature == 1
    hr  =   ROCKS(:,10)';
    ta  =   ROCKS(:,11)';
    tb  =   ROCKS(:,12)';
    cp  =   ROCKS(:,13)';
    T   =   ROCKS(:,14)';
    kk  =   ROCKS(:,15)';
end
if Powerlaw == 1
    n   =   ROCKS(:,16);
    Q   =   ROCKS(:,17);
end
if Diffusion_creep == 1
    dCr =   ROCKS(:,18);
    Gr  =   ROCKS(:,19);    
    m   =   ROCKS(:,20);
    Qd  =   ROCKS(:,21);
end
    
variable_lambda=0;      % 1: initial lambda; 2: equalized whole column; 3: equalized top 4 km
lambda_bottom=0.95;
lambda_increase=0.00;
OH_const = 50e6;
Lith = 100e3;

Tm              =   550+273;        % Temperature: +273 for conversion to Kelvin
geoT            =   25;             % Geotherm: ?C per kilometer depth
p_init          =   (Tm-273)/geoT*1000*9.81*2700; % Pressure at depth of 600?C (9.81 = gravity constant | 2700 = density of overlying rock mass)


% Geometry      Phase       x-min       x-max       y-min       y-max
SETUP   =   [   5           0           Lx          0           Ly      ];
      
% Plastic weakeing thresholds
w_1 = 0.1;
w_2 = 1.0;

% Parameters
gravity_y       =   0;
gravity_x       =   0.0;
SecYear         =   3600*24*365.25;
% p_init          =   1e5;

% Surface process
SedimentationStyle  =   1;          % 1 = sed and ero || 2 = only sed || 3 = only ero || 4 = linear below sealevel
SurfaceCoeff        =  [1e-6];   % Coefficient for [SEDIMENTATION EROSION] FOR DIFFERENT SURFACE PROCESSES, OTHERWISE ONLY ONE COEFFICIENT
SurfaceBC_left      =   -1;       % For free slip: -1
SurfaceBC_right     =   -1;       % For free slip: -1
Surface_dt          =   1;         % makes surface process every xx timestep
time_old            =   0;
SediMarker          =  [ ];       % Marker type for sedimentation
ErosMarker          =   1;          % Marker type for erosion
SedChange           =   2e6;        % Marker change for sedimentation in years
WaterLevel          =   0;          % If WaterLevel == 0 it is switched off !!
SediRate            =   0;      % Linear sedimentation rate in m/yr
surface_nodes       =   5;          % Times x nodes along surface
surface_init        =   10000;
surface_x           =   [0:(Lx/(nx-1))/surface_nodes:Lx];
surface_y           =   [surface_init*ones(1,length(surface_x))];
surface_smoother    =   3;


% Cutoff viscosities
eta_max = 1e24;
eta_min = 1e14;
% eta_bingham = 5e18;

if Temperature==1
    % Temperature boundary conditions
    A_thick     =   10e3; % Thickness of sticky-air
    C_thick     =   33e3; % Thickness of Crust
    L_thick     =   Lith; % Thickness of Lithosphere
    Moho_temp   =   660;  % Temperature at Moho in ?C
    L_A_temp    =   1270+Lith/2000; % Temperature at Lithosphere/Asthenosphere boundary in ?C
    M_grad      =   0.5;  % Temperature gradient in Mantle (background)
    Pert_beg    =   490e3;
    Pert_end    =   510e3;
    Pert_add    =   0e3; % km difference for Lithosphere/Asthenosphere depth
    top_T       =   273;    % temperature in case of shear model (in Kelvin)
    Ext_T       =   0.996;
    if Ext_T == 1
        bottom_T    =   L_A_temp+273+M_grad*(Ly-A_thick-L_thick)/1e3;
    elseif Ext_T == 2
        bottom_T    =   top_T;
    else
        bottom_T    =   ((L_A_temp+273)+M_grad*(Ly-L_thick-A_thick)/1000+(Ly/(ny-1))/1000*M_grad/2) - Ext_T*((L_A_temp+273)+M_grad*(Ly-L_thick-A_thick)/1000-(Ly/(ny-1))/1000*M_grad/2);
    end
    shear_heat  =   0.99;
    ra_heat     =   1;
    adiab_heat  =   1;
    TPdep_dens  =   0;
end

% Velocity boundary conditions

%   L T T T T T T T T T T T RR      T T T T T T T T T T T T T   
%   L                       RR      L                       R
%   L      x-velocity       RR      L      y-velocity       R
%   L                       RR      L                       R
%   L                       RR      B B B B B B B B B B B B B
%   L B B B B B B B B B B B RR      B B B B B B B B B B B B B


BC_left_init     =  0;  % no need if mirrored
BC_right_init    =  0;  % no need if mirrored
BC_top_init      =  0.00025*1e-3;
BC_bottom_init   = -0.00025*1e-3;
bound            = {'freeslip','noslip','velocity','external','mixed','mirror'};
top_BC           = bound(3);
bottom_BC        = bound(3);
left_BC          = bound(6);
right_BC         = bound(6);

% Inversion
inversion = [100e6 100e6 150e6 150e6];     % Time in Ma for velocity inversion

% Incoming marker type
mark_top    = 1;
mark_bottom = 17;

%=============================================
% Initialize matrices and coordinates
initialize_beg;
%=============================================

% Initial marker pattern
pattern_marker  =   [   ];     % Types of markers with pattern (check marker phases !!!!)
pattern_type    =   {'horizontal','vertical','kaki'};   % Stripes or kaki
pattern_type    =   pattern_type(1);
pattern_xdim    =   2e3;    % thickness of vertical stripes
pattern_ydim    =   2e3;    % thickness of horizontal stripes

%=============================================
% Initialize phase and temperature distribution of marker 
marker_distribution;
%=============================================

% Initial time and iteration setup
time        = 0;
t_beg       = 1;
dt_value    = 0.2;          % threshold value for maximal timestep, 0.1 = 10% movement of dx or dy
Ddt         = 100*SecYear;  % Initial timestep to pre-stress
short_dt    = 100*SecYear;  % Short initial timestep
dt_max      = 100*SecYear;  % Timestep
n_short_dt  = 1;            % Number of short initial timesteps
miniter     = 1;             % minimal number of iterations
maxiter     = 1;             % maximal number of iterations
error       = 1;             % 1 = Error 1 (average nodal velocity change); 2 = Error 2 (largest nodal velocity change)
vel_res     = 1e-14;         % sum(velocity) change fot iter break
strain_rate_smoother = 0.87;

% Load breakpoint file if necessary
if exist(char('Breakpoint.mat'))
    load Breakpoint.mat
    t_beg    = timestep+1;
end

%=========================== START TIME LOOP ==========================
%======================================================================

for timestep = t_beg:2e6
    tic

%     if timestep>100
%         left_BC          = bound(1);
%         right_BC         = bound(1);
%     end
    
    dt=dt_max;
    if timestep <= n_short_dt
        dt = short_dt;  
    end
    
    velocity_inversion;
    
    for niter = 1:maxiter

        fprintf('Timestep: %d\n', timestep);
        fprintf('Iteration: %d\n', niter);
        
        %=============================================
        % Reload old values and initialize matrices
        initialize_iter;
        %=============================================
        
        %=============================================
        % Calculate power-law and brittle viscosities
        viscosity_calculation_wet;
        %=============================================
        
        %=============================================
        % Fill nodal values from marker information
        marker_to_nodes; % dispatched loops, better for desktop
        toc, fprintf('for updating nodal values!\n'); 
        %=============================================   
        
        if Temperature == 1            
            % Boundary conditions to interpolated T
            % Upper BC
            Temp((ny+1)+1:(ny+1):(nx-1)*(ny+1)+1)   =   top_T*2 - Temp((ny+1)+2:(ny+1):(nx-1)*(ny+1)+2);
            % Lower BC
            if Ext_T == 1 || Ext_T == 2
                Temp(2*(ny+1):(ny+1):(nx)*(ny+1))       =   bottom_T'*2 - Temp(2*(ny)+1:(ny+1):(nx)*(ny+1)-1);
            else
                Temp(2*(ny+1):(ny+1):(nx)*(ny+1))       =   bottom_T + Ext_T*Temp(2*(ny)+1:(ny+1):(nx)*(ny+1)-1);
            end
            % Left BC
            Temp(1:ny+1)                            =   Temp(ny+2:2*(ny+1));
            % Right BC
            Temp(nx*(ny+1)+1:(nx+1)*(ny+1))         =   Temp((nx-1)*(ny+1)+1:(nx)*(ny+1));
        end
        
        %=====================================================
        
        %=============================================
        % Fill and solve matrix for stokes
        stokes_direct_solver;
        toc, fprintf('for solving the matrix!\n')
        %=============================================
        
        Vx    =   S(1:(ny+1)*(nx+1));
        Vy    =   S(1+(ny+1)*(nx+1):2*(ny+1)*(nx+1));
        P     =   S(2*(ny+1)*(nx+1)+1:end).*kcont;
        
        Vx2d    =   reshape(Vx,ny+1,nx+1);
        Vy2d    =   reshape(Vy,ny+1,nx+1);
        P2d     =   reshape(P,ny+1,nx+1);
                
        %======================================================================
        % redistributing S
        
        %   1--6--11--16
        %   |  |  |    |
        %   2--7--12--17
        %   |  |  |    |
        %   3--8--13--18
        %   |  |  |    |
        %   4--9--14--19
        %   |  |  |    |
        %   5-10--15--20
        
        Vx2d_s   =   (Vx2d(1:end-1,1:end-1)+Vx2d(2:end,1:end-1))./2;
        Vy2d_s   =   (Vy2d(1:end-1,1:end-1)+Vy2d(1:end-1,2:end))./2;
        P2d_p    =   P2d;
        %     T2dn     =   T2d;
        
        % define optimal timestep
        Ddt = dt;
        
        Vx_max  =   max(max(abs(Vx2d_s)));
        Vy_max  =   max(max(abs(Vy2d_s)));
        
        if dt_value*1/(Vx_max/dx + Vy_max/dy) < Ddt
            Ddt  =   dt_value*1/(Vx_max/dx + Vy_max/dy);
        end
                
        %=============================================
        % Calculating strain rates and stresses
        strain_rate_stress_calculation;
        toc, fprintf('for strain/stress calculation and interpolation to markers!\n');
             fprintf('========= Ddt = %d years =============================\n',Ddt/SecYear);
        %=============================================
        
        fprintf('========= Min Pressure = %d MPa ===========================\n',min(P)/1e6);
        fprintf('========= Max Pressure = %d MPa ===========================\n',max(P)/1e6);
   
        %=============================================
        % Calculating strain rates and stresses
        vel_res_calculation;
        if timestep>1
            fprintf('========= VELOCITY ERROR %d = %d =========\n\n',error,iterations(error,niter+maxiter*(timestep-1)));
            fprintf('========= Pressure change on node = %d =========\n\n',dP_node(timestep));
        end
        %=============================================

        % Plotting if running on desktop
        if brutus==0 && timestep>1
            E2nd_s_2d   =   reshape(E2nd_s,ny+1,nx+1);
            eta_s_2d    =   reshape(eta_s,ny+1,nx+1);
            strain_2d   =   reshape(strain,ny+1,nx+1);
            strainv_2d   =   reshape(strainv,ny+1,nx+1);
            lambda_s_2d   =   reshape(lambda_s,ny+1,nx+1);
            
            figure(1), clf
            colormap jet
            
            subplot(231)
            pcolor(x_Vx(1:end-1),y_Vy(1:end-1),log10(E2nd_s_2d(1:end-1,1:end-1)))%, caxis([-14 -10])
            shading interp
            colorbar
            axis image, axis ij
            title('E2nd [1/s]')
            
            subplot(232)
            pcolor(x_Vx(1:end-1),y_Vy(1:end-1),log10(eta_s_2d(1:end-1,1:end-1))), caxis(log10([eta_min eta_max]))
            shading interp
            colorbar
            axis image, axis ij
            title(['\eta_{s} [Pa.s] ',    num2str(Ddt/SecYear)])
            drawnow

            subplot(233)
            pcolor(x_Vx(1:end-1),y_Vy(1:end-1),(strain_2d(1:end-1,1:end-1)))
            shading interp
            colorbar
            axis image, axis ij
            title(['Strain [-] ',    num2str(Ddt/SecYear)])
            drawnow
            
            subplot(234)
            pcolor(x_Vx(1:end-1),y_Vy(1:end-1),(T2nd_s_2d(1:end-1,1:end-1)))
            shading interp
            colorbar
            axis image, axis ij
            title(['\tau_{s} [Pa] ',    num2str(Ddt/SecYear)])
            drawnow

            subplot(235)
            pcolor(x_Vx(1:end-1),y_Vy(1:end-1),Vx2d_s*SecYear)
            shading interp
            colorbar
            axis image, axis ij
            title(['V_{x} [m/yr] ',    num2str(Ddt/SecYear)])
            drawnow

            subplot(236)
            pcolor(x_Vx(1:end-1),y_Vy(1:end-1),Vy2d_s*SecYear)
            shading interp
            colorbar
            axis image, axis ij
            title(['V_{y} [m/yr] ',    num2str(Ddt/SecYear)])
            drawnow

            
            figure(444)
            subplot(211), semilogy((niter-1)/maxiter+(timestep-1),iterations(1,niter+maxiter*(timestep-1)),'ko'), hold on
            title('Error 1: with sticky-air')
            subplot(212), semilogy((niter-1)/maxiter+(timestep-1),iterations(2,niter+maxiter*(timestep-1)),'ko'), hold on
            title('Error 2: with sticky-air')
        end
        
        % Exiting the iteration loop
        if (niter>=miniter && iterations(error,niter+maxiter*(timestep-1))<vel_res) || timestep==1
            break;
        end
    end

    %==========================================================================
    dt = Ddt;
    
    % Calculating accumulated plastic/viscous strain and grain size
    Strain_GrainSize_GSE;
    %==========================================================================

    %=============================================
    % Calculating temperature
    subgrid_diffusion_stresses;
    stress_rotation;
    % Average deviatoric stress
    Sigma_d(1,timestep)     =   sum((Txxm.^2 + Txym.^2).^0.5)./length(T2ndm);
    Visc_d(1,timestep)      =   sum(eta_effm)./length(eta_effm);
    ind=find((Im==2 | Im==3) & ind_plast==0);
    D_qtz_d(1,timestep)     =   sum(dm(ind))/length(ind);
    def_mode_qtz_d(1,timestep)     =   sum(def_modem(ind))/length(ind);
    ind=find((Im==4 | Im==5) & ind_plast==0);
    D_plg_d(1,timestep)     =   sum(dm(ind))/length(ind);
    def_mode_plg_d(1,timestep)     =   sum(def_modem(ind))/length(ind);
    %=============================================
    
    %=============================================
    % Calculating temperature
    if Temperature == 1
        temperature_direct_solver;
        toc, fprintf('for temperature calculation!\n');
    end
    %=============================================


    %=============================================
    % Move markers
    move_marker;
%     move_surface;
    if Beam_function == 1
        beam_function;
    end

    toc, fprintf('for moving makers!\n');
    %=============================================

    time = time + dt;

    %======================================================================
    % Visualization begin
    %======================================================================


    if brutus == 0  figure(3), clf
        
        for i=1:20
            Phase1  =   find(Im == i);
            marksize = 5;
            hold on
            plot(xm(Phase1),ym(Phase1),'.','Color',Colormaps(i,:),'MarkerSize',marksize)
        end
      
        
%         plot(surface_x,surface_y,'r'), hold on
        
        axis image, axis ij
        title(['Time: ',num2str((time)/SecYear/1e6),' Ma'])
%         plot(x2d_b,y2d_b,'k',x2d_b',y2d_b','k')
%         quiver(x2d_b,y2d_b,Vx2db,Vy2db,'r','LineWidth',1.5)
        E2nd_s_2d   =   reshape(E2nd_s,ny+1,nx+1);
        eta_s_2d    =   reshape(eta_s,ny+1,nx+1);
        T2nd_s_2d   =   reshape(T2nd_s,ny+1,nx+1);
                
        hold off
                
%         figure(8)
%         plot(surface_x,surface_y,'b'), hold off
%         axis ij
        
        figure(4), clf
        colormap jet
        
        subplot(311)
        pcolor(x_Vx(1:end-1),y_Vy(1:end-1),log10(E2nd_s_2d(1:end-1,1:end-1)))%, caxis([-14 -10])
        shading interp
        colorbar
        axis image, axis ij
        title('E2nd ')
        
        subplot(312)
        pcolor(x_Vx(1:end-1),y_Vy(1:end-1),log10(eta_s_2d(1:end-1,1:end-1))), caxis(log10([eta_min eta_max]))
        shading interp
        colorbar
        axis image, axis ij
        title(['\eta_{s}   ',    num2str(dt/SecYear)])
        
        subplot(313)
        pcolor(x_Vx(1:end-1),y_Vy(1:end-1),(T2nd_s_2d(1:end-1,1:end-1)))
        shading interp
        colorbar
        axis image, axis ij
        title(['T2nd'])
        drawnow
    end

    %==========================================================================
    % Visualization end
    %==========================================================================

    %=============================================
    % Delete markers out of grid
%     outgrid_marker;
    %=============================================

    %=============================================
    % New markers from sides
%     incoming_marker;
%     toc, fprintf('for outgoing/incoming markers!\n');
    %=============================================
    
    %=============================================
    % Surface process
    if  Surface_diffusion == 1
%         if WaterLevel > 0            
%             if SedimentationStyle == 4
%                 ind=find(Im==2 & ym>WaterLevel);
%                 if round(time/SecYear/SedChange)<time/SecYear/SedChange && length(SediMarker)==2
%                     ii=1;
%                 else
%                     ii=2;
%                 end
%                 Im(ind) = SediMarker(ii);
%                 
%                 ind=surface_y>WaterLevel;
%                 surface_y(ind)=WaterLevel;
%             end
%             if SedimentationStyle == 5
%                 pf_x1=Lx/2-time/SecYear*0.0025;
%                 pf_x2=Lx/2+time/SecYear*0.0025;
%                 
%                 ind=find(-(surface_x-pf_x1)*tand(30)+WaterLevel<surface_y & surface_x<=Lx/2);
%                 surface_y(ind)=-(surface_x(ind)-pf_x1)*tand(30)+WaterLevel;
%                 ind=find((surface_x-pf_x2)*tand(30)+WaterLevel<surface_y & surface_x>Lx/2);
%                 surface_y(ind)=(surface_x(ind)-pf_x2)*tand(30)+WaterLevel;
%                 ind=surface_y<WaterLevel;
%                 surface_y(ind)=WaterLevel;
%             end        
%         end
        
        surface_calculation;
        toc, fprintf('for surface process!\n');
    end
    %=============================================
    work_ratem  = (Txxm.^2)./eta_effm + (Txym.^2)./eta_effm;
    workm       = workm + work_ratem.*dt;


    fprintf('dt   = %d years\n', dt/SecYear)
    fprintf('Time = %d years\n\n', time/SecYear)

    toc, fprintf('for complete timestep\n===============\nMarknum: %d\n===============\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n', marknum);    

    %======================== SAVE OUTPUT FILE ============================
    if mod(timestep,create_output)==0 || timestep==2 || timestep==3
        
        fname = ['Shearing_',num2str(1e6+timestep),'.mat'];
        cd OUTPUT
        save([char(fname)]   ,'Vx2d_s','Vy2d_s','E2nd_s','T2nd_s','xm','Txx','Txy',...
            'ym','Im','eta_s','P2d','Lx','Ly','strain','strainv','straind','lambda_s','rho_s','dm','Sigma_d','Visc_d',...
            'nx','ny','time','dt','SecYear','x','y','iterations','Pm','def_mode','D_grain','Dd_grain','work','D_qtz_d','D_plg_d','def_mode_qtz_d','def_mode_plg_d')
        cd ..
    end

    % Save breakpoint file
    if mod(timestep,create_breakpoint)==0
        fname = ['Breakpoint1.mat'];
        save(char(fname))
%         break
        movefile('Breakpoint1.mat','Breakpoint.mat');
    end
    if time/SecYear >= 1e6 && time/SecYear < 1e6+1e2
        fname = ['Breakpoint_10g.mat'];
        save(char(fname))
%         break;
    end
    if time/SecYear >= 2e6 && time/SecYear < 2e6+1e2
        fname = ['Breakpoint_20g.mat'];
        save(char(fname))
        break;
    end

end
