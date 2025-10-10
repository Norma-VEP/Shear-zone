% composing LT RT

if timestep == 1
    matrix_T_eq                     = eta_p;
    
    % Fill up matrix
    kT_ind      =   find(matrix_T_eq ~= 0);
    kT_leftBC   =   1:(ny+1);
    kT_rightBC  =   (ny+1)*(nx)+1:(ny+1)*(nx+1);
    kT_topBC    =   (ny+1)+1:(ny+1):(ny+1)*(nx-1)+1;
    kT_bottomBC =   2*(ny+1):(ny+1):(ny+1)*(nx);
end

num1 = 1; num2 = length(kT_ind);
% Temperature
iiT(num1:num2)   =   kT_ind;
jjT(num1:num2)   =   kT_ind;
ssT(num1:num2)   =  -(k_vy(kT_ind)/dy^2 + k_vy(kT_ind-1)/dy^2 + k_vx(kT_ind)/dx^2 + k_vx(kT_ind-(ny+1))/dx^2) - rhocp_p(kT_ind)/dt;
num1 = num2+1; num2 = num2+length(kT_ind);
iiT(num1:num2)   =   kT_ind;
jjT(num1:num2)   =   kT_ind+1;
ssT(num1:num2)   =   k_vy(kT_ind)/dy^2;
num1 = num2+1; num2 = num2+length(kT_ind);
iiT(num1:num2)   =   kT_ind;
jjT(num1:num2)   =   kT_ind-1;
ssT(num1:num2)   =   k_vy(kT_ind-1)/dy^2;
num1 = num2+1; num2 = num2+length(kT_ind);
iiT(num1:num2)   =   kT_ind;
jjT(num1:num2)   =   kT_ind+(ny+1);
ssT(num1:num2)   =   k_vx(kT_ind)/dx^2;
num1 = num2+1; num2 = num2+length(kT_ind);
iiT(num1:num2)   =   kT_ind;
jjT(num1:num2)   =   kT_ind-(ny+1);
ssT(num1:num2)   =   k_vx(kT_ind-(ny+1))/dx^2;
num1 = num2+1; num2 = num2+length(kT_leftBC);
RT(kT_ind)       =   -Temp(kT_ind).*rhocp_p(kT_ind)./dt;
if shear_heat ~= 0
    RT(kT_ind)       =   RT(kT_ind) - shear_heat*(Txx(kT_ind).^2./eta_p(kT_ind) + Txy(kT_ind-(ny+1)-1).^2./(4.*eta_s(kT_ind-(ny+1)-1)) + Txy(kT_ind-(ny+1)).^2./(4.*eta_s(kT_ind-(ny+1))) + Txy(kT_ind).^2./(4.*eta_s(kT_ind)) + Txy(kT_ind-1).^2./(4.*eta_s(kT_ind-1)));
end
if ra_heat == 1
    RT(kT_ind)       =   RT(kT_ind) - Hr(kT_ind);
end
if adiab_heat == 1
    RT(kT_ind)       =   RT(kT_ind) - (Temp(kT_ind).*Ta(kT_ind).*rho_p(kT_ind).*(gravity_x * (Vx(kT_ind-(ny+1)) + Vx(kT_ind))./2 + gravity_y * (Vy(kT_ind-1) + Vy(kT_ind))./2));
end

% Left BC
iiT(num1:num2)   =   kT_leftBC;
jjT(num1:num2)   =   kT_leftBC;
ssT(num1:num2)   =   1;
num1 = num2+1; num2 = num2+length(kT_leftBC);
iiT(num1:num2)   =   kT_leftBC;
jjT(num1:num2)   =   kT_leftBC+(ny+1);
ssT(num1:num2)   =  -1;
num1 = num2+1; num2 = num2+length(kT_rightBC);
RT(kT_leftBC)    =   0;
% Right BC
iiT(num1:num2)   =   kT_rightBC;
jjT(num1:num2)   =   kT_rightBC;
ssT(num1:num2)   =   1;
num1 = num2+1; num2 = num2+length(kT_rightBC);
iiT(num1:num2)   =   kT_rightBC;
jjT(num1:num2)   =   kT_rightBC-(ny+1);
ssT(num1:num2)   =  -1;
num1 = num2+1; num2 = num2+length(kT_topBC);
RT(kT_rightBC)   =   0;
% Top BC
iiT(num1:num2)   =   kT_topBC;
jjT(num1:num2)   =   kT_topBC;
ssT(num1:num2)   =   0.5;
num1 = num2+1; num2 = num2+length(kT_topBC);
iiT(num1:num2)   =   kT_topBC;
jjT(num1:num2)   =   kT_topBC+1;
ssT(num1:num2)   =   0.5;
num1 = num2+1; num2 = num2+length(kT_bottomBC);
RT(kT_topBC)     =    top_T;
% Bottom BC
if Ext_T == 1 || Ext_T == 2
    iiT(num1:num2)   =   kT_bottomBC;
    jjT(num1:num2)   =   kT_bottomBC;
    ssT(num1:num2)   =   0.5;
    num1 = num2+1; num2 = num2+length(kT_topBC);
    iiT(num1:num2)   =   kT_bottomBC;
    jjT(num1:num2)   =   kT_bottomBC-1;
    ssT(num1:num2)   =   0.5;
else
    iiT(num1:num2)   =   kT_bottomBC;
    jjT(num1:num2)   =   kT_bottomBC;
    ssT(num1:num2)   =   1;
    num1 = num2+1; num2 = num2+length(kT_topBC);
    iiT(num1:num2)   =   kT_bottomBC;
    jjT(num1:num2)   =   kT_bottomBC-1;
    ssT(num1:num2)   =  -Ext_T;
end
RT(kT_bottomBC)  =   bottom_T;

LT = sparse(iiT,jjT,ssT,(nx+1)*(ny+1),(nx+1)*(ny+1));

T_new = LT\RT;

% Computing DT
DT  =   T_new - Temp;

%%%   SUBGRID DIFFUSION
DTsubgrid   =   zeros((ny+1)*(nx+1),1);
wtsubgrid   =   zeros((ny+1)*(nx+1),1);

j = fix((xm+dx/2)/dx)+1;
i = fix((ym+dy/2)/dy)+1;

ind = j<1;
j(ind) = 1;
ind = j>nx-1;
j(ind) = nx-1;
ind = i<1;
i(ind) = 1;
ind = i>ny-1;
i(ind) = ny-1;

dxm =   xm - x_P(j);
dym =   ym - y_P(i);

ind_T = (j-1).*(ny+1)+i;

T_old     =      Temp(ind_T)'           .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
    +Temp(ind_T+(ny+1))'    .* (dxm/dx)       .* (1 - (dym/dy)) ...
    +Temp(ind_T+1)'         .* (1 - (dxm)/dx) .* (dym/dy) ...
    +Temp(ind_T+(ny+1)+1)'  .* (dxm/dx)       .* (dym/dy);

deltaT      =   Tm - T_old;
deltaT_new  =   deltaT .* exp((-2*km*dt)./(rhocpm * dx^2));
Tdelta      =   deltaT_new - deltaT;
Tm          =   Tm + Tdelta;

wt                          =   (1-abs(xm-x_P(j+1))/dx).*(1-abs(ym-y_P(i+1))/dy);
DTsubgrid(ind_T+(ny+1)+1)   =   DTsubgrid(ind_T+(ny+1)+1)' + Tdelta .* rhocpm .* wt;
wtsubgrid(ind_T+(ny+1)+1)   =   wtsubgrid(ind_T+(ny+1)+1)' + wt .* rhocpm;

ind = wtsubgrid > 0;
DTsubgrid(ind)  = DTsubgrid(ind)./wtsubgrid(ind);

% COMPUTE remaining temp changes on nodes
DT = DT-DTsubgrid;

% interpolating DT to markers
j = fix((xm+dx/2)/dx)+1;
i = fix((ym+dy/2)/dy)+1;

ind = j<1;
j(ind) = 1;
ind = j>nx;
j(ind) = nx;
ind = i<1;
i(ind) = 1;
ind = i>ny;
i(ind) = ny;

dxm =   xm - x_P(j);
dym =   ym - y_P(i);

ind_T = (j-1).*(ny+1)+i;

dTnodal     =    DT(ind_T)'         .* (1 - dxm/dx) .* (1 - dym/dy) ...
    +DT(ind_T+(ny+1))'  .* (dxm/dx)     .* (1 - dym/dy) ...
    +DT(ind_T+1)'       .* (1 - dxm/dx) .* (dym/dy) ...
    +DT(ind_T+(ny+1)+1)'.* (dxm/dx)     .* (dym/dy);

Tm  =   Tm + dTnodal;


% MATCH MARKER AND NODAL T FOR THE FIRST TIMESTEP

if timestep==1
    % INTERPOLATE NEW TEMPERATURE FROM NODES TO MARKERS
    j = fix((xm+dx/2)/dx)+1;
    i = fix((ym+dy/2)/dy)+1;
    
    ind = j<1;
    j(ind) = 1;
    ind = j>nx;
    j(ind) = nx;
    ind = i<1;
    i(ind) = 1;
    ind = i>ny;
    i(ind) = ny;
    
    dxm =   xm - x_P(j);
    dym =   ym - y_P(i);
    
    Tm  =    T_new(ind_T)'          .* (1 - dxm/dx) .* (1 - dym/dy) ...
        +T_new(ind_T+(ny+1))'   .* (dxm/dx)     .* (1 - dym/dy) ...
        +T_new(ind_T+1)'        .* (1 - dxm/dx) .* (dym/dy) ...
        +T_new(ind_T+(ny+1)+1)' .* (dxm/dx)     .* (dym/dy);
    
    %==========================================================================
    %   figure(5);clf;pcolor(T2d); shading interp; colorbar;axis ij image;pause(3)
    %         a=1000
    %         pause(3);
end
