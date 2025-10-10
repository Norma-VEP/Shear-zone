% Strain rates
ind = [1:(ny+1)*nx]';
Exx(ind+(ny+1))     =   (Vx(ind+(ny+1)) - Vx(ind))/dx;
Exx(1:ny+1)         =   Exx(ny+2:(ny+1)*2);
Exx(end-ny:end)     =   Exx(end-(ny+1)-ny:end-(ny+1));

ind = find(eta_s~=0);
Exy(ind)            =   0.5*((Vx(ind+1) - Vx(ind))/dy + (Vy(ind+(ny+1)) - Vy(ind))/dx);
w(ind)              =   0.5*((Vy(ind+(ny+1)) - Vy(ind))/dx - (Vx(ind+1) - Vx(ind))/dy);

% Delta tao
dTxx   =   (2*eta_p .* Exx - Toxx) .* mu_p*Ddt./(eta_p+mu_p*Ddt);
dTxy   =   (2*eta_s .* Exy - Toxy) .* mu_s*Ddt./(eta_s+mu_s*Ddt);

% Set ghost points to zero
dTxx(1:ny+1)         =   dTxx(ny+2:(ny+1)*2);
dTxx(end-ny:end)     =   dTxx(end-(ny+1)-ny:end-(ny+1));
dTxx(1:ny+1:end-ny)  =   dTxx(2:(ny+1):end-ny+1);
dTxx(ny+1:ny+1:end)  =   dTxx(ny:ny+1:end-1);
ind = isnan(dTxy) == 1;
dTxy(ind) = 0;

% Tao
Txx     =  (dTxx + Toxx);
Txy     =  (dTxy + Toxy);

% Sigma
Sxx     =   Txx - P;
Sxy     =   Txy;

T2nd_s  =   zeros((ny+1)*(nx+1),1);
E2nd_s  =   zeros((ny+1)*(nx+1),1);

ind = find(eta_s ~= 0);
E2nd_s(ind)   =   (((Exx(ind) + Exx(ind+1) + Exx(ind+(ny+1)) + Exx(ind+(ny+1)+1))/4).^2 + Exy(ind).^2).^0.5;
T2nd_s(ind)   =   (((Txx(ind) + Txx(ind+1) + Txx(ind+(ny+1)) + Txx(ind+(ny+1)+1))/4).^2 + Txy(ind).^2).^0.5;

T2nd_p  =   zeros((ny+1)*(nx+1),1);
E2nd_p  =   zeros((ny+1)*(nx+1),1);

ind = find(eta_p ~= 0);
T2nd_p(ind)   =   (Txx(ind).^2 + ((Txy(ind) + Txy(ind-1) + Txy(ind-(ny+1)) + Txy(ind-(ny+1)-1))/4).^2).^0.5;
E2nd_p(ind)   =   (Exx(ind).^2 + ((Exy(ind) + Exy(ind-1) + Exy(ind-(ny+1)) + Exy(ind-(ny+1)-1))/4).^2).^0.5;

% Set ghost points to value
if strcmp('mirror',right_BC)==0
    P(1:ny+1)         =   P(ny+2:(ny+1)*2);
    P(end-ny:end)     =   P(end-(ny+1)-ny:end-(ny+1));
end
P(1:ny+1:end-ny)  =   P(2:(ny+1):end-ny+1);
P(ny+1:ny+1:end)  =   2.*P(ny:ny+1:end-1)-P(ny-1:ny+1:end-2);

%----------------- ADD MARKERS TO EMPTY CELLS ---------------------
% Find empty cells

% Interpolate stresses and strain rates from nodes to markers
clear j i
j = fix((xm+dx/2)/dx)+1;
i = fix((ym+dy/2)/dy)+1;

dxm =   xm - x_P(j);
dym =   ym - y_P(i);

ind_eq = (j-1).*(ny+1)+i;

% for nodal points
Exxn      =    Exx(ind_eq)'     .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
    +Exx(ind_eq+(ny+1))'        .* (dxm/dx)       .* (1 - (dym/dy)) ...
    +Exx(ind_eq+1)'             .* (1 - (dxm)/dx) .* (dym/dy) ...
    +Exx(ind_eq+(ny+1)+1)'      .* (dxm/dx)       .* (dym/dy);

dtxxm   =    dTxx(ind_eq)'      .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
    +dTxx(ind_eq+(ny+1))'       .* (dxm/dx)       .* (1 - (dym/dy)) ...
    +dTxx(ind_eq+1)'            .* (1 - (dxm)/dx) .* (dym/dy) ...
    +dTxx(ind_eq+(ny+1)+1)'     .* (dxm/dx)       .* (dym/dy);

Pm      =    P(ind_eq)'         .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
    +P(ind_eq+(ny+1))'          .* (dxm/dx)       .* (1 - (dym/dy)) ...
    +P(ind_eq+1)'               .* (1 - (dxm)/dx) .* (dym/dy) ...
    +P(ind_eq+(ny+1)+1)'        .* (dxm/dx)       .* (dym/dy);

% for basic nodes
clear j i ind_eq
j = fix((xm)/dx)+1;
i = fix((ym)/dy)+1;

dxm =   xm - x_Vx(j);
dym =   ym - y_Vy(i);

ind_eq = (j-1).*(ny+1)+i;

Exyn    =    Exy(ind_eq)'      .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
    +Exy(ind_eq+(ny+1))'       .* (dxm/dx)       .* (1 - (dym/dy)) ...
    +Exy(ind_eq+1)'            .* (1 - (dxm)/dx) .* (dym/dy) ...
    +Exy(ind_eq+(ny+1)+1)'     .* (dxm/dx)       .* (dym/dy);

E2ndn   =    E2nd_s(ind_eq)'   .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
    +E2nd_s(ind_eq+(ny+1))'    .* (dxm/dx)       .* (1 - (dym/dy)) ...
    +E2nd_s(ind_eq+1)'         .* (1 - (dxm)/dx) .* (dym/dy) ...
    +E2nd_s(ind_eq+(ny+1)+1)'  .* (dxm/dx)       .* (dym/dy);

dtxym    =    dTxy(ind_eq)'    .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
    +dTxy(ind_eq+(ny+1))'      .* (dxm/dx)       .* (1 - (dym/dy)) ...
    +dTxy(ind_eq+1)'           .* (1 - (dxm)/dx) .* (dym/dy) ...
    +dTxy(ind_eq+(ny+1)+1)'    .* (dxm/dx)       .* (dym/dy);

wm      =    w(ind_eq)'        .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
    +w(ind_eq+(ny+1))'         .* (dxm/dx)       .* (1 - (dym/dy)) ...
    +w(ind_eq+1)'              .* (1 - (dxm)/dx) .* (dym/dy) ...
    +w(ind_eq+(ny+1)+1)'       .* (dxm/dx)       .* (dym/dy);

% 2nd strain rate invariant on marker
Exxm    =   (Txxm)./(2.*eta_effm) + dtxxm./(2*Ddt.*mum);
Exym    =   (Txym)./(2.*eta_effm) + dtxym./(2*Ddt.*mum);
E2ndm   =   (Exxm.^2 + Exym.^2).^0.5;

% Set negativ marker stresses to 0 that Yield cannot be 0 !!!!
ind=Pm<0;
Pm(ind)=0;

% Correcting marker strain rates
E2ndm    =  E2ndn.*(E2ndm./E2ndn).^strain_rate_smoother;