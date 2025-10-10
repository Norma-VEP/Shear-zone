% SUBGRID DIFFUSION OF STRESSES -- START

dtoxxn = zeros((ny+1)*(nx+1),1);
dtoxyn = zeros((ny+1)*(nx+1),1);
wts    = zeros((ny+1)*(nx+1),1);
wtp    = zeros((ny+1)*(nx+1),1);

dsubgrid = 1;

sdm     = eta_effm./mum;
sdiff   = -dsubgrid*dt./sdm;
ind=sdiff<-30;
sdiff(ind)=-30;

sdiff = 1-exp(sdiff);

j = fix((xm+dx/2)/dx)+1;
i = fix((ym+dy/2)/dy)+1;

dxm =   xm - x_P(j);
dym =   ym - y_P(i);

ind_eq = (j-1).*(ny+1)+i;

toxxm    =    Toxx(ind_eq)'     .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
    +Toxx(ind_eq+(ny+1))'       .* (dxm/dx)       .* (1 - (dym/dy)) ...
    +Toxx(ind_eq+1)'            .* (1 - (dxm)/dx) .* (dym/dy) ...
    +Toxx(ind_eq+(ny+1)+1)'     .* (dxm/dx)       .* (dym/dy);

% for basic nodes
clear j i ind_eq
j = fix((xm)/dx)+1;
i = fix((ym)/dy)+1;

dxm =   xm - x_Vx(j);
dym =   ym - y_Vy(i);

ind_eq = (j-1).*(ny+1)+i;

toxym   =    Toxy(ind_eq)'      .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
    +Toxy(ind_eq+(ny+1))'       .* (dxm/dx)       .* (1 - (dym/dy)) ...
    +Toxy(ind_eq+1)'            .* (1 - (dxm)/dx) .* (dym/dy) ...
    +Toxy(ind_eq+(ny+1)+1)'     .* (dxm/dx)       .* (dym/dy);

% Update stresses according to subgrid part
dtoxxm  = toxxm-Txxm;
dtoxxm  = dtoxxm.*sdiff;
Txxm    = Txxm+dtoxxm;

dtoxym  = toxym-Txym;
dtoxym  = dtoxym.*sdiff;
Txym    = Txym+dtoxym;

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
    dtoxyn(ind_move(fil))   =   dtoxyn(ind_move(fil))+(dtoxym(fil).*wtijm(fil));
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
    dtoxxn(ind_move(fil))   =   dtoxxn(ind_move(fil))+dtoxxm(fil)*wtijm(fil);
    wtp(ind_move(fil))      =   wtp(ind_move(fil)) + wtijm(fil);
end

clear wtijm

% Divide by weights (i,j)
ind=wts>0;
dtoxyn(ind) = dtoxyn(ind)./wts(ind);

ind=wtp>0;
dtoxxn(ind) = dtoxxn(ind)./wtp(ind);


dTxx = dTxx - dtoxxn;
dTxy = dTxy - dtoxyn;

% SUBGRID DIFFUSION OF STRESSES -- END

% Updating stresses for markers
% Interpolate stresses and strain rates from nodes to markers
clear j i
j = fix((xm+dx/2)/dx)+1;
i = fix((ym+dy/2)/dy)+1;

dxm =   xm - x_P(j);
dym =   ym - y_P(i);

ind_eq = (j-1).*(ny+1)+i;

% for nodal points
dtxxm   =    dTxx(ind_eq)'      .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
    +dTxx(ind_eq+(ny+1))'       .* (dxm/dx)       .* (1 - (dym/dy)) ...
    +dTxx(ind_eq+1)'            .* (1 - (dxm)/dx) .* (dym/dy) ...
    +dTxx(ind_eq+(ny+1)+1)'     .* (dxm/dx)       .* (dym/dy);

% for basic nodes
clear j i ind_eq
j = fix((xm)/dx)+1;
i = fix((ym)/dy)+1;

dxm =   xm - x_Vx(j);
dym =   ym - y_Vy(i);

ind_eq = (j-1).*(ny+1)+i;

dtxym    =    dTxy(ind_eq)'    .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
    +dTxy(ind_eq+(ny+1))'      .* (dxm/dx)       .* (1 - (dym/dy)) ...
    +dTxy(ind_eq+1)'           .* (1 - (dxm)/dx) .* (dym/dy) ...
    +dTxy(ind_eq+(ny+1)+1)'    .* (dxm/dx)       .* (dym/dy);


% 
% Txxm = Txxm+dtxxm;
% Txym = Txym+dtxym;

