
xmABCD1   =   surface_x;
ymABCD1   =   surface_y;

% define indices for Vx nodes
j = fix(xmABCD1/dx)+1;
i = fix((ymABCD1+dy/2)/dy)+1;

ind = j<1;
j(ind) = 1;
ind = j>nx;
j(ind) = nx;
ind = i<1;
i(ind) = 1;
ind = i>ny;
i(ind) = ny;

ind_move = (j-1).*(ny+1)+i;

% define weights for 4 Vx nodes around the marker
wtijm        =   (1-(xmABCD1-x_Vx(j))/dx).*(1-(ymABCD1-y_Vx(i))/dy);
wtij1m       =   ((xmABCD1-x_Vx(j))/dx).*(1-(ymABCD1-y_Vx(i))/dy);
wti1jm       =   (1-(xmABCD1-x_Vx(j))/dx).*((ymABCD1-y_Vx(i))/dy);
wti1j1m      =   ((xmABCD1-x_Vx(j))/dx).*((ymABCD1-y_Vx(i))/dy);

% define Vx for the marker
VxmABCD1      =   Vx(ind_move)'.*wtijm +Vx(ind_move+(ny+1))'.*wtij1m +Vx(ind_move+1)'.*wti1jm +Vx(ind_move+(ny+1)+1)'.*wti1j1m;

% define indices for Vy nodes
j = fix((xmABCD1+dx/2)/dx)+1;
i = fix(ymABCD1/dy)+1;

ind = j<1;
j(ind) = 1;
ind = j>nx;
j(ind) = nx;
ind = i<1;
i(ind) = 1;
ind = i>ny;
i(ind) = ny;

ind_move = (j-1).*(ny+1)+i;

% define weights for 4 Vy nodes around the marker
wtijm        =   (1-(xmABCD1-x_Vy(j))/dx).*(1-(ymABCD1-y_Vy(i))/dy);
wtij1m       =   ((xmABCD1-x_Vy(j))/dx).*(1-(ymABCD1-y_Vy(i))/dy);
wti1jm       =   (1-(xmABCD1-x_Vy(j))/dx).*((ymABCD1-y_Vy(i))/dy);
wti1j1m      =   ((xmABCD1-x_Vy(j))/dx).*((ymABCD1-y_Vy(i))/dy);

% define Vy for the marker
VymABCD1      =   Vy(ind_move)'.*wtijm +Vy(ind_move+(ny+1))'.*wtij1m +Vy(ind_move+1)'.*wti1jm +Vy(ind_move+(ny+1)+1)'.*wti1j1m;

% MOVE MARKER
xmABCD2 = xmABCD1 + VxmABCD1*dt/2;
ymABCD2 = ymABCD1 + VymABCD1*dt/2;

% ~~~~~~~~~~~~~ 2nd ~~~~~~~~~~~~~~~~~
% define indices for Vx nodes
j = fix(xmABCD2/dx)+1;
i = fix((ymABCD2+dy/2)/dy)+1;

ind = j<1;
j(ind) = 1;
ind = j>nx;
j(ind) = nx;
ind = i<1;
i(ind) = 1;
ind = i>ny;
i(ind) = ny;

ind_move = (j-1).*(ny+1)+i;

% define weights for 4 Vx nodes around the marker
wtijm        =   (1-(xmABCD2-x_Vx(j))/dx).*(1-(ymABCD2-y_Vx(i))/dy);
wtij1m       =   ((xmABCD2-x_Vx(j))/dx).*(1-(ymABCD2-y_Vx(i))/dy);
wti1jm       =   (1-(xmABCD2-x_Vx(j))/dx).*((ymABCD2-y_Vx(i))/dy);
wti1j1m      =   ((xmABCD2-x_Vx(j))/dx).*((ymABCD2-y_Vx(i))/dy);

% define Vx for the marker
VxmABCD2      =   Vx(ind_move)'.*wtijm +Vx(ind_move+(ny+1))'.*wtij1m +Vx(ind_move+1)'.*wti1jm +Vx(ind_move+(ny+1)+1)'.*wti1j1m;

% define indices for Vy nodes
j = fix((xmABCD2+dx/2)/dx)+1;
i = fix(ymABCD2/dy)+1;

ind = j<1;
j(ind) = 1;
ind = j>nx;
j(ind) = nx;
ind = i<1;
i(ind) = 1;
ind = i>ny;
i(ind) = ny;

ind_move = (j-1).*(ny+1)+i;

% define weights for 4 Vy nodes around the marker
wtijm        =   (1-(xmABCD2-x_Vy(j))/dx).*(1-(ymABCD2-y_Vy(i))/dy);
wtij1m       =   ((xmABCD2-x_Vy(j))/dx).*(1-(ymABCD2-y_Vy(i))/dy);
wti1jm       =   (1-(xmABCD2-x_Vy(j))/dx).*((ymABCD2-y_Vy(i))/dy);
wti1j1m      =   ((xmABCD2-x_Vy(j))/dx).*((ymABCD2-y_Vy(i))/dy);

% define Vy for the marker
VymABCD2      =   Vy(ind_move)'.*wtijm +Vy(ind_move+(ny+1))'.*wtij1m +Vy(ind_move+1)'.*wti1jm +Vy(ind_move+(ny+1)+1)'.*wti1j1m;

% MOVE MARKER
xmABCD3 = xmABCD1 + VxmABCD2*dt/2;
ymABCD3 = ymABCD1 + VymABCD2*dt/2;

% ~~~~~~~~~~~~~ 3rd ~~~~~~~~~~~~~~~~~
% define indices for Vx nodes
j = fix(xmABCD3/dx)+1;
i = fix((ymABCD3+dy/2)/dy)+1;

ind = j<1;
j(ind) = 1;
ind = j>nx;
j(ind) = nx;
ind = i<1;
i(ind) = 1;
ind = i>ny;
i(ind) = ny;

ind_move = (j-1).*(ny+1)+i;

% define weights for 4 Vx nodes around the marker
wtijm        =   (1-(xmABCD3-x_Vx(j))/dx).*(1-(ymABCD3-y_Vx(i))/dy);
wtij1m       =   ((xmABCD3-x_Vx(j))/dx).*(1-(ymABCD3-y_Vx(i))/dy);
wti1jm       =   (1-(xmABCD3-x_Vx(j))/dx).*((ymABCD3-y_Vx(i))/dy);
wti1j1m      =   ((xmABCD3-x_Vx(j))/dx).*((ymABCD3-y_Vx(i))/dy);

% define Vx for the marker
VxmABCD3      =   Vx(ind_move)'.*wtijm +Vx(ind_move+(ny+1))'.*wtij1m +Vx(ind_move+1)'.*wti1jm +Vx(ind_move+(ny+1)+1)'.*wti1j1m;

% define indices for Vy nodes
j = fix((xmABCD3+dx/2)/dx)+1;
i = fix(ymABCD3/dy)+1;

ind = j<1;
j(ind) = 1;
ind = j>nx;
j(ind) = nx;
ind = i<1;
i(ind) = 1;
ind = i>ny;
i(ind) = ny;

ind_move = (j-1).*(ny+1)+i;

% define weights for 4 Vy nodes around the marker
wtijm        =   (1-(xmABCD3-x_Vy(j))/dx).*(1-(ymABCD3-y_Vy(i))/dy);
wtij1m       =   ((xmABCD3-x_Vy(j))/dx).*(1-(ymABCD3-y_Vy(i))/dy);
wti1jm       =   (1-(xmABCD3-x_Vy(j))/dx).*((ymABCD3-y_Vy(i))/dy);
wti1j1m      =   ((xmABCD3-x_Vy(j))/dx).*((ymABCD3-y_Vy(i))/dy);

% define Vy for the marker
VymABCD3      =   Vy(ind_move)'.*wtijm +Vy(ind_move+(ny+1))'.*wtij1m +Vy(ind_move+1)'.*wti1jm +Vy(ind_move+(ny+1)+1)'.*wti1j1m;

% MOVE MARKER
xmABCD4 = xmABCD1 + VxmABCD3*dt;
ymABCD4 = ymABCD1 + VymABCD3*dt;

% ~~~~~~~~~~~~~ 4th ~~~~~~~~~~~~~~~~~
% define indices for Vx nodes
j = fix(xmABCD4/dx)+1;
i = fix((ymABCD4+dy/2)/dy)+1;

ind = j<1;
j(ind) = 1;
ind = j>nx;
j(ind) = nx;
ind = i<1;
i(ind) = 1;
ind = i>ny;
i(ind) = ny;

ind_move = (j-1).*(ny+1)+i;

% define weights for 4 Vx nodes around the marker
wtijm        =   (1-(xmABCD4-x_Vx(j))/dx).*(1-(ymABCD4-y_Vx(i))/dy);
wtij1m       =   ((xmABCD4-x_Vx(j))/dx).*(1-(ymABCD4-y_Vx(i))/dy);
wti1jm       =   (1-(xmABCD4-x_Vx(j))/dx).*((ymABCD4-y_Vx(i))/dy);
wti1j1m      =   ((xmABCD4-x_Vx(j))/dx).*((ymABCD4-y_Vx(i))/dy);

% define Vx for the marker
VxmABCD4      =   Vx(ind_move)'.*wtijm +Vx(ind_move+(ny+1))'.*wtij1m +Vx(ind_move+1)'.*wti1jm +Vx(ind_move+(ny+1)+1)'.*wti1j1m;

% define indices for Vy nodes
j = fix((xmABCD4+dx/2)/dx)+1;
i = fix(ymABCD4/dy)+1;

ind = j<1;
j(ind) = 1;
ind = j>nx;
j(ind) = nx;
ind = i<1;
i(ind) = 1;
ind = i>ny;
i(ind) = ny;

ind_move = (j-1).*(ny+1)+i;

% define weights for 4 Vy nodes around the marker
wtijm        =   (1-(xmABCD4-x_Vy(j))/dx).*(1-(ymABCD4-y_Vy(i))/dy);
wtij1m       =   ((xmABCD4-x_Vy(j))/dx).*(1-(ymABCD4-y_Vy(i))/dy);
wti1jm       =   (1-(xmABCD4-x_Vy(j))/dx).*((ymABCD4-y_Vy(i))/dy);
wti1j1m      =   ((xmABCD4-x_Vy(j))/dx).*((ymABCD4-y_Vy(i))/dy);

% define Vy for the marker
VymABCD4      =   Vy(ind_move)'.*wtijm +Vy(ind_move+(ny+1))'.*wtij1m +Vy(ind_move+1)'.*wti1jm +Vy(ind_move+(ny+1)+1)'.*wti1j1m;


VxmS     =   (1/6)*(VxmABCD1+2*VxmABCD2+2*VxmABCD3+VxmABCD4);
VymS     =   (1/6)*(VymABCD1+2*VymABCD2+2*VymABCD3+VymABCD4);

surface_x_new = surface_x + VxmS*dt;
surface_y_new = surface_y + VymS*dt;

surface_y = spline(surface_x_new,surface_y_new,surface_x);
surface_y(1,:)=smooth(surface_y,surface_smoother);

    %%%%%%%%%%%%%%%%%%%%%
% figure(8),    plot(surface_x,surface_y,'r'), hold on

