%======================================================================
% delete markers outside the grid
%     if timestep > 1
erase = find(xm < x(1) | xm > x(end) | ym < y(1) | ym > y(end));
Im(erase)          =  -1;
xm(erase)          =   0;
ym(erase)          =   0;
etam(erase)        =   1;
mum(erase)         =   1;
wm(erase)          =   0;
rhom0(erase)        =   1;
phim(erase)        =   1;
Cm(erase)          =   1;
Txxm(erase)        =   1;
Txym(erase)        =   1;
txxm(erase)        =   1;
txym(erase)        =   1;
dtxxm(erase)        =   1;
dtxym(erase)        =   1;
Exxm(erase)        =   1;
%         Exxn(erase)        =   0;
%         Exyn(erase)        =   0;
E2ndm(erase)       =   1;
E2ndn(erase)       =   1;
Exym(erase)        =   1;
Pm(erase)          =   1;
eta_effm(erase)    =   1;
alpha(erase)       =   1;
lambdam(erase)     =   0;
strainm(erase)     =   0;
strainvm(erase)     =   0;
workm(erase)    =   0;
if Temperature == 1
    if ra_heat == 1
        hrm(erase)         =   0;
    end
    if adiab_heat == 1 && TPdep_dens == 0
        tam(erase)         =   0;
    end
    if TPdep_dens == 1
        tam(erase)         =   0;
        tbm(erase)         =   0;
    end
    Tm(erase)          =   273;
    km(erase)          =   1;
    cpm(erase)         =   1;
    rhocpm0(erase)      =   1;
end
if Powerlaw == 1
    nm(erase)          =   1;
    Qm(erase)          =   1;
end
if Diffusion_creep == 1
    dm(erase)     =   0;
    ddm(erase)    =   0;
    straindm(erase)     =   0;
    mm(erase)     =   0;
    eta_diffm(erase)     =   0;
    def_modem(erase)     =   0;
    etadm(erase)     =   0;
    Qdm(erase)      =   1;
end

%     end

% marknum = numel(xm);
