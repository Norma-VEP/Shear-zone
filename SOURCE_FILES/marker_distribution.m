% Defining values of markers
marker = 0;

load SETUP.mat Granite_corr

for im = 1:nym
    for jm = 1:nxm
        marker  =   marker+1;
        xm(marker)  =   Lx/nxm/2 + Lx/nxm*(jm-1) + (rand-0.5)*Lx/nxm;
        ym(marker)  =   Ly/nym/2 + Ly/nym*(im-1) + (rand-0.5)*Ly/nym;
        
                
%         for iS = 1:size(SETUP,1)
%             % Rock1
%             if xm(marker) > SETUP(iS,2) && xm(marker) < SETUP(iS,3) && ym(marker) > SETUP(iS,4) && ym(marker) < SETUP(iS,5)
%                 rocktype = SETUP(iS,1);
%             end
%         end
        
        % =============== BUILD YOUR OWN GEOMETRY HERE =================
        if Granite_corr(im,jm) == 1
            rocktype = 5;
        elseif Granite_corr(im,jm) == 2
            rocktype = 3;
        elseif Granite_corr(im,jm) == 3
            rocktype = 6;
        end 
        % ==============================================================

        
        if ~isempty(find(pattern_marker==Phase(rocktype))==1)
            % Horizontal stripes
            if (round(ym(marker)/pattern_ydim)<ym(marker)/pattern_ydim) && (strcmp('horizontal',pattern_type)==1 || strcmp('kaki',pattern_type)==1)
                rocktype=rocktype+1;
            end
            % Vertical stripes            
            if (round(xm(marker)/pattern_xdim)<xm(marker)/pattern_xdim) && (strcmp('vertical',pattern_type)==1 || strcmp('kaki',pattern_type)==1)
                rocktype=rocktype+1;
            end
        end
        
        % put values to markers
        etam(marker)    =   eta(rocktype);
        mum(marker)     =   mu(rocktype);
        rhom0(marker)   =   rho(rocktype);
        Cm(marker)      =   C(rocktype);
        phim(marker)    =   phi(rocktype);
        lambdam(marker) =   lambda(rocktype);
        Im(marker)      =   Phase(rocktype);
        if Temperature == 1
            if ra_heat == 1
                hrm(marker)     =   hr(rocktype);
            end
            if adiab_heat == 1 && TPdep_dens == 0
                tam(marker)     =   ta(rocktype);
            end
            if TPdep_dens == 1
                tam(marker)     =   ta(rocktype);
                tbm(marker)     =   tb(rocktype);
            end
            Tm(marker)      =   T(rocktype);
            km(marker)      =   kk(rocktype);
            cpm(marker)     =   cp(rocktype);
            rhocpm0(marker) =   rho(rocktype)*cp(rocktype);
        end
        if Powerlaw == 1
            nm(marker)      =   n(rocktype);
            Qm(marker)      =   Q(rocktype);
        end
        if Diffusion_creep == 1
            dm(marker)          =   Gr(rocktype);   % grain size
            mm(marker)          =   m(rocktype);    % diffusion creep coefficient
            etadm(marker)       =   dCr(rocktype);  % diffusion creep viscosity
            Qdm(marker)         =   Qd(rocktype);  % diffusion creep viscosity
        end

%         if lambdam(marker)>0.8 && (xm(marker)<15e3 || xm(marker)>85e3)
%             lambdam(marker)=0.8;
%         end
    end
end

clear Granite_corr

% Set initial temperature distribution
if Temperature == 1
    for marker=1:marknum
        
        % PERTUBATION
        value1 = (sin((xm(marker)-Pert_beg)/((Pert_end-Pert_beg)/2)*pi-pi/2)+1)/2;
        if (xm(marker)<Pert_beg || xm(marker)>Pert_end)
            value1=0;
        end
        
        % DEEPER MANTLE
        Tm(marker)=273+L_A_temp+M_grad*(ym(marker)-A_thick-L_thick)/1e3;
        
        % LITHOSPHERIC MANTLE
        if ym(marker)<=A_thick+L_thick-Pert_add*value1
            Tm(marker)=273+Moho_temp+(L_A_temp-(Pert_add/1e3*value1*M_grad)-Moho_temp)/(L_thick-C_thick)*(ym(marker)-A_thick-C_thick+Pert_add*value1);
        end
        
        % CRUST
        if ym(marker)<=A_thick+C_thick
            Tm(marker)=273+Moho_temp/C_thick*(ym(marker)-A_thick);
        end
        
        % STICKY-AIR
        if ym(marker)<A_thick
            Tm(marker)=top_T;
        end
        
        %         % OCEANIC LITHOSPHERE
        %         if ym(marker)>12e3 && ym(marker)<105e3+12e3
        %             Tm_O(marker)=1400.*erf((ym(marker)-12e3)./(2*sqrt(1e-6*50e6*(SecYear))))+273;
        %         elseif ym(marker)<=12e3
        %             Tm_O(marker)=273;
        %         end
        %
        %         % SETUP
        %         if xm(marker)<=400e3
        %             Tm(marker)=Tm_O(marker);
        %         elseif xm(marker)>=500e3
        %             Tm(marker)=Tm_C(marker);
        %         end
        %         if xm(marker)>400e3 && xm(marker)<500e3
        % %             depthm=12e3-((xm(maker)-400e3)/50e3)*4e3;
        %             Tm(marker)=Tm_O(marker)+(((xm(marker)-400e3)/100e3)*(Tm_C(marker)-Tm_O(marker)));
        %         end
        
        
    end
   
          
    
end

