%% Visualisation VISJON-LP

clear all%, clf

mkdir FIGURES

cd COLORMAPS
    mat = dir('*.mat');
    for i=1:length(mat)
        load(mat(i).name);
    end    
cd ..

figure('position',[300 300 850 650])

for num = [100:100:1000000]
    
    cd OUTPUT
    if exist(['Shearing_',num2str(1e6+num),'.mat'])
        load (['Shearing_',num2str(1e6+num),'.mat']);
        cd ..
    else
        cd ..
        break;
    end
    
 clf, subplot(221)
    
        for i=1:14
            Phase1  =   find(Im == i);
            marksize = 0.1;
            hold on
            plot(xm(Phase1).*100,ym(Phase1).*100,'.','Color',Colormaps(i,:),'MarkerSize',marksize)
        end
    shading interp
    axis image, axis ij
    colorbar
    colormap(subplot(221),gray(3))
    set(colorbar,'YTick',0:0.5:1,'YTickLabel',{'Biotite','Quartz','Plagioclase'});
    title(['Composition'])
%     plot(x2d_b,y2d_b,'k',x2d_b',y2d_b','k')
%     quiver(x2d_b,y2d_b,Vx2db,Vy2db,'r','LineWidth',1.5)
    hold on
%     comp=['Himalaya_',num2str(num+1e6)];
%     print('-r150', '-djpeg', comp)
%     ind=find(rho_s<2000);
%     eta_s(ind)=NaN;
    E2nd_s_2d   =   reshape(E2nd_s,ny+1,nx+1);
    eta_s_2d    =   reshape(eta_s,ny+1,nx+1);
    T2nd_s_2d   =   reshape(T2nd_s,ny+1,nx+1);
    rho_s_2d    =   reshape(rho_s,ny+1,nx+1);
%     Temp_2d     =   reshape(Temp,ny+1,nx+1);
    strain_2d   =   reshape(strain,ny+1,nx+1);
    D2d        =   reshape(D_grain,ny+1,nx+1);
    def_mode_2d        =   reshape(def_mode,ny+1,nx+1);
        
    plot([0 5 5 0 0],[0 0 5 5 0],'k')
    xlabel('x (cm)')
    ylabel('y (cm)')

    
%     contour(x,y,(Temp_2d(1:end-1,1:end-1)-273),[100:200:1300],'w','LineWidth',1)
%     plot([0 0 1000e3 1000e3 0],[0 670e3 670e3 0 0],'k-')
%     axis([0 1e6 0 670e3])
    
    hold off
    
    subplot(222)
%     figure(2), clf
        colormap(subplot(222),flip(batlow,1))
    
%     ind = find(rho_s_2d < 200);
%     E2nd_s_2d(ind)  = 1e-14;   
%     T2nd_s_2d(ind)  = 1e20;   
%     eta_s_2d(ind)   = 1e21;   
    
%         subplot(311)
    pcolor(x.*100,y.*100,log10(eta_s_2d(1:end-1,1:end-1))), caxis([16 23])
            colormap(subplot(222),flip(batlow,1))

    shading interp
    colorbar
    ylabel(colorbar,'log_{10} \eta (Pa.s)','FontSize',12,'Rotation',90);
    hColourbar.Label.Position(1) = 3;
    axis image, axis ij
    title('Viscosity')
    hold on
    contour(x.*100,y.*100,((rho_s_2d(1:end-1,1:end-1))),[2750 2900],'k')
%     fill([0 0 500e3 500e3 0],[0 10e3 10e3 0 0],'w','EdgeColor','w'), hold on
%     plot([0 0 1000e3 1000e3 0],[0 670e3 670e3 0 0],'k-')
        plot([0 5 5 0 0],[0 0 5 5 0],'k')
    xlabel('x (cm)')
    ylabel('y (cm)')

caxis([16 23])
%     fill([0 x Lx 0],[0 surface_y(1:5:end) 0 0],'w','EdgeColor','w','FaceColor','w','LineWidth',0.001)

%     
%     subplot(223)
%     pcolor(x,y,T2nd_s_2d(1:end-1,1:end-1))%, caxis([18 24])
%     shading interp
%     colorbar
%     axis image, axis ij
%     title('\tau [Pa]')
%     hold on
% %     fill([0 x Lx 0],[0 surface_y(1:5:end) 0 0],'w','EdgeColor','w','FaceColor','w','LineWidth',0.001)
% %    
% %         subplot(132)
%     figure(2), clf
    
%     ind = find(rho_s_2d < 200);
%     E2nd_s_2d(ind)  = 1e-14;   
%     T2nd_s_2d(ind)  = 1e20;   
%     eta_s_2d(ind)   = 1e21;   
    
       subplot(223)
    pcolor(x.*100,y.*100,log10(D2d(1:end-1,1:end-1))), caxis([-6 -3])
        colormap(subplot(223),flip(davos,1))

    shading interp
    colorbar
    ylabel(colorbar,'log_{10} d (m)','FontSize',12,'Rotation',90);
    hColourbar.Label.Position(1) = 3;
    axis image, axis ij
    title('Grain size')
    axis image, axis ij
    hold on
    contour(x.*100,y.*100,((rho_s_2d(1:end-1,1:end-1))),[2750 2900],'k')
%     fill([0 0 500e3 500e3 0],[0 10e3 10e3 0 0],'w','EdgeColor','w'), hold on
%     plot([0 0 1000e3 1000e3 0],[0 670e3 670e3 0 0],'k-')
    plot([0 5 5 0 0],[0 0 5 5 0],'k')
    xlabel('x (cm)')
    ylabel('y (cm)')
caxis([-6 -3])
%     subplot(224)
%     pcolor(x,y,log10(E2nd_s_2d(1:end-1,1:end-1)))%, caxis([-17 -12])
%     shading interp
%     colorbar
%     axis image, axis ij
%     title('\dotepsilon [1/s]')
%     hold on
% %     fill([0 x Lx 0],[0 surface_y(1:5:end) 0 0],'w','EdgeColor','w','FaceColor','w','LineWidth',0.001)
% %      
subplot(224)
    pcolor(x.*100,y.*100,(def_mode_2d(1:end-1,1:end-1))), caxis([-1 1])
        colormap(subplot(224),flip(hawaii,1))

    shading interp
    colorbar
%     ylabel(colorbar,'B  ','FontSize',12,'Rotation',90);
    set(colorbar,'YTick',-1:1:1,'YTickLabel',{'Brittle','Disl. creep','Diff. creep',});
    axis image, axis ij
    title('Deformation mechanism')
    hold on
    contour(x.*100,y.*100,((rho_s_2d(1:end-1,1:end-1))),[2750 2900],'k')
%     fill([0 0 500e3 500e3 0],[0 10e3 10e3 0 0],'w','EdgeColor','w'), hold on
%     plot([0 0 1000e3 1000e3 0],[0 670e3 670e3 0 0],'k-')
    plot([0 5 5 0 0],[0 0 5 5 0],'k')
    xlabel('x (cm)')
    ylabel('y (cm)')
caxis([-1 1])

sgtitle(['Time: ',num2str((time)/SecYear/1e6),' Ma'])

%     
cd FIGURES
comp=['CompVisc_vid_',num2str(num+1e6)];
print('-r150', '-djpeg', comp)
% %
cd ..

end