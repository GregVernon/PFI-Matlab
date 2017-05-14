function [FIG,MOV] = VizFlowField(cCOORD,pCOORD,PSI,OMEGA,U,V,FIG,VIZ,MOV,IMG,pCycle)
%% UNPACK VARIABLES
% Unpack pCOORD
NX = pCOORD.NX;
NY = pCOORD.NY;
xMin = pCOORD.xMin;
xMax = pCOORD.xMax;
yMin = pCOORD.yMin;
yMax = pCOORD.yMax;
xx = pCOORD.xx;
yy = pCOORD.yy;

% Unpack cCOORD
pX = cCOORD.pX;
pY = cCOORD.pY;
bx = cCOORD.bx;
by = cCOORD.by;
dmux = cCOORD.dmux;
dmuy = cCOORD.dmuy;
detx = cCOORD.detx;
dety = cCOORD.dety;
eta = cCOORD.eta;
mu = cCOORD.mu;
%% COMPUTE CARTESIAN VELOCITY
[pU,pV,VEL] = PARA_VEL(U,V,cCOORD,pCOORD);
%%
OMEGA = gather(OMEGA);
PSI = gather(PSI);
U = gather(U);
V = gather(V);
VEL = gather(VEL);
%%
if isempty(FIG)
    %%
    if VIZ.plotStateVar == true
        % Vorticity
        FIG.plotState.fig1 = figure('units','normalized','position',[0 0 1 1],'visible','on');
        FIG.plotState.oPlot = pcolor(xx,yy,OMEGA);
        FIG.plotState.oPlot.FaceColor = 'interp';
        FIG.plotState.oPlot.EdgeColor = 'none';
        axis('square');colorbar
        title('OMEGA (Vorticity)')
        
        
        % Stream Function
        FIG.plotState.fig2 = figure('units','normalized','position',[0 0 1 1],'visible','on');
        FIG.plotState.pPlot = pcolor(xx,yy,PSI);
        FIG.plotState.pPlot.FaceColor = 'interp';
        FIG.plotState.pPlot.EdgeColor = 'none';
        axis('square'); colorbar;
        title('PSI (Stream-Function)')
        
        % Velocity
        FIG.plotState.fig3 = figure('units','normalized','position',[0 0 1 1],'visible','on');
        hold on;
        FIG.plotState.qvPlot = quiver(xx,yy,U,V);
        [FIG.plotState.pvPlot] = pcolor(xx,yy,VEL);
        FIG.plotState.pvPlot.FaceColor = 'interp';
        FIG.plotState.pvPlot.EdgeColor = 'none';
        axis('square');axis([xMin xMax yMin yMax]);
        title('VEL (Velocity)')
        % Movie Output
        if MOV.plotStateVar == true
            MOV.mov1 = VideoWriter('PFI_PlotState.mp4','MPEG-4');
            MOV.mov1.Quality = 100;
            open(MOV.mov1);
        end
        % Image Output
    end
    %%
    if VIZ.plotFlow == true
        FIG.plotFlow.fig1 = figure('visible','off','units','normalized','outerposition',[0 0 2 2],'menubar','none');
        hold on;
        inView = (pX >= -2 & pX <=16 & pY >= -9 & pY <=9);
        [inVrow,inVcol] = find(inView);
        Vrow = (min(inVrow):max(inVrow)+1);
        Vcol = (min(inVcol)-1:max(inVcol)+1);
        FIG.plotFlow.ax1 = gca();
        FIG.plotFlow.aVPlot = plot(bx,by,'k','linewidth',1);
        %         [~,FIG.plotFlow.cVPlot] = contourf(pX,pY,VEL,ceil(100),'linestyle','none');
        [FIG.plotFlow.cVPlot] = pcolor(pX(Vrow,Vcol),pY(Vrow,Vcol),VEL(Vrow,Vcol));
        [~,FIG.plotFlow.cPlot] = contour(pX(Vrow,Vcol),pY(Vrow,Vcol),PSI(Vrow,Vcol),'-k','LevelList',sort([linspace(min(min(PSI(inView))),max(max(PSI(inView))),20) 0]));
        %         FIG.plotFlow.mVPlot = mesh(pX,pY,zeros(size(pX)),'edgecolor','k','facecolor','none');
        % FIG.plotFlow.qVPlot = quiver(pX,pY,pU,pV,.15,'color','k','ShowArrowHead','on');
        %         axis('equal')
        %         colormap jet
        %         colorbar
        if verLessThan('matlab','8.4')
            GROOT = groot;
            scrsz = GROOT.ScreenSize;
            set(FIG.plotFlow.ax1,'Units','pixels','position',scrsz);
            set(FIG.plotFlow.ax1,'CLim',[-max(max(abs(VEL))) max(max(abs(VEL)))]/1);
            set(FIG.plotFlow.cVPlot,'FaceColor','interp');
            set(FIG.plotFlow.cVPlot,'EdgeColor','none');
        else
            FIG.plotFlow.ax1.CLim = [min(min((VEL(inView)))) max(max((VEL(inView))))]/1;
            FIG.plotFlow.cVPlot.FaceColor = 'interp';
            FIG.plotFlow.cVPlot.EdgeColor = 'none';
        end
        axis([-2 16 -9 9])
        % Movie Output
        if MOV.plotFlow == true
            MOV.mov2 = VideoWriter('PFI_PlotFlow_2.mp4','MPEG-4');
            open(MOV.mov2);
        end
    end
else
    if VIZ.plotStateVar == true
        FIG.plotState.oPlot.CData = OMEGA;
        FIG.plotState.pPlot.CData = PSI;
        FIG.plotState.qvPlot.UData = U;
        FIG.plotState.qvPlot.VData = V;
        FIG.plotState.pvPlot.CData = VEL;
    end
    
    if VIZ.plotFlow == true
        inView = (pX >= -2 & pX <=16 & pY >= -9 & pY <=9);
        [inVrow,inVcol] = find(inView);
        Vrow = (min(inVrow):max(inVrow)+1);
        Vcol = (min(inVcol)-1:max(inVcol)+1);
        FIG.plotFlow.cPlot.ZData = PSI(Vrow,Vcol);
        FIG.plotFLOW.cPlot.LevelList = sort([linspace(min(min(PSI(inView))),max(max(PSI(inView))),20) 0]);
        FIG.plotFlow.cVPlot.CData = VEL(Vrow,Vcol);
        if verLessThan('matlab','8.4')
            set(FIG.plotFlow.ax1,'CLim',[-max(max(abs(OMEGA))) max(max(abs(OMEGA)))]/1);
        else
            FIG.plotFlow.ax1.CLim = [min(min((VEL(inView)))) max(max((VEL(inView))))]/1;
        end
        % FIG.qVPlot.UData = pU;
        % FIG.qVPlot.VData = pV;
    end
end

if MOV.plotStateVar == true
    currFrame = getframe(FIG.fig1);
    if verLessThan('matlab','8.4')
        set(FIG.plotState.fig1,'Visible','off');
    else
        FIG.plotState.fig1.Visible = 'on';
    end
    writeVideo(MOV.mov1,currFrame);
end

if MOV.plotFlow == true
    if verLessThan('matlab','8.4')
        GROOT = groot;
        scrsz = GROOT.ScreenSize;
        set(FIG.plotFlow.ax1,'Units','pixels','position',scrsz);
        currFrame = getframe(FIG.plotFlow.ax1);
        %         set(FIG.plotFlow.fig1,'Visible','off');
    else
        currFrame = getframe(FIG.plotFlow.ax1);
        FIG.plotFlow.fig1.Visible = 'on';
    end
    writeVideo(MOV.mov2,currFrame);
end

% Image Output
if IMG.plotFlow == true
    print(['FlowPlot_' num2str(pCycle)],'-dpng','-r300');
end
drawnow