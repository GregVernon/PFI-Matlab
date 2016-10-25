function [FIG,MOV] = VizFlowField(cCOORD,pCOORD,U,V,FIG,VIZ,MOV)
%% UNPACK VARIABLES
% Unpack pCOORD
NX = pCOORD.NX;
NY = pCOORD.NY;
xMin = pCOORD.xMin;
xMax = pCOORD.xMax;
yMin = pCOORD.yMin;
yMax = pCOORD.yMax;

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
if isempty(FIG)
    if VIZ.plotStateVar == true
        % Vorticity
        FIG.fig1 = figure('units','normalized','position',[0 0 1 1],'visible','off');
        subplot(131)
        FIG.oPlot = pcolor(xx,yy,OMEGA);
        FIG.oPlot.EdgeColor = 'none';
        axis('square');colorbar
        
        % Stream Function
        subplot(132)
        FIG.pPlot = pcolor(xx,yy,PSI);
        FIG.pPlot.EdgeColor = 'none';
        axis('square'); colorbar;
        
        % Velocity
        subplot(133); hold on;
        FIG.qPlot = quiver(xx,yy,U,V);
        [~,FIG.cPlot] = contour(xx(jj,ii),yy(jj,ii),PSI(jj,ii),100);
        axis('square');axis([xMin xMax yMin yMax]);
        
        % Movie Output
        if MOV.plotStateVar == true
            MOV.mov1 = VideoWriter('PFI_PlotState.mp4','MPEG-4');
            MOV.mov1.Quality = 100;
            open(MOV.mov1);
        end
    end
    
    if VIZ.plotFlow == true
        FIG.plotFlow.fig1 = figure('visible','off','units','normalized','outerposition',[0 0 2 2],'menubar','none');
        hold on;
        FIG.plotFlow.ax1 = gca();
        FIG.plotFlow.aVPlot = plot(bx,by,'k','linewidth',1);
        [~,FIG.plotFlow.cVPlot] = contourf(pX,pY,VEL,ceil(100),'linestyle','none');
        FIG.plotFlow.mVPlot = mesh(pX,pY,zeros(size(pX)),'edgecolor','k','facecolor','none');
        % FIG.plotFlow.qVPlot = quiver(pX,pY,pU,pV,.15,'color','k','ShowArrowHead','on');
%         axis('equal')
        colormap jet
%         colorbar
        if verLessThan('matlab','8.4')
            GROOT = groot;
            scrsz = GROOT.ScreenSize;
            set(FIG.plotFlow.ax1,'Units','pixels','position',scrsz);
            set(FIG.plotFlow.ax1,'CLim',[-max(max(abs(VEL))) max(max(abs(VEL)))]/1);
        else
            FIG.plotFlow.ax1.CLim = [-max(max(abs(VEL))) max(max(abs(VEL)))]/10;
        end
        axis([-2 16 -2 8])
        % Movie Output
        if MOV.plotFlow == true
            MOV.mov2 = VideoWriter('PFI_PlotFlow.mp4','MPEG-4');
            open(MOV.mov2);
        end
    end
else
    if VIZ.plotStateVar == true
        FIG.oPlot.CData = OMEGA;
        FIG.pPlot.CData = PSI;
        FIG.qPlot.UData = U;
        FIG.qPlot.VData = V;
        FIG.cPlot.ZData = PSI;
    end
    
    if VIZ.plotFlow == true
        FIG.plotFlow.cVPlot.ZData = VEL;
        if verLessThan('matlab','8.4')
            set(FIG.plotFlow.ax1,'CLim',[-max(max(abs(VEL))) max(max(abs(VEL)))]/1);
        else
            FIG.plotFlow.ax1.CLim = [-max(max(abs(VEL))) max(max(abs(VEL)))]/10;
        end
        % FIG.qVPlot.UData = pU;
        % FIG.qVPlot.VData = pV;
    end
end

if MOV.plotStateVar == true
    currFrame = getframe(FIG.fig1);
    if verLessThan('matlab','8.4')
        set(FIG.fig1,'Visible','off');
    else
        FIG.fig1.Visible = 'off';
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
        FIG.plotFlow.fig1.Visible = 'off';
    end
    writeVideo(MOV.mov2,currFrame);
end
drawnow