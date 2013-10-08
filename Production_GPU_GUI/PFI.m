function PFI(SimType)
close all
clc
if SimType
    % Continue from file
    [Fname,Fpath]=uigetfile('*.*');
    DATA=load(strcat(Fpath,Fname));
    xwidth=DATA(end-1,1);
    yheight=DATA(end-1,2);
    k=DATA(end-1,3);
    ibl=DATA(end-1,5);
    re=DATA(end,1);
    Dt=DATA(end,6);
    tol=DATA(end,7);
    xmin=DATA(end,8);
    xmax=DATA(end,9);
    ymin=DATA(end,10);
    ymax=DATA(end,11);
    a=DATA(end,12);
else
    Fname='';
    Fpath='';
    xwidth=200;
    yheight=400;
    ibl=20;
    re=700;
    tol=1e-5;
    xmin=-20;
    xmax=20;
    ymin=1;
    ymax=11;
    a=0;
    Dt=delta_t(xwidth,yheight,xmin,xmax,ymin,ymax,re);
end
hf=figure('MenuBar','none','Name','Parabola Flow Interactive','Position',...
    [0,0,1024,700],'Units','pixels','NumberTitle','off','Resize','off'...
    ,'Color','w','UserData',{Fpath,Fname});
movegui('center')
GridProp=uipanel(hf,'Title','Grid Properties','FontSize',12,'Units',...
    'pixels','Position',[10,590,770/3,110],'TitlePosition','CenterTop',...
    'BorderType','none','BackgroundColor',[.8,.8,.8]);
FluidProp=uipanel(hf,'Title','Fluid Properties','FontSize',12,'Units',...
    'pixels','Position',[20+770/3,590,770/3,110],'TitlePosition',...
    'CenterTop','BorderType','none','BackgroundColor',[.8,.8,.8]);
FileProp=uipanel(hf,'Title','File Settings','FontSize',12,'Units',...
    'pixels','Position',[30+2*770/3,590,770/3,110],'TitlePosition',...
    'CenterTop','BorderType','none','BackgroundColor',[.8,.8,.8]);

Xmin=uicontrol(GridProp,'Style','PushButton','Units','pixels',...
    'Position',[10,70,100,20],'String',{strcat('Xmin = ',num2str(xmin))},'FontSize',8,...
    'UserData',xmin,'CallBack',@XminCallBack,'BackgroundColor','w');
Xmax=uicontrol(GridProp,'Style','PushButton','Units','pixels',...
    'Position',[10,50,100,20],'String',{strcat('XMax = ',num2str(xmax))},'FontSize',8,...
    'UserData',xmax,'CallBack',@XmaxCallBack,'BackgroundColor','w');
Ymin=uicontrol(GridProp,'Style','PushButton','Units','pixels',...
    'Position',[10,30,100,20],'String',{strcat('Ymin = ',num2str(ymin))},'FontSize',8,...
    'UserData',ymin,'CallBack',@YminCallBack,'BackgroundColor','w');
Ymax=uicontrol(GridProp,'Style','PushButton','Units','pixels',...
    'Position',[10,10,100,20],'String',{strcat('Ymax = ',num2str(ymax))},'FontSize',8,...
    'UserData',ymax,'CallBack',@YmaxCallBack,'BackgroundColor','w');
Nx=uicontrol(GridProp,'Style','PushButton','Units','pixels','Position',...
    [770/3-110,70,100,20],'String',{strcat('Nx = ',num2str(xwidth))},'FontSize',8,'UserData',...
    xwidth,'CallBack',@NxCallBack,'BackgroundColor','w');
My=uicontrol(GridProp,'Style','PushButton','Units','pixels','Position',...
    [770/3-110,50,100,20],'String',{strcat('My = ',num2str(yheight))},'FontSize',8,'UserData',...
    yheight,'CallBack',@MyCallBack,'BackgroundColor','w');
Xskew=uicontrol(GridProp,'Style','PushButton','Units','pixels',...
    'Position',[770/3-110,30,100,20],'String','Xskew = 1','FontSize',...
    8,'UserData',1,'CallBack',@XskewCallBack,'BackgroundColor','w');
Yskew=uicontrol(GridProp,'Style','PushButton','Units','pixels',...
    'Position',[770/3-110,10,100,20],'String','Yskew = 1','FontSize',8,...
    'UserData',1,'CallBack',@YskewCallBack,'BackgroundColor','w');

A=uicontrol(FluidProp,'Style','PushButton','Units','pixels','Position',...
    [10,70,100,20],'String',{strcat('A = ',num2str(a))},'FontSize',8,'UserData',a,...
    'CallBack',@ACallBack,'BackgroundColor','w');
Re=uicontrol(FluidProp,'Style','PushButton','Units','pixels','Position',...
    [10,50,770/3-20,20],'String',{strcat('Re = ',num2str(re))},'FontSize',8,'UserData',re,...
    'CallBack',@ReCallBack,'BackgroundColor','w');


dt=uicontrol(FluidProp,'Style','PushButton','Units','pixels','Position',...
    [10,30,770/3-20,20],'String',{strcat('Dt = ',num2str(Dt))},'FontSize',...
    8,'UserData',Dt,'CallBack',@dtCallBack,'BackgroundColor','w');
Tol=uicontrol(FluidProp,'Style','PushButton','Units','pixels',...
    'Position',[10,10,770/3-20,20],'String',{strcat('Tolerance = ',num2str(tol))},'FontSize',8,...
    'UserData',tol,'CallBack',@TolCallBack,'BackgroundColor','w');
IBL=uicontrol(FluidProp,'Style','PushButton','Units','pixels',...
    'Position',[770/3-110,70,100,20],'String',{strcat('IBL = ',num2str(ibl))},'FontSize',8,...
    'UserData',ibl,'CallBack',@IBLCallBack,'BackgroundColor','w');

report=uicontrol(FileProp,'Style','PushButton','Units','pixels',...
    'Position',[10,70,770/3-20,20],'String','Report: 100 ITER ',...
    'FontSize',8,'UserData',100,'CallBack',@reportCallBack...
    ,'BackgroundColor','w');
PSAVE=uicontrol(FileProp,'Style','PushButton','Units','pixels',...
    'Position',[10,50,770/3-20,20],'String','FileWrite: 1000 ITER',...
    'FontSize',8,'UserData',1000,'CallBack',@PSAVECallBack...
    ,'BackgroundColor','w');
Ot=uicontrol(FileProp,'Style','PushButton','Units','pixels','Position',...
    [10,30,770/3-20,20],'String','Iterations: 1000000','FontSize',8,...
    'UserData',1000000,'CallBack',@OtCallBack,'BackgroundColor','w');
filename=uicontrol(FileProp,'Style','PushButton','Units','pixels',...
    'Position',[10,10,770/3-20,20],'String',{strcat('FileName :','pfRe7a0p00')},...
    'FontSize',8,'UserData',{'pfiRe7A0p00'},'CallBack',@filenameCallBack...
    ,'BackgroundColor','w');

AFCJET=uipanel(hf,'Title','Jet Properties','TitlePosition','CenterTop','FontSize',12,'Units','Pixels','Position',[10,488,770/3,100],'BackgroundColor',[.8,.8,.8]);
JETTEST=uicontrol(AFCJET,'Style','checkbox','String','Include Jets','Position',[20,60,80,20],'BackgroundColor',[.8,.8,.8],'CallBack',@JETTESTCallBack);
JETLIST=uicontrol(AFCJET,'Style','popupmenu','Units','Pixels','String',{'JET 1'},'Position',[140,50,50,30],'CallBack',@JETLISTCallBack,'Enable','off');
JETADD=uicontrol(AFCJET,'Style','PushButton','Units','Pixels','String','ADD JET','Position',[200,60,50,20],'CallBack',@JETADDCallBack,'Enable','off');
JETSTART=uicontrol(AFCJET,'Style','PushButton','Units','Pixels','String','Start= 110','UserData',110,'Position',[10,30,100,20],'BackgroundColor','w','CallBack',@JETSTARTCallBack,'Enable','off');
JETEND=uicontrol(AFCJET,'Style','PushButton','Units','Pixels','String','End= 114','UserData',114,'Position',[10,10,100,20],'BackgroundColor','w','CallBack',@JETENDCallBack,'Enable','off');
JETAMP=uicontrol(AFCJET,'Style','PushButton','Units','Pixels','String','Amplitude= 1','UserData',1,'Position',[145,30,100,20],'BackgroundColor','w','CallBack',@JETAMPCallBack,'Enable','off');
JETFRQ=uicontrol(AFCJET,'Style','PushButton','Units','Pixels','String','Frequency= 1','UserData',1,'Position',[145,10,100,20],'BackgroundColor','w','CallBack',@JETFRQCallBack,'Enable','off');

COMPGRID=axes('Parent',hf,'Units','Pixels','Position',[40,40,300,300],...
    'XGrid','off','YGrid','off','ZGrid','off','Box','on');
PHYSGRID=axes('Parent',hf,'Units','Pixels','Position',[684,40,300,300],...
    'XGrid','off','YGrid','off','ZGrid','off','Box','on');


[PGrid,CGrid]=GridCompute;

COMPVIEW=uibuttongroup(hf,'Units','Pixels','Position',...
    [80,0,225,20],'SelectionChangeFCN',@COMPAXIS...
    ,'BackgroundColor','w');
CalcGrid=uicontrol(COMPVIEW,'Units','Pixels','Style','Radio','String',...
    'Calculation Grid','Position',[105,2,110,15],'Value',1,...
    'BackgroundColor','w');
IndxGrid=uicontrol(COMPVIEW,'Units','Pixels','Style','Radio','String',...
    'Grid Indexing','Position',[5,2,90,15],'Value',0,...
    'BackgroundColor','w');



PHYSVIEW=uibuttongroup(hf,'Units','Pixels','Position',...
    [724,0,225,20],'SelectionChangeFCN',@PHYSAXIS...
    ,'BackgroundColor','w');
NormView=uicontrol(PHYSVIEW,'Units','Pixels','Style','Radio','String',...
    'Rendered View','Position',[105,2,110,15],'Value',0,...
    'BackgroundColor','w');
FullView=uicontrol(PHYSVIEW,'Units','Pixels','Style','Radio','String',...
    'Full View','Position',[5,2,90,15],'Value',1,'BackgroundColor','w');

RUN=uicontrol(hf,'Units','Pixels','Position',[462,0,100,40],'String',...
    'RUN','FontSize',24,'CallBack',@RUNCallBack);

VARPAR=uipanel(hf,'Title','Variable Parameters','TitlePosition',...
    'CenterTop','FontSize',12,'Units','Pixels','Position',...
    [770/3+20,488,770/3,100],'BackgroundColor',[.8,.8,.8]);
NOVAR=uicontrol(VARPAR,'Style','togglebutton','Units','pixels',...
    'Position',[140,30,50,15],'String','NONE','BackgroundColor',...
    [.2,.2,.2],'Max',true,'Min',false,'Value',true,'Enable','on',...
    'CallBack',@NOVARCallBack,'ForegroundColor','w');
VARMAP=uicontrol(VARPAR,'Units','Pixels','Position',[140,2,50,20],...
    'String','EDIT','BackgroundColor',[.8,.8,.8],'CallBack',@VARMAPCallBack);

PARAM=uibuttongroup(VARPAR,'Units','Pixels','Position',[2,2,192,80],...
    'BorderType','none','BackgroundColor',[.8,.8,.8]);
VAR_A=uicontrol('Parent',PARAM,'Style','checkbox','Units','pixels',...
    'Position',[0,62,130,15],'String','Circulation Parameter',...
    'BackgroundColor',[.8,.8,.8],'Enable','off');
VAR_RE=uicontrol('Parent',PARAM,'Style','checkbox','Units','pixels',...
    'Position',[0,42,130,15],'String','Reynolds Number',...
    'BackgroundColor',[.8,.8,.8],'Enable','off');
VAR_Y=uicontrol('Parent',PARAM,'Style','checkbox','Units','pixels',...
    'Position',[0,22,130,15],'String','Y-Displacement','BackgroundColor',...
    [.8,.8,.8],'Enable','off');
VAR_J=uicontrol('Parent',PARAM,'Style','checkbox','Units','pixels',...
    'Position',[0,2,130,15],'String','Jet Parameters','BackgroundColor',...
    [.8,.8,.8],'Enable','off');

MISC=uipanel(hf,'Title','Miscellaneous','TitlePosition','CenterTop','FontSize',12,'Units','Pixels','Position',[30+2*770/3,488,770/3,100],'BackgroundColor',[.8,.8,.8]);
UVELFFT=uicontrol(MISC,'Style','Checkbox','String','UVEL FFT','Units','Pixels','Position',[183,60,70,20],'BackgroundColor',[.8,.8,.8],'CallBack',@UVELFFTCallBack);
UVEL_YData=[2:20]';
UVEL_YData=num2str(UVEL_YData);
UVEL_Y=uicontrol(MISC,'Style','listbox','String',UVEL_YData,'Value',1,'Min',18,'Max',19,'Units','Pixels','Position',[183,0,72,60],'BackgroundColor','w','Enable','off','CallBack',@UVELFFTCallBack);
AVG_DATA=uicontrol(MISC,'Style','Checkbox','String','Output Average Data','Units','Pixels','Position',[2,65,130,15],'BackgroundColor',[.8 .8 .8],'CallBack',@AVG_DATACallBack);
AVG_VAR=uicontrol(MISC,'Style','popupmenu','String',{'Select Variable','All','Vorticity','Stream-Function','U-Velocity','V-Velocity'},'Position',[2,50,100,15],'Enable','Off','CallBack',@AVG_VARCallBack);
AVG_YRG=uicontrol(MISC,'Style','popupmenu','String',{'Select Y-Range','Full','2:IBL','Custom'},'Position',[2,30,100,15],'Enable','Off','CallBack',@AVG_YRGCallBack);
AVG_XRG=uicontrol(MISC,'Style','popupmenu','String',{'Select X-Range','Full','Custom'},'Position',[2,10,100,15],'Enable','Off','CallBack',@AVG_XRGCallBack);

% MEMSIZE = uipanel(hf,'Title','Memory Statistics','TitlePosition','CenterTop','Fontsize',12,'Units','Pixels','Position',[400,100,220,120],'BackgroundColor','w');
% memsize = ((xwidth+1)+((yheight+1)*(xwidth+1))*7+((yheight+1)*(xwidth+1))*4)*8; %double precision
% gpuInfo = gpuDevice;
% MEMSIZE_GPU = uicontrol(MEMSIZE,'Style','Text','String',{[num2str(gpuInfo.TotalMemory/(1024^2)) ' Total MB on GPU'];[num2str(gpuInfo.FreeMemory/(1024^2)) ' Free MB on GPU'];[num2str(memsize/1024/1024) ' MB Required (Est.)']},'Units','Pixels','Position',[0,20,210,60],'FontSize',12,'BackgroundColor','w');

%% Figure 2
hf2=figure('MenuBar','none','Name','Variable Parameters',...
    'Position',[0,0,800,600],'Units','pixels','NumberTitle',...
    'off','Resize','off','Color','w','Visible','off');
data=[get(A,'UserData') get(Re,'UserData') 0 {'Tolerance'},...
                                                      get(Tol,'UserData')];
TBL=uitable(hf2,'Position',[0,50,802,552],'ColumnName',...
    {'Circulation Parameter','Modified Reynolds Number','Y-Position',...
    'Increment Criteria','Incremement Value'},'ColumnFormat',...
    {'numeric','numeric','numeric',{'Tolerance','Iteration'},'numeric'},...
    'ColumnEditable',[true true true true true],'Data',data);
ADDROW=uicontrol(hf2,'Position',[400,0,100,50],'String','Add Row',...
    'CallBack',@ADDROWCallBack);
%% Computation Methods
PFI_VER=uibuttongroup(hf,'Units','Pixels','Position',[10,433,125,55],...
    'Title','Method','FontSize',12,'TitlePosition','CenterTop',...
    'BorderType','none','BackgroundColor',[.8,.8,.8]);


g = gpuDevice;
if g.DeviceSupported
    PFI_GPU=uicontrol(PFI_VER,'Units','Pixels','Style','Radio','String',...
        'CUDA-GPU Compute','Position',[1,20,125,20],'BackgroundColor',...
        [.8,.8,.8]);
end

PFI_NRM=uicontrol(PFI_VER,'Units','Pixels','Style','Radio','String',...
    'Serial CPU Compute','Position',[1,0,125,20],'BackgroundColor',...
    [.8,.8,.8]);
if g.DeviceSupported
    set(PFI_GPU,'Value',1)
    set(PFI_NRM,'Value',0)
else
    set(PFI_NRM,'Value',0)
end

%% Call Back Functions
    function AVG_YRGCallBack
        if get(AVG_YRG,'Value')==1;  % Nothing Selected Yet...
            set(AVG_YRG,'UserData',0)
        elseif get(AVG_YRG,'Value')==2;  % Full Range
            yheight1 = get(My,'UserData')+1;
            set(AVG_YRG,'UserData',yheight1)
        elseif get(AVG_YRG,'Value')==3;  % IBL is Range
            yheight1 = get(IBL,'UserData');
            set(AVG_YRG,'UserData',yheight1)
        elseif get(AVG_YRG,'Value')==4;  % Custom Range
            yheight1 = inputdlg({'Enter Min Y-Value';'Enter Max Y-Value'},'Enter Custom Range',2);
        end
    end
    function AVG_DATACallBack(varargin)
        xwidth = get(Nx,'UserData')+1;
        yheight = get(My,'UserData')+1;
        if get(AVG_DATA,'Value')
            if get(AVG_VAR,'Value')==1  % Nothing Selected Yet...
                NumVar = 0;               
            elseif get(AVG_VAR,'Value')==2 % All Variables
                NumVar = 4;
            else                            % Only One Variable
                NumVar = 1;
            end
            
%             memsize = (xwidth+(yheight*xwidth)*7+(yheight*xwidth)*4 + (yheight1*xwidth1)*4)*8; %double precision
%             gpuInfo = gpuDevice;
%             set(MEMSIZE_GPU,'String',{[num2str(gpuInfo.TotalMemory/(1024^2)) ' Total MB on GPU'];[num2str(gpuInfo.FreeMemory/(1024^2)) ' Free MB on GPU'];[num2str(memsize/1024/1024) ' MB Required (Est.)']});
%         
%         else
%             memsize = (xwidth+(yheight*xwidth)*7+(yheight*xwidth)*4)*8; %double precision
%             gpuInfo = gpuDevice;
%             set(MEMSIZE_GPU,'String',{[num2str(gpuInfo.TotalMemory/(1024^2)) ' Total MB on GPU'];[num2str(gpuInfo.FreeMemory/(1024^2)) ' Free MB on GPU'];[num2str(memsize/1024/1024) ' MB Required (Est.)']}); 
        end
        
    end

    function UVELFFTCallBack(varargin)
        if get(UVELFFT,'Value')
            set(UVEL_Y,'Enable','on')
            IBLmax=get(IBL,'UserData');
            UVEL_YData=[2:IBLmax]';
            UVEL_YData=num2str(UVEL_YData);
            set(UVEL_Y,'String',UVEL_YData)
        else
            set(UVEL_Y,'Enable','off')
        end
        
    end

    function JETFRQCallBack(varargin)
        selJet=get(JETLIST,'Value');
        selJetFRQ=get(JETFRQ,'UserData');
        tempJet=inputdlg('Enter Frequency',char({['Jet ' num2str(selJet)]}));
        selJetFRQ(selJet)=str2double(tempJet{1});
        set(JETFRQ,'UserData',selJetFRQ,'String',{['Frequency= ' num2str(selJetFRQ(selJet))]})
    end

    function JETAMPCallBack(varargin)
        selJet=get(JETLIST,'Value');
        selJetAMP=get(JETAMP,'UserData');
        tempJet=inputdlg('Enter Amplitude',char({['Jet ' num2str(selJet)]}));
        selJetAMP(selJet)=str2double(tempJet{1});
        set(JETAMP,'UserData',selJetAMP,'String',{['Amplitude= ' num2str(selJetAMP(selJet))]})
    end

    function JETENDCallBack(varargin)
        selJet=get(JETLIST,'Value');
        selJetEND=get(JETEND,'UserData');
        tempJet=inputdlg('Enter End Location',char({['Jet ' num2str(selJet)]}));
        selJetEND(selJet)=str2double(tempJet{1});
        set(JETEND,'UserData',selJetEND,'String',{['End= ' num2str(selJetEND(selJet))]})
        GridCompute
    end

    function JETSTARTCallBack(varargin)
        selJet=get(JETLIST,'Value');
        selJetSTART=get(JETSTART,'UserData');
        tempJet=inputdlg('Enter Start Location',char({['Jet ' num2str(selJet)]}));
        selJetSTART(selJet)=str2double(tempJet{1});
        set(JETSTART,'UserData',selJetSTART,'String',{['Start= ' num2str(selJetSTART(selJet))]})
        GridCompute
    end

    function JETLISTCallBack(varargin)
        selJet=get(JETLIST,'Value');
        JETSTARTdata=get(JETSTART,'UserData');
        JETENDdata=get(JETEND,'UserData');
        JETAMPdata=get(JETAMP,'UserData');
        JETFRQdata=get(JETFRQ,'UserData');
        set(JETSTART,'String',{['Start= ' num2str(JETSTARTdata(selJet))]})
        set(JETEND,'String',{['End= ' num2str(JETENDdata(selJet))]})
        set(JETAMP,'String',{['Amplitude= ' num2str(JETAMPdata(selJet))]})
        set(JETFRQ,'String',{['Frequency= ' num2str(JETFRQdata(selJet))]})
    end

    function JETADDCallBack(varargin)
        currJets=get(JETLIST,'String');
        NumJets=length(currJets);
        currJets(NumJets+1)={['JET ' num2str(NumJets+1)]};
        set(JETLIST,'String',currJets)
        set(JETLIST,'Value',NumJets+1)
        
        JETSTARTdata=get(JETSTART,'UserData');
        JETENDdata=get(JETEND,'UserData');
        JETAMPdata=get(JETAMP,'UserData');
        JETFRQdata=get(JETFRQ,'UserData');
        
        JETSTARTdata(NumJets+1)=JETSTARTdata(NumJets);
        JETENDdata(NumJets+1)=JETENDdata(NumJets);
        JETAMPdata(NumJets+1)=JETAMPdata(NumJets);
        JETFRQdata(NumJets+1)=JETFRQdata(NumJets);
        
        set(JETSTART,'String',{['Start= ' num2str(JETSTARTdata(NumJets+1))]},'UserData',JETSTARTdata)
        set(JETEND,'String',{['End= ' num2str(JETENDdata(NumJets+1))]},'UserData',JETENDdata)
        set(JETAMP,'String',{['Amplitude= ' num2str(JETAMPdata(NumJets+1))]},'UserData',JETAMPdata)
        set(JETFRQ,'String',{['Frequency= ' num2str(JETFRQdata(NumJets+1))]},'UserData',JETFRQdata)
    end

    function JETTESTCallBack(varargin)
        if get(JETTEST,'Value')
            set(JETLIST,'Enable','on')
            set(JETADD,'Enable','on')
            set(JETSTART,'Enable','on')
            set(JETEND,'Enable','on')
            set(JETAMP,'Enable','on')
            set(JETFRQ,'Enable','on')
            if get(NOVAR,'Value')
                set(VAR_J,'Value',0,'Enable','off')
            else
                set(VAR_J,'Enable','on')
            end
        else
            set(JETLIST,'String',{'JET 1'},'Enable','off')
            set(JETADD,'Enable','off')
            set(JETSTART,'Enable','off')
            set(JETEND,'Enable','off')
            set(JETAMP,'Enable','off')
            set(JETFRQ,'Enable','off') 
            set(VAR_J,'Value',0,'Enable','off')
        end
    end

    function NOVARCallBack(varargin)
        if get(NOVAR,'Value')
            set(VAR_A,'Value',0,'Enable','off')
            set(VAR_RE,'Value',0,'Enable','off')
            set(VAR_Y,'Value',0,'Enable','off')
            set(VAR_J,'Value',0,'Enable','off')
            set(NOVAR,'BackgroundColor',[.2,.2,.2],'ForegroundColor','w')
        else
            set(VAR_A,'Enable','on')
            set(VAR_RE,'Enable','on')
            set(VAR_Y,'Enable','on')
            if get(JETTEST,'Value')
                set(VAR_J,'Enable','on')
            else
                set(VAR_J,'Value',0,'Enable','off')
            end
            set(NOVAR,'BackgroundColor',[.8,.8,.8],'ForegroundColor','k')
        end
    end

    function VARMAPCallBack(varargin)
        VAL=logical([get(VAR_A,'Value') get(VAR_RE,'Value')...
                                                  get(VAR_Y,'Value') 1 1]);
        set(TBL,'ColumnEditable',VAL)
        set(hf2,'Visible','on')
        movegui(hf2,'center')        
    end

    function ADDROWCallBack(varargin)
        data=get(TBL,'Data');
        data{end+1,1}=data{end,1};
        data{end,2}=data{end-1,2};
        data{end,3}=data{end-1,3};
        data{end,4}=data{end-1,4};
        data{end,5}=data{end-1,5};
        set(TBL,'Data',data);
    end


    function COMPAXIS(varargin)
        XWidth=get(Nx,'UserData')+1;
        YHeight=get(My,'UserData')+1;
        xmin=get(Xmin,'UserData');
        xmax=get(Xmax,'UserData');
        ymin=get(Ymin,'UserData');
        ymax=get(Ymax,'UserData');
                
        if get(IndxGrid,'Value')
            comp_x=1:XWidth;
            comp_y=1:YHeight;
            [comp_x,comp_y]=meshgrid(comp_x,comp_y);
            set(CGrid,'XData',comp_x,'YData',comp_y)
            set(COMPGRID,'YDir','reverse')
        elseif get(CalcGrid,'Value')
            clear comp_x comp_y
            dx=abs(xmax-xmin)/(XWidth-1);
            comp_x(1)=xmin;
            comp_x(XWidth)=xmax;
            for ii=2:XWidth-1
                comp_x(ii)=xmin+dx*(ii-1);  %#ok<AGROW>
            end

            dy=abs(ymax-ymin)/(YHeight-1);
            comp_y(1)=ymin;
            comp_y(YHeight)=ymax;
            for ii=2:YHeight-1
                comp_y(ii)=ymin+dy*(ii-1);  %#ok<AGROW>
            end
            [comp_x,comp_y]=meshgrid(comp_x,comp_y);
            set(CGrid,'XData',comp_x,'YData',comp_y)
            set(COMPGRID,'YDir','normal')
        end
    end

    function PHYSAXIS(varargin)
        if get(NormView,'Value')
            set(hf,'CurrentAxes',PHYSGRID)
            axis equal
            axis([-2 16 -2 8])
        elseif get(FullView,'Value')
            set(hf,'CurrentAxes',PHYSGRID)
            axis auto
            axis equal
        end
    end


    function [PGrid,CGrid]=GridCompute(varargin)
        nx=get(Nx,'UserData');
        my=get(My,'UserData');
        xmin=get(Xmin,'UserData');
        xmax=get(Xmax,'UserData');
        ymin=get(Ymin,'UserData');
        ymax=get(Ymax,'UserData');
        xskew=get(Xskew,'UserData');
        yskew=get(Yskew,'UserData');
        
        currJets=get(JETLIST,'String');
        NumJets=length(currJets);
        JetStartLocs=get(JETSTART,'UserData');
        JetEndLocs=get(JETEND,'UserData');
        
        XWidth=nx+1;
        YHeight=my+1;
        dx=abs(xmax-xmin)/(nx);
        comp_x(1)=xmin;
        comp_x(nx+1)=xmax;
        for ii=2:nx
            comp_x(ii)=xmin+dx*(ii-1); %#ok<AGROW>
        end

        dy=abs(ymax-ymin)/(my);
        comp_y(1)=ymin;
        comp_y(my+1)=ymax;
        for ii=2:my
            comp_y(ii)=ymin+dy*(ii-1); %#ok<AGROW>
        end
        [comp_x,comp_y]=meshgrid(comp_x,comp_y);
        [COLX,COLY]=meshgrid(1:XWidth,1:YHeight);
        COL(:,:,1) = COLX/(abs(max(max(COLX)))-min(min(COLX))); 
        COL(:,:,3) = COLY/(abs(max(max(COLY))-min(min(COLY))));
        COL(:,:,2) = .5*ones(size(COLX));
        COL(1:get(IBL,'UserData'),:,:)=0;
        
        
       
        
        CGrid=mesh(COMPGRID,comp_x,comp_y,zeros(size(comp_x)),COL);
        view(COMPGRID,[0,90])
        set(COMPGRID,'XGrid','off','YGrid','off','ZGrid','off','Box','on')
        title(COMPGRID,'Computational Map','FontSize',18)
        set(hf,'CurrentAxes',COMPGRID)
        axis equal


        x=zeros(1,XWidth);
        y=zeros(1,YHeight);
        for ii=1:XWidth
            if (ii<=abs(xmin/(xmax-xmin))*(XWidth-1))
                x(ii)=(-1/2)*abs(xmax-xmin)*abs((ii-((XWidth-1)/2)-1)/...
                                                     ((XWidth-1)/2))^xskew;
            else
                x(ii)=(1/2)*abs(xmax-xmin)*abs((ii-((XWidth-1)/2)-1)/...
                                                    ((XWidth-1)/2))^xskew;
            end
        end
        for ii=1:YHeight
            y(ii)=ymin+abs(ymax-ymin)*((ii-1)/(YHeight-1))^yskew;
        end
        [mu,eta]=meshgrid(x,y);
        x=1/2*(mu.^2-eta.^2);
        y=mu.*eta;
        
        cla(PHYSGRID)
        PGrid=mesh(PHYSGRID,x,y,zeros(size(x)),COL);
        set(PHYSGRID,'nextplot','add')
        if get(JETTEST,'Value')
            for ii=1:NumJets
               plot(PHYSGRID,x(1,JetStartLocs(ii):JetEndLocs(ii)),y(1,JetStartLocs(ii):JetEndLocs(ii)),'b','LineWidth',2);
            end
        end
        view(PHYSGRID,[0,90])
        set(PHYSGRID,'XGrid','off','YGrid','off','ZGrid','off','Box','on')
        title(PHYSGRID,'Physical Map','FontSize',18)
        if exist('NormView','var')
            PHYSAXIS
        end
    end

    function XminCallBack(varargin)
            DefAns=num2str(get(Xmin,'UserData'));
            answer = inputdlg('Enter Value for Xmin','Xmin',1,{DefAns});
            set(Xmin,'UserData',str2double(answer),'String',...
                                                strcat('Xmin = ',(answer)))
            Dt=delta_t(get(Nx,'UserData'),get(My,'UserData'),...
                get(Xmin,'UserData'),get(Xmax,'UserData'),...
                get(Ymin,'UserData'),get(Ymax,'UserData'),...
                get(Re,'UserData'));
            set(dt,'UserData',Dt,'String',strcat('dt = ',num2str(Dt)))
            [PGrid,CGrid]=GridCompute;
            COMPAXIS;
            PHYSAXIS;
    end
    function XmaxCallBack(varargin)
            DefAns=num2str(get(Xmax,'UserData'));
            answer = inputdlg('Enter Value for Xmax','Xmax',1,{DefAns});
            set(Xmax,'UserData',str2double(answer),'String',...
                                                strcat('Xmax = ',(answer)))
            Dt=delta_t(get(Nx,'UserData'),get(My,'UserData'),...
                get(Xmin,'UserData'),get(Xmax,'UserData'),...
                get(Ymin,'UserData'),get(Ymax,'UserData'),...
                get(Re,'UserData'));
            set(dt,'UserData',Dt,'String',strcat('dt = ',num2str(Dt)))                                            
            [PGrid,CGrid]=GridCompute;
            COMPAXIS;
            PHYSAXIS;
    end
    function YminCallBack(varargin)
            DefAns=num2str(get(Ymin,'UserData'));
            answer = inputdlg('Enter Value for Ymin','Ymin',1,{DefAns});
            set(Ymin,'UserData',str2double(answer),'String',...
                                                strcat('Ymin = ',(answer)))
            Dt=delta_t(get(Nx,'UserData'),get(My,'UserData'),...
                get(Xmin,'UserData'),get(Xmax,'UserData'),...
                get(Ymin,'UserData'),get(Ymax,'UserData'),...
                get(Re,'UserData'));
            set(dt,'UserData',Dt,'String',strcat('dt = ',num2str(Dt)))                                            
            [PGrid,CGrid]=GridCompute;
            COMPAXIS;
            PHYSAXIS;
    end
    function YmaxCallBack(varargin)
            DefAns=num2str(get(Ymax,'UserData'));
            answer = inputdlg('Enter Value for Ymax','Ymax',1,{DefAns});
            set(Ymax,'UserData',str2double(answer),'String',...
                                                strcat('Ymax = ',(answer)))
            Dt=delta_t(get(Nx,'UserData'),get(My,'UserData'),...
                get(Xmin,'UserData'),get(Xmax,'UserData'),...
                get(Ymin,'UserData'),get(Ymax,'UserData'),...
                get(Re,'UserData'));
            set(dt,'UserData',Dt,'String',strcat('dt = ',num2str(Dt)))                                            
            [PGrid,CGrid]=GridCompute;
            COMPAXIS;
            PHYSAXIS;
    end
    function NxCallBack(varargin)
            DefAns=num2str(get(Nx,'UserData'));
            answer = inputdlg('Enter Value for Nx','Nx',1,{DefAns});
            set(Nx,'UserData',str2double(answer),'String',...
                                                  strcat('Nx = ',(answer)))
            Dt=delta_t(get(Nx,'UserData'),get(My,'UserData'),...
                get(Xmin,'UserData'),get(Xmax,'UserData'),...
                get(Ymin,'UserData'),get(Ymax,'UserData'),...
                get(Re,'UserData'));
            set(dt,'UserData',Dt,'String',strcat('dt = ',num2str(Dt)))
            
            xwidth = get(Nx,'UserData')+1;
            yheight = get(My,'UserData')+1;
            
%             if get(AVG_DATA,'Value')
%                 memsize = (xwidth+(yheight*xwidth)*7+(yheight*xwidth)*4*2)*8; %double precision
%                 gpuInfo = gpuDevice;
%                 set(MEMSIZE_GPU,'String',{[num2str(gpuInfo.TotalMemory/(1024^2)) ' Total MB on GPU'];[num2str(gpuInfo.FreeMemory/(1024^2)) ' Free MB on GPU'];[num2str(memsize/1024/1024) ' MB Required (Est.)']});
%         
%             else
%                 memsize = (xwidth+(yheight*xwidth)*7+(yheight*xwidth)*4)*8; %double precision
%                 gpuInfo = gpuDevice;
%                 set(MEMSIZE_GPU,'String',{[num2str(gpuInfo.TotalMemory/(1024^2)) ' Total MB on GPU'];[num2str(gpuInfo.FreeMemory/(1024^2)) ' Free MB on GPU'];[num2str(memsize/1024/1024) ' MB Required (Est.)']}); 
%             end
        
            [PGrid,CGrid]=GridCompute;
            COMPAXIS;
            PHYSAXIS;
    end
    function MyCallBack(varargin)
            DefAns=num2str(get(My,'UserData'));
            answer = inputdlg('Enter Value for My','My',1,{DefAns});
            set(My,'UserData',str2double(answer),'String',...
                                                  strcat('My = ',(answer)))
            set(IBL,'UserData',ceil(get(My,'UserData')/20),'String',...
                strcat('IBL = ',num2str(ceil(get(My,'UserData')/20))))
            
            Dt=delta_t(get(Nx,'UserData'),get(My,'UserData'),...
                get(Xmin,'UserData'),get(Xmax,'UserData'),...
                get(Ymin,'UserData'),get(Ymax,'UserData'),...
                get(Re,'UserData'));
            set(dt,'UserData',Dt,'String',strcat('dt = ',num2str(Dt)))   
            
            xwidth = get(Nx,'UserData')+1;
            yheight = get(My,'UserData')+1;
            
%             if get(AVG_DATA,'Value')
%                 memsize = (xwidth+(yheight*xwidth)*7+(yheight*xwidth)*4*2)*8; %double precision
%                 gpuInfo = gpuDevice;
%                 set(MEMSIZE_GPU,'String',{[num2str(gpuInfo.TotalMemory/(1024^2)) ' Total MB on GPU'];[num2str(gpuInfo.FreeMemory/(1024^2)) ' Free MB on GPU'];[num2str(memsize/1024/1024) ' MB Required (Est.)']});
%         
%             else
%                 memsize = (xwidth+(yheight*xwidth)*7+(yheight*xwidth)*4)*8; %double precision
%                 gpuInfo = gpuDevice;
%                 set(MEMSIZE_GPU,'String',{[num2str(gpuInfo.TotalMemory/(1024^2)) ' Total MB on GPU'];[num2str(gpuInfo.FreeMemory/(1024^2)) ' Free MB on GPU'];[num2str(memsize/1024/1024) ' MB Required (Est.)']}); 
%             end
            
            [PGrid,CGrid]=GridCompute;
            COMPAXIS;
            PHYSAXIS;
    end
    function XskewCallBack(varargin)
            DefAns=num2str(get(Xskew,'UserData'));
            answer = inputdlg('Enter Value for Xskew','Xskew',1,{DefAns});
            set(Xskew,'UserData',str2double(answer),'String',...
                                               strcat('Xskew = ',(answer)))
            [PGrid,CGrid]=GridCompute;
            COMPAXIS;
            PHYSAXIS;
    end
    function YskewCallBack(varargin)
            DefAns=num2str(get(Yskew,'UserData'));
            answer = inputdlg('Enter Value for Yskew','Yskew',1,{DefAns});
            set(Yskew,'UserData',str2double(answer),'String',...
                                               strcat('Yskew = ',(answer)))
            [PGrid,CGrid]=GridCompute;
            COMPAXIS;
            PHYSAXIS;
    end
    function ACallBack(varargin)
            DefAns=num2str(get(A,'UserData'));
            answer = inputdlg('Enter Value for A','A',1,{DefAns});
            set(A,'UserData',str2double(answer),'String',...
                                                   strcat('A = ',(answer)))
    end
    function ReCallBack(varargin)
            DefAns=num2str(get(Re,'UserData'));
            answer = inputdlg('Enter Value for Re','Re',1,{DefAns});
            set(Re,'UserData',str2double(answer),'String',...
                                                  strcat('Re = ',(answer)))
            Dt=delta_t(get(Nx,'UserData'),get(My,'UserData'),...
                get(Xmin,'UserData'),get(Xmax,'UserData'),...
                get(Ymin,'UserData'),get(Ymax,'UserData'),...
                get(Re,'UserData'));
            set(dt,'UserData',Dt,'String',strcat('dt = ',num2str(Dt)))                                              
    end
    function dtCallBack(varargin)
            DefAns=num2str(get(dt,'UserData'));
            answer = inputdlg('Enter Value for Dt','Dt',1,{DefAns});
            set(dt,'UserData',str2double(answer),'String',...
                                                  strcat('Dt = ',(answer)))
    end
    function TolCallBack(varargin)
            DefAns=num2str(get(Tol,'UserData'));
            answer = inputdlg('Enter Value for Tolerance','Tolerance',1,...
                                                                {DefAns});
            set(Tol,'UserData',str2double(answer),'String',...
                                                strcat('TOL = ',(answer)))
    end
    function IBLCallBack(varargin)
            DefAns=num2str(get(IBL,'UserData'));
            answer = inputdlg('Enter Value for IBL','IBL',1,{DefAns});
            set(IBL,'UserData',str2double(answer),'String',...
                                                strcat('IBL = ',(answer)))
    end
    function reportCallBack(varargin)
            DefAns=num2str(get(report,'UserData'));
            answer=inputdlg('Enter Value for Report','Report',1,{DefAns});
            set(report,'UserData',str2double(answer),'String',...
                                    strcat('Report : ',(answer),' ITER'))
    end
    function PSAVECallBack(varargin)
            DefAns=num2str(get(PSAVE,'UserData'));
            answer = inputdlg('Enter Value for FileWrite','FileWrite',1,...
                                                                 {DefAns});
            set(PSAVE,'UserData',str2double(answer),'String',...
                                strcat('FileWrite : ',(answer),' ITER'))
    end
    function OtCallBack(varargin)
            DefAns=num2str(get(Ot,'UserData'));
            answer = inputdlg('Enter Value for Iterations','Iterations',...
                                                               1,{DefAns});
            set(Ot,'UserData',str2double(answer),'String',...
                                        strcat('Iterations : ',(answer)))
    end
    function filenameCallBack(varargin)
            DefAns=get(filename,'UserData');
            answer = inputdlg('Enter FileName','FileName',1,DefAns);
            set(filename,'UserData',answer,'String',{strcat('FileName : '...
                                                              ,answer{1})})
    end


    function RUNCallBack(varargin)
        if get(PFI_NRM,'Value')
                PFI_CPU(SimType,get(Xmin,'UserData'),get(Xmax,'UserData'),...
                get(Ymin,'UserData'),get(Ymax,'UserData'),...
                get(Nx,'UserData'),get(My,'UserData'),...
                get(Re,'UserData'),get(A,'UserData'),get(Ot,'UserData'),...
                get(report,'UserData'),get(dt,'UserData'),...
                get(Tol,'UserData'),get(IBL,'UserData'),...
                cell2mat(get(filename,'UserData')),get(TBL,'Data'),...
                get(VAR_RE,'Value'),get(VAR_A,'Value'),...
                get(VAR_Y,'Value'),get(VAR_J,'Value'),...
                get(JETTEST,'Value'),get(JETLIST,'String'),get(JETSTART,'UserData'),...
                get(JETEND,'UserData'),get(JETFRQ,'UserData'),...
                get(JETAMP,'UserData'),get(PSAVE,'UserData'),...
                get(UVELFFT,'Value'),get(UVEL_Y,'Value'),get(AVG_DATA,'Value'),get(hf,'UserData'));
        elseif get(PFI_GPU,'Value')
                PFI_CUDA(SimType,get(Xmin,'UserData'),get(Xmax,'UserData'),...
                get(Ymin,'UserData'),get(Ymax,'UserData'),...
                get(Nx,'UserData'),get(My,'UserData'),...
                get(Re,'UserData'),get(A,'UserData'),get(Ot,'UserData'),...
                get(report,'UserData'),get(dt,'UserData'),...
                get(Tol,'UserData'),get(IBL,'UserData'),...
                cell2mat(get(filename,'UserData')),get(TBL,'Data'),...
                get(VAR_RE,'Value'),get(VAR_A,'Value'),...
                get(VAR_Y,'Value'),get(VAR_J,'Value'),...
                get(JETTEST,'Value'),get(JETLIST,'String'),get(JETSTART,'UserData'),...
                get(JETEND,'UserData'),get(JETFRQ,'UserData'),...
                get(JETAMP,'UserData'),get(PSAVE,'UserData'),...
                get(UVELFFT,'Value'),get(UVEL_Y,'Value'),get(AVG_DATA,'Value'),get(hf,'UserData'));
                
%             end
         end
    end


    
end