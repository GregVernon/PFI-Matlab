function PFI_CPU(SimType,Xmin,Xmax,Ymin,Ymax,Nx,My,Re,A,Ot,report,dt,Tol,IBL,...
                        filename,VARPAR,dRe,dA,dY,dJ,JETTEST,JETLIST,JETSTART,JETEND,JETFRQ,JETAMP,psave,UVELFFT,UVEL_Y,AVG_DATA,oldfilename)    
    %% Load Data
    k=0;
    kerr=1;
    tick=0;	
    ct=0;
    Nx=Nx-1;
    My=My-1;
    
    LL=length(filename);
    filename(LL+1)='P';
    NumOut=floor(Ot/psave);

    %% Generate Grid
    dx=abs(Xmax-Xmin)/(Nx+1);
    x(1)=Xmin;
    x(Nx+2)=Xmax;
    for ii=2:Nx+1
        x(ii)=Xmin+dx*(ii-1); %#ok<AGROW>
    end
    dx2=2*dx;
    dxx=dx*dx;

    dy=abs(Ymax-Ymin)/(My+1);
    y(1)=Ymin;
    y(My+2)=Ymax;
    for ii=2:My+1
        y(ii)=Ymin+dy*(ii-1); %#ok<AGROW>
    end
    dy2=2*dy;
    dyy=dy*dy;

    xx=repmat(x,My+2,1);
    yy=repmat(y',1,Nx+2);
    ii=2:Nx+1;
    jj=2:My+1;
    d(1:My,1:Nx)=sqrt(xx(jj,ii).^2+yy(jj,ii).^2);
    d2(1:My,1:Nx)=xx(2:My+1,2:Nx+1).^2+yy(2:My+1,2:Nx+1).^2;
    d3(jj,ii)=sqrt(xx(jj,ii).^2+yy(jj,ii).^2);
    d6=x.^2+y(1)^2;
    dip1(1:My,1:Nx)=sqrt(xx(jj,ii+1).^2+yy(jj,ii).^2);
    dim1(1:My,1:Nx)=sqrt(xx(jj,ii-1).^2+yy(jj,ii).^2);
    djp1(1:My,1:Nx)=sqrt(xx(jj,ii).^2+yy(jj+1,ii).^2);
    djm1(1:My,1:Nx)=sqrt(xx(jj,ii).^2+yy(jj-1,ii).^2);

    Kappa2=(dx/dy)^2;
    KappaA = 1.0/(2.0*(1+Kappa2));
    Rc = Re*dx;

    Cx = dt/dx;
    Cx2 = .5*Cx;
    Cy = dt/dy;
    Cy2 = .5*Cy;

    if Cx>Cy
        C=Cx;
    else
        C=Cy;
    end

    alphaX = dt/(dxx*Re);
    alphaY = dt/(dyy*Re);
    alpha = 2*alphaX + 2*alphaY;

    COURSTR=[{['Courant Number C =',num2str(C)],...
        ' Must be less than 1'},{['The Cell Reynolds Number Rc =',...
        num2str(Rc)],['MUST BE LESS THAN 4/C =', num2str(4/C)]},...
       {'Grid Spacing dx, dy, dt:',mat2str([dx dy dt])},{'CONTINUE?'}];
    COUR=questdlg(COURSTR);
    if strcmpi(COUR,'Yes')
        close all
        
        Progress=figure('Name',['Re = ',num2str(Re),' A = ',num2str(A)],'MenuBar','none','Position',...
                                                            [0,0,500,600]);
        movegui(Progress,'center')
        str1=sprintf(['Iteration','\t','kPsi','\t','Omega-Tol','\t','Psi-Tol','\t',...
                                                            'RunTime(s)']);
        STR1=[{'Report'};{str1}];
        Col1=uicontrol(Progress,'Style','Edit','Position',[0,40,300,560],...
                                        'String',STR1,'Max',1000,'Min',1);
        STR2='Incremental Files';
        Col2=uicontrol(Progress,'Style','Edit','Position',...
                        [300,40,200,560],'String',STR2,'Max',1000,'Min',1);
        
        pause(0.1)
        SNAP=uicontrol(Progress,'Style','ToggleButton','Position',[400,5,100,30],'Max',true,'Min',false,'Value',false,'String','SNAPSHOT');
    else exit
    end
    %% Initialize
    if SimType % Continue from old file
        DATA=load(char(strcat(oldfilename(1),oldfilename(2))));
        Omega=DATA(1:My+2,1:end);
        Psi=DATA(1*(My+2)+1:2*(My+2),1:end);
        u=DATA(2*(My+2)+1:3*(My+2),1:end);
        v=DATA(3*(My+2)+1:4*(My+2),1:end);
    else
        v=zeros(My+2,Nx+2);
        u=zeros(My+2,Nx+2);
        Psi=zeros(My+2,Nx+2);
        Omega=zeros(My+2,Nx+2);
        for jj=1:Nx+2
            for ii=1:My+2
                v(ii,jj)=-(y(ii)-1)/sqrt(x(jj)^2+y(ii)^2);
            end
        end
        for jj=1:Nx+2
            for ii=2:My+2
                u(ii,jj)=(x(jj)+A)/sqrt(x(jj)^2+y(ii)^2);
            end
        end
        for jj=1:Nx+2
            for ii=1:My+2
                Psi(ii,jj)=(x(jj)+A)*(y(ii)-1);
            end
        end
        for ii=1:Nx+2
                Omega(1,ii) = (7.0*Psi(1,ii)-8.0*Psi(2,ii)+Psi(3,ii))...
                                                           /(2.0*dyy)/d6(1,ii);
        end
    end
    %% Compute
    
    avgdata = 0;
    
    time=tic;
    numJets=length(JETLIST);
    while k<Ot
        drawnow
        if get(SNAP,'Value')
            set(SNAP,'Value',false)
            SNAP_EXE(STR2,length(y),length(x),Xmin,Xmax,Ymin,Ymax,gather(u),gather(v));
        end
        k=k+1;
        t=k*dt;
        
        if AVG_DATA
            avgdata = (avgdata*(k-1)+[Omega;Psi;u;v])/k;
        end
        
        if (dA||dRe||dY||dJ)
            if(strcmpi(VARPAR{CYC,4},'Tolerance') && VARPAR{CYC,5}>OmTol) ||...
                (strcmpi(VARPAR{CYC,4},'Iteration') && VARPAR{CYC,5}<(k-klast))
            klast=k;
%             Reold=Re;
            CYC=CYC+1;
            A=VARPAR{CYC,1};
%             Re=VARPAR{CYC,2};
%             if Re~=Reold
%                 
%             end
            end
        end
        
        Omega0=Omega;
        %% Vectorized OmegaCalc
        ii=2:Nx+1;
        jj=2:My+1;
        Omega(jj,ii) = Omega0(jj,ii).*(1-alpha./(d(jj-1,ii-1).^2)) + ...
        Omega0(jj,ii+1).*(-Cx2.*u(jj,ii+1).*dip1(jj-1,ii-1) + alphaX)./(d(jj-1,ii-1).^2) + ...
        Omega0(jj,ii-1).*( Cx2.*u(jj,ii-1).*dim1(jj-1,ii-1) + alphaX)./(d(jj-1,ii-1).^2) + ...
        Omega0(jj+1,ii).*(-Cy2.*v(jj+1,ii).*djp1(jj-1,ii-1) + alphaY)./(d(jj-1,ii-1).^2) + ...
        Omega0(jj-1,ii).*( Cy2.*v(jj-1,ii).*djm1(jj-1,ii-1) + alphaY)./(d(jj-1,ii-1).^2);
        %% Vectorized PsiCalc
        PsiTol = 1;
        a=My:-1:1;
        b=My+1:-1:2;
        c=My+2:-1:3;
        dd=1:Nx;
        e=2:Nx+1;
        f=3:Nx+2;
        % Jacobian Solver Method
        while (abs(PsiTol) > Tol)
            Psi0 = Psi;
            Psi(My+1:-1:2,2:Nx+1) = KappaA.*(dxx.*Omega(b,e).*d2(a,dd) + Psi0(b,f) + ...
                Psi0(b,dd) + Kappa2.*(Psi0(c,e) + Psi0(a,e)));
            PsiTol=max(max(abs(Psi-Psi0)));
        end

        %% Vectorized Boundary Conditions
        if JETTEST
            ia = JETSTART(1);
            ib = JETEND(1);
            c0 = JETAMP(1);
            freq = JETFRQ(1);
            amewa = (ia-ceil((Nx+2)/2))*dx;
            f = sin(2*pi*freq*t);
             % Lower
             %% Non Vectorized
%             for ii=1:Nx+2
%                 if ii>=ia && ii<=ib
%                     Psi(1,ii)=(-c0*(-.5*amewa*sqrt(amewa^2+1)-.5*sinh(amewa)+.5*x(ii)*sqrt(x(ii).^2+1)+.5*sinh(x(ii))))*f;
%                 
%                 elseif ii>ib
%                     Psi(1,ii)=Psi(1,ib);
%                 else
%                     Psi(1,ii)=0;
%                 end
%                 Omega(1,ii)=(7*Psi(1,ii)-8*Psi(2,ii)+Psi(3,ii))/(2*dyy)/d6(ii);
%                 u(1,ii)=0;
%                 v(1,ii)=0;
%                 if ii>ia && ii<ib
%                     v(1,ii)=c0*f;
%                 end
%                 if ii>=(ia-1) && ii<=(ib+1)
%                     Omega(1,ii)=Omega(1,ii)+v(1,ii+1)*sqrt(x(ii+1).^2+1)-v(1,ii-1)*sqrt(x(ii-1).^2+1)/(2*dx)/d6(ii);
%                 end
%                 
%             end
            %% Vectorized
            Psi(1,1:Nx+2) = 0;
            Psi(1,ia:ib) = (-c0.*(-.5*amewa*sqrt(amewa.^2+1)-.5*sinh(amewa)+.5*x(ia:ib).*sqrt(x(ia:ib).^2+1)+.5*sinh(x(ia:ib))))*f;
            Psi(1,ib+1:end) = Psi(1,ib);
            
            Omega(1,:) = (7*Psi(1,:)-8*Psi(2,:)+Psi(3,:))./(2*dyy)./d6(1,:);
            u(1,1:Nx+2) = 0;
            v(1,1:Nx+2) = 0;
            v(1,ia+1:ib-1) = c0*f;
            
            ii = ia-1:ib+1;
            Omega(1,ii) = Omega(1,ii)+v(1,ii+1).*sqrt(x(ii+1).^2+1)-v(1,ii-1).*sqrt(x(ii-1).^2+1)./(2*dx)./d6(1,ii);
        else
            Psi(1,1:Nx+2) = 0;
            Omega(1,1:Nx+2) = (7.0*Psi(1,1:Nx+2)-8.0*Psi(2,1:Nx+2)+...
                                                          Psi(3,1:Nx+2))...
                                                            /(2.0*dyy)./d6;
            u(1,1:Nx+2) = 0;
            v(1,1:Nx+2) = 0;
        end
        
        % upper
        jj=My+2;
        Omega(jj,1:Nx+2) = 0;
        Psi(jj,1:Nx+2) = (x+A)*(y(jj)-1);
        u(jj,1:Nx+2) = (x+A)./sqrt(x.*x+y(jj)*y(jj));
        v(jj,1:Nx+2) =-(y(jj)-1)./sqrt(x.*x+y(jj)*y(jj));

        % Side BC
        % jj<IBL
        ii=1;
        Omega(2:IBL,1)= Omega(2:IBL,ii+1);
        Psi(2:IBL,ii) = Psi(2:IBL,ii+1) ;
        u(2:IBL,ii) = u(2:IBL,ii+1);
        v(2:IBL,ii) = v(2:IBL,ii+1);
        ii=Nx+2;
        Omega(2:IBL,ii) = Omega(2:IBL,ii-1);
        Psi(2:IBL,ii) = Psi(2:IBL,ii-1) ;
        u(2:IBL,ii) = u(2:IBL,ii-1);
        v(2:IBL,ii) = v(2:IBL,ii-1);

        % jj>IBL
        ii=1;
        Omega(IBL+1:My+1,ii) = 0.0;
        Psi(IBL+1:My+1,ii) = (x(ii)+A)*(y(IBL+1:My+1)-1);
        u(IBL+1:My+1,ii) = (x(ii)+A)./...
                            sqrt(x(ii)*x(ii)+y(IBL+1:My+1).*y(IBL+1:My+1));
        v(IBL+1:My+1,ii) =-(y(IBL+1:My+1)-1)./...
                           sqrt(x(ii).*x(ii)+y(IBL+1:My+1).*y(IBL+1:My+1));
        ii=Nx+2;
        Omega(IBL+1:My+1,ii) = 0.0;
        Psi(IBL+1:My+1,ii) = (x(ii)+A)*(y(IBL+1:My+1)-1);
        u(IBL+1:My+1,ii) = (x(ii)+A)./...
                            sqrt(x(ii)*x(ii)+y(IBL+1:My+1).*y(IBL+1:My+1));
        v(IBL+1:My+1,ii) =-(y(IBL+1:My+1)-1)./...
                            sqrt(x(ii)*x(ii)+y(IBL+1:My+1).*y(IBL+1:My+1));
        %% Vectorized VelCalc
        jj=2:Nx+1;
        ii=2:My+1;
        u(ii,jj)=(Psi(ii+1,jj)-Psi(ii-1,jj))./(dy2)./d3(ii,jj);
        v(ii,jj) = -(Psi(ii,jj+1)-Psi(ii,jj-1))./(dx2)./d3(ii,jj);
        %% Check Max Value Change
        OmTol=max(max(abs(Omega-Omega0)));
        PsiTol=max(max(abs(Psi-Psi0)));
        
        %% Modify Filename %%
        tick=tick+1;
        if tick==psave
            if k~=Ot
                tick=0;
                ct=ct+1;
%                 disp(datestr(clock))
                [STR2]=WriteFile(Nx,My,Psi,Omega,...
                    u,v,k,kerr,IBL,Re,OmTol,...
                    PsiTol,dx,dy,dt,ct,Tol,Xmin,Xmax,Ymin,Ymax,...
                                                   A,filename,NumOut,STR2);
                set(Col2,'String',STR2)
                drawnow expose update
            end
        end
        %% Output Iteration %%
        if UVELFFT
            if k==1
                kk=1;
                usavect=1;
                UVEL=zeros(10000,Nx+2);
                UVEL(k,:)=gather(u(UVEL_Y+1,1:Nx+2))';
            else
                kk=kk+1;
                if k==usavect*10000
                    UVEL(kk,:)=gather(u(UVEL_Y+1,1:Nx+2))';
                    usavect=usavect+1;
                    kk=0;
                    for ii=1:Nx+2
                        USAVE=UVEL(:,ii);
                        if k==10000
                            save(strcat(filename,'_Xeq_',num2str(ii)),'USAVE','-ascii','-double','-tabs');
                        else
                            save(strcat(filename,'_Xeq_',num2str(ii)),'USAVE','-append','-ascii','-double','-tabs');
                        end
                    end
                    UVEL=zeros(10000,Nx+2);
                else
                    UVEL(kk,:)=gather(u(UVEL_Y+1,1:Nx+2))';
                end
            end
        end
        
        if k==kerr
            STR1=[STR1;{[num2str(k,'%i'),'   ', num2str(OmTol,...
                '%3e'),'   ' num2str(PsiTol,'%3e'),'   ',...
                                        num2str(floor(toc(time)),'%i')]}];
            set(Col1,'String',STR1)
            drawnow expose update
%             disp(datestr(clock))
%             disp([k,OmTol,PsiTol,kerr])
            kerr=kerr+report;
        end
        if k==Ot
            ct=ct+1;
            STR1=[STR1;{[filename,' Finished']}];
            STR2=[STR2;{[filename,' Finished']}];
            set(Col1,'String',STR1)
            set(Col2,'String',STR2)
            drawnow expose update
%             disp(datestr(clock))
            [STR2]=WriteFile(Nx,My,Psi,Omega,u,...
                v,k,kerr,IBL,Re,OmTol,PsiTol,dx,...
                dy,dt,ct,Tol,Xmin,Xmax,Ymin,Ymax,A,filename,NumOut,STR2);      
            if AVG_DATA
                WriteAvgFile(Nx,My,avgdata,k,kerr,IBL,Re,...
                    OmTol,PsiTol,dx,dy,dt,ct,Tol,Xmin,...
                    Xmax,Ymin,Ymax,A,filename,NumOut);
            end
        end
    end
end


function [STR2]=WriteFile(Nx,My,Psi,Omega,u,v,k,kerr,IBL,Re,OmTol,...
                        PsiTol,dx,dy,dt,ct,Tol,Xmin,Xmax,Ymin,Ymax,A,...
                        filename,NumOut,STR2) %#ok<*INUSD,*INUSL>

OUTIN=[Nx+2,My+2,k,kerr,IBL];
OUTIN(1,length(OUTIN)+1:Nx+2)=zeros(1,Nx+2-length(OUTIN));

OUTDP=[Re,OmTol,PsiTol,dx,dy,dt,Tol,Xmin,Xmax,Ymin,Ymax,A];
OUTDP(1,length(OUTDP)+1:Nx+2)=zeros(1,Nx+2-length(OUTDP));

RESULT=[Omega;Psi;u;v;OUTIN;OUTDP]; %#ok<NASGU>

fill=length(num2str(NumOut));
incc=num2str(ct);
fp=length(incc)-1;
fileend(1:fill)='0';
fileend(end-fp:end)=incc(1:end);

STR2=[STR2;{strcat(filename,fileend)}];

save(strcat(filename,fileend),'RESULT','-ascii','-double','-tabs')
end

function WriteAvgFile(Nx,My,avgdata,k,kerr,IBL,Re,OmTol,PsiTol,dx,dy,dt,...
                              ct,Tol,Xmin,Xmax,Ymin,Ymax,A,filename,NumOut)
                
    OUTIN=[Nx+2,My+2,k,kerr,IBL];
    OUTIN(1,length(OUTIN)+1:Nx+2)=zeros(1,Nx+2-length(OUTIN));

    OUTDP=[Re,OmTol,PsiTol,dx,dy,dt,Tol,Xmin,Xmax,Ymin,Ymax,A];
    OUTDP(1,length(OUTDP)+1:Nx+2)=zeros(1,Nx+2-length(OUTDP));

    RESULT=[avgdata;OUTIN;OUTDP]; %#ok<NASGU>
    save(strcat(filename,'AVGDATA'),'RESULT','-ascii','-double','-tabs')
                
end


function SNAP_EXE(STR2,YHeight,XWidth,xmin,xmax,ymin,ymax,U,V)
    done=0;
    while done==0;
        VIS1=questdlg('What would you like to see?','Visualization','Frame Preview','Velocity Profile','Cancel','Frame Preview');
        if strcmpi(VIS1,'Cancel')
            done=1;
             % Back to Computations
        else
            VIS2=questdlg('I want to analyze:','Visualization','The Current Iteration','An Output File','The Current Iteration');
            if strcmpi(VIS2,'The Current Iteration')
                if strcmpi(VIS1,'Velocity Profile')
                    list={[2:YHeight]'};
                    VIS3=listdlg('ListString',num2str(list{1}),'Name','UVEL Plot','PromptString','Select the value of Y to plot','SelectionMode','Single','InitialValue',1,'OKString','Plot','ListSize',[160 200]);
                    figure();
                    plot(1:size(U,2),U(VIS3+1,:))
                    title(['Velocity Profile @Y= ',num2str(VIS3+1)])
                    xlabel('X-Position')
                    ylabel('Velocity')
                else
                    %% FlowField
                    xskew=1;
                    yskew=1;
                    %   
                    x=zeros(1,XWidth);
                    y=zeros(1,YHeight);
                    for ii=1:XWidth
                        if (ii<=abs(xmin/(xmax-xmin))*(XWidth-1))
                            x(ii)=(-1/2)*abs(xmax-xmin)*abs((ii-((XWidth-1)/2)-1)/((XWidth-1)/2))^xskew;
                        else
                            x(ii)=(1/2)*abs(xmax-xmin)*abs((ii-((XWidth-1)/2)-1)/((XWidth-1)/2))^xskew;
                        end
                    end
                    for ii=1:YHeight
                        y(ii)=ymin+abs(ymax-ymin)*((ii-1)/(YHeight-1))^yskew;
                    end
                   
                    [mu,eta]=meshgrid(x,y);
                    x=1/2*(mu.^2-eta.^2);
                    y=mu.*eta;
                    bx=1/2*(mu(1,:).^2-1);
                    by=mu(1,:);

                    r=sqrt(x.^2+y.^2);
        %             tic;
                    dmux=zeros(YHeight,XWidth);
                    dmuy=zeros(YHeight,XWidth);
                    detx=zeros(YHeight,XWidth);
                    dety=zeros(YHeight,XWidth);

                    dmux(2:YHeight,1:round((XWidth/2))-1)=1/2.*(x(2:YHeight,1:round((XWidth/2))-1)./r(2:YHeight,1:round((XWidth/2))-1)+1)./(mu(2:YHeight,1:round((XWidth/2))-1));
                    dmux(2:YHeight,round((XWidth/2))+1:XWidth)=1/2.*(x(2:YHeight,round((XWidth/2))+1:XWidth)./r(2:YHeight,round((XWidth/2))+1:XWidth)+1)./(mu(2:YHeight,round((XWidth/2))+1:XWidth));
                    dmuy(2:YHeight,1:round((XWidth/2))-1)=1/2.*(y(2:YHeight,1:round((XWidth/2))-1)./r(2:YHeight,1:round((XWidth/2))-1))./mu(2:YHeight,1:round((XWidth/2))-1) ;
                    dmuy(2:YHeight,round((XWidth/2))+1:XWidth)=1/2.*(y(2:YHeight,round((XWidth/2))+1:XWidth)./r(2:YHeight,round((XWidth/2))+1:XWidth))./mu(2:YHeight,round((XWidth/2))+1:XWidth) ;
                    detx(2:YHeight,1:round((XWidth/2))-1)=1/2.*(x(2:YHeight,1:round((XWidth/2))-1)./r(2:YHeight,1:round((XWidth/2))-1)-1)./(eta(2:YHeight,1:round((XWidth/2))-1));
                    detx(2:YHeight,round((XWidth/2))+1:XWidth)=1/2.*(x(2:YHeight,round((XWidth/2))+1:XWidth)./r(2:YHeight,round((XWidth/2))+1:XWidth)-1)./(eta(2:YHeight,round((XWidth/2))+1:XWidth));
                    dety(2:YHeight,1:round((XWidth/2))-1)=1/2.*(y(2:YHeight,1:round((XWidth/2))-1)./r(2:YHeight,1:round((XWidth/2))-1))./eta(2:YHeight,1:round((XWidth/2))-1) ;
                    dety(2:YHeight,round((XWidth/2))+1:XWidth)=1/2.*(y(2:YHeight,round((XWidth/2))+1:XWidth)./r(2:YHeight,round((XWidth/2))+1:XWidth))./eta(2:YHeight,round((XWidth/2))+1:XWidth) ;

                    % average mu=0 line:
                    dmux(1:YHeight,round(XWidth/2))=(dmux(1:YHeight,round((XWidth/2))-1)+dmux(1:YHeight,round((XWidth/2))+1))/2;
                    dmuy(1:YHeight,round(XWidth/2))=(dmuy(1:YHeight,round((XWidth/2))-1)+dmuy(1:YHeight,round((XWidth/2))+1))/2;
                    detx(1:YHeight,round(XWidth/2))=(detx(1:YHeight,round((XWidth/2))-1)+detx(1:YHeight,round((XWidth/2))+1))/2;
                    dety(1:YHeight,round(XWidth/2))=(dety(1:YHeight,round((XWidth/2))-1)+dety(1:YHeight,round((XWidth/2))+1))/2;

                    %
                    U_new=zeros(YHeight,XWidth);
                    V_new=zeros(YHeight,XWidth);

                    U_new(2:YHeight,1:XWidth)=(U(2:YHeight,1:XWidth).*dety(2:YHeight,1:XWidth)-V(2:YHeight,1:XWidth).*dmuy(2:YHeight,1:XWidth)).*sqrt(mu(2:YHeight,1:XWidth).^2+eta(2:YHeight,1:XWidth).^2);  %Vx
                    V_new(2:YHeight,1:XWidth)=-(U(2:YHeight,1:XWidth).*detx(2:YHeight,1:XWidth)-V(2:YHeight,1:XWidth).*dmux(2:YHeight,1:XWidth)).*sqrt(mu(2:YHeight,1:XWidth).^2+eta(2:YHeight,1:XWidth).^2); %Vy 

                    clear s
                    %   Creates V-vctr from Vx and Vy
                    s=zeros(size(U_new,1),XWidth);
                    s(1:size(U_new,1),1:round(XWidth/2))=sqrt(U_new(1:size(U_new,1),1:round(XWidth/2)).^2+V_new(1:size(U_new,1),1:round(XWidth/2)).^2);
                    s(1:size(U_new,1),round(XWidth/2):XWidth)=sqrt(U_new(1:size(U_new,1),round(XWidth/2):XWidth).^2+V_new(1:size(U_new,1),round(XWidth/2):XWidth).^2).*sign(U_new(1:size(U_new,1),round(XWidth/2):XWidth));
        %             toc
                    hf=figure('Units','pixels','visible','off','Color','w');
                    hfaxis=axes();
                    set(hf,'CurrentAxes',hfaxis)
                    hold on
                    plot(bx,by,'k','linewidth',1) %contour plot
                    contourf(x,y,s,ceil(YHeight/2),'linestyle','none'); 
                    colormap(jet(round(YHeight/2)))
                    scale=.15;
                    quiver(x,y,U_new,V_new,scale,'color','k','ShowArrowHead','off')
                    axis('equal')
                    axis([-2 16 -2 8])
                    set(hf,'visible','on')
                    %% End Flowfield
                end
            else
                %% Load File
                VIS5=listdlg('ListString',{STR2{2:end}},'Name','PFI','PromptString','Select Data File','SelectionMode','Single');
                fid=fopen(STR2{VIS5+1});
                FlowData2=textscan(fid,'%f',(XWidth*YHeight*4+XWidth*2));
                fclose(fid);
                FlowData2=cell2mat(FlowData2);
                kk=0;
                FlowData=zeros(YHeight*4+2,XWidth);
                for ii=1:YHeight*4+2;
                    for jj=1:XWidth
                        kk=kk+1;
                        FlowData(ii,jj)=FlowData2(kk);
                    end
                end      
                % Delete the last two rows from Data file to give us only lab data
                for delete_rows=1:2
                    FlowData(size(FlowData,1),:)=[];
                end
                %
                U=FlowData(2*YHeight+1:3*YHeight,:);
                V=FlowData(3*YHeight+1:4*YHeight,:);
                %% End Load
                if strcmpi(VIS1,'Velocity Profile')
                    list={[2:YHeight]'};
                    VIS3=listdlg('ListString',num2str(list{1}),'Name','UVEL Plot','PromptString','Select the value of Y to plot','SelectionMode','Single','InitialValue',1,'OKString','Plot','ListSize',[160 200]);
                    figure();
                    plot(1:size(U,2),U(VIS3+1,:))
                    title(['Velocity Profile @Y= ',num2str(VIS3+1)])
                    xlabel('X-Position')
                    ylabel('Velocity')
                else
                                        %% FlowField
                    xskew=1;
                    yskew=1;
                    %   
                    x=zeros(1,XWidth);
                    y=zeros(1,YHeight);
                    for ii=1:XWidth
                        if (ii<=abs(xmin/(xmax-xmin))*(XWidth-1))
                            x(ii)=(-1/2)*abs(xmax-xmin)*abs((ii-((XWidth-1)/2)-1)/((XWidth-1)/2))^xskew;
                        else
                            x(ii)=(1/2)*abs(xmax-xmin)*abs((ii-((XWidth-1)/2)-1)/((XWidth-1)/2))^xskew;
                        end
                    end
                    for ii=1:YHeight
                        y(ii)=ymin+abs(ymax-ymin)*((ii-1)/(YHeight-1))^yskew;
                    end
                   
                    [mu,eta]=meshgrid(x,y);
                    x=1/2*(mu.^2-eta.^2);
                    y=mu.*eta;
                    bx=1/2*(mu(1,:).^2-1);
                    by=mu(1,:);

                    r=sqrt(x.^2+y.^2);
        %             tic;
                    dmux=zeros(YHeight,XWidth);
                    dmuy=zeros(YHeight,XWidth);
                    detx=zeros(YHeight,XWidth);
                    dety=zeros(YHeight,XWidth);

                    dmux(2:YHeight,1:round((XWidth/2))-1)=1/2.*(x(2:YHeight,1:round((XWidth/2))-1)./r(2:YHeight,1:round((XWidth/2))-1)+1)./(mu(2:YHeight,1:round((XWidth/2))-1));
                    dmux(2:YHeight,round((XWidth/2))+1:XWidth)=1/2.*(x(2:YHeight,round((XWidth/2))+1:XWidth)./r(2:YHeight,round((XWidth/2))+1:XWidth)+1)./(mu(2:YHeight,round((XWidth/2))+1:XWidth));
                    dmuy(2:YHeight,1:round((XWidth/2))-1)=1/2.*(y(2:YHeight,1:round((XWidth/2))-1)./r(2:YHeight,1:round((XWidth/2))-1))./mu(2:YHeight,1:round((XWidth/2))-1) ;
                    dmuy(2:YHeight,round((XWidth/2))+1:XWidth)=1/2.*(y(2:YHeight,round((XWidth/2))+1:XWidth)./r(2:YHeight,round((XWidth/2))+1:XWidth))./mu(2:YHeight,round((XWidth/2))+1:XWidth) ;
                    detx(2:YHeight,1:round((XWidth/2))-1)=1/2.*(x(2:YHeight,1:round((XWidth/2))-1)./r(2:YHeight,1:round((XWidth/2))-1)-1)./(eta(2:YHeight,1:round((XWidth/2))-1));
                    detx(2:YHeight,round((XWidth/2))+1:XWidth)=1/2.*(x(2:YHeight,round((XWidth/2))+1:XWidth)./r(2:YHeight,round((XWidth/2))+1:XWidth)-1)./(eta(2:YHeight,round((XWidth/2))+1:XWidth));
                    dety(2:YHeight,1:round((XWidth/2))-1)=1/2.*(y(2:YHeight,1:round((XWidth/2))-1)./r(2:YHeight,1:round((XWidth/2))-1))./eta(2:YHeight,1:round((XWidth/2))-1) ;
                    dety(2:YHeight,round((XWidth/2))+1:XWidth)=1/2.*(y(2:YHeight,round((XWidth/2))+1:XWidth)./r(2:YHeight,round((XWidth/2))+1:XWidth))./eta(2:YHeight,round((XWidth/2))+1:XWidth) ;

                    % average mu=0 line:
                    dmux(1:YHeight,round(XWidth/2))=(dmux(1:YHeight,round((XWidth/2))-1)+dmux(1:YHeight,round((XWidth/2))+1))/2;
                    dmuy(1:YHeight,round(XWidth/2))=(dmuy(1:YHeight,round((XWidth/2))-1)+dmuy(1:YHeight,round((XWidth/2))+1))/2;
                    detx(1:YHeight,round(XWidth/2))=(detx(1:YHeight,round((XWidth/2))-1)+detx(1:YHeight,round((XWidth/2))+1))/2;
                    dety(1:YHeight,round(XWidth/2))=(dety(1:YHeight,round((XWidth/2))-1)+dety(1:YHeight,round((XWidth/2))+1))/2;

                    %
                    U_new=zeros(YHeight,XWidth);
                    V_new=zeros(YHeight,XWidth);

                    U_new(2:YHeight,1:XWidth)=(U(2:YHeight,1:XWidth).*dety(2:YHeight,1:XWidth)-V(2:YHeight,1:XWidth).*dmuy(2:YHeight,1:XWidth)).*sqrt(mu(2:YHeight,1:XWidth).^2+eta(2:YHeight,1:XWidth).^2);  %Vx
                    V_new(2:YHeight,1:XWidth)=-(U(2:YHeight,1:XWidth).*detx(2:YHeight,1:XWidth)-V(2:YHeight,1:XWidth).*dmux(2:YHeight,1:XWidth)).*sqrt(mu(2:YHeight,1:XWidth).^2+eta(2:YHeight,1:XWidth).^2); %Vy 

                    clear s
                    %   Creates V-vctr from Vx and Vy
                    s=zeros(size(U_new,1),XWidth);
                    s(1:size(U_new,1),1:round(XWidth/2))=sqrt(U_new(1:size(U_new,1),1:round(XWidth/2)).^2+V_new(1:size(U_new,1),1:round(XWidth/2)).^2);
                    s(1:size(U_new,1),round(XWidth/2):XWidth)=sqrt(U_new(1:size(U_new,1),round(XWidth/2):XWidth).^2+V_new(1:size(U_new,1),round(XWidth/2):XWidth).^2).*sign(U_new(1:size(U_new,1),round(XWidth/2):XWidth));
        %             toc
                    hf=figure('Units','pixels','visible','off','Color','w');
                    hfaxis=axes();
                    set(hf,'CurrentAxes',hfaxis)
                    hold on
                    plot(bx,by,'k','linewidth',1) %contour plot
                    contourf(x,y,s,ceil(YHeight/2),'linestyle','none'); 
                    colormap(jet(round(YHeight/2)))
                    scale=1;
                    quiver(x,y,U_new,V_new,scale,'color','k','ShowArrowHead','off')
                    axis('equal')
                    axis([-2 16 -2 8])
                    set(hf,'visible','on')
                    %% End Flowfield

                end
            end
            VIS4=questdlg('What would you like to do?','PFI','Back To Computations','Make Another Plot','Back To Computations');
            if strcmpi(VIS4,'Back To Computations')
                done=1;
            end
        end
    end
end