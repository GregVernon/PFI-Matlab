%%
clear
close all
clc
%% Geometry Definition
NX = 201;
NY = 401;
xMin = -20;
xMax = 20;
yMin = 1;
yMax = 11;
% Define Parabolic Coordinates
x = linspace(xMin,xMax,NX);
y = linspace(yMin,yMax,NY);
dx = x(2) - x(1);
dx2=2*dx;
dxx=dx*dx;

dy = y(2) - y(1);
dy2=2*dy;
dyy=dy*dy;

[xx,yy] = meshgrid(x,y);
% Compute Mapping "Helper" Variables
d = sqrt(xx.^2 + yy.^2);
di = 1./d;
d2=d.^2;

ii=2:NX-1;
jj=2:NY-1;
dip1 = sqrt(xx(jj,ii+1).^2+yy(jj,ii).^2);
dim1 = sqrt(xx(jj,ii-1).^2+yy(jj,ii).^2);
djp1 = sqrt(xx(jj,ii).^2+yy(jj+1,ii).^2);
djm1 = sqrt(xx(jj,ii).^2+yy(jj-1,ii).^2);

D2=xx.^2+yy.^2;
D3(jj,ii)=sqrt(xx(jj,ii).^2+yy(jj,ii).^2);

d6=x.^2+y(1)^2;

Kappa2=(dx/dy)^2;
KappaA = 1.0/(2.0*(1+Kappa2));

% Compute Mapping: Parabolic -> Cartesian
pCOORD.NX = NX;
pCOORD.NY = NY;
pCOORD.xMin = xMin;
pCOORD.xMax = xMax;
pCOORD.yMin = yMin;
pCOORD.yMax = yMax;
cCOORD = Para2Cart(pCOORD);

%% Flow Definition
Re = 700;
Visc=1/Re;
Rc = Re*dx;
IBL = ceil(.05*NY);%40;
A = 2;
%% Solution Control
tol = 1e-5;
tmax = 100;
dt_max = 1/(10*Re);
method = 'linfactor';

%% Results Options
RES.OutputData = false;
dBuffer = 1e2;
if RES.OutputData == true
    RES.VelProfile{1}.saveData = true;
    RES.VelProfile{1}.fileName = 'SurfaceVelocity.mat';
    RES.VelProfile{1}.ii = 1:NX;
    RES.VelProfile{1}.jj = 2;
    RES.VelProfile{1}.DATA = zeros(length(RES.VelProfile{1}.jj),length(RES.VelProfile{1}.ii),dBuffer);
    
    RES.VelProfile{2}.saveData = true;
    RES.VelProfile{2}.fileName = 'VelocityProfile.mat';
    RES.VelProfile{2}.ii = ceil(NX/2);
    RES.VelProfile{2}.jj = 1:2*IBL;
    RES.VelProfile{2}.DATA = zeros(length(RES.VelProfile{2}.jj),length(RES.VelProfile{2}.ii),dBuffer);
end
dCycle = 0;

dtOUT = 0.01;
tNext = 0;
MOV.plotStateVar = false;
MOV.plotFlow = true;

%% Visualization Options
VIZ.plotStateVar = false;
VIZ.plotFlow = true;
FIG = [];
%% Initialization
tstart = tic;
if strcmpi(method,'jacobi')
    %     continue
else
    % Assemble the Coefficient Matrix
    
    disp('Assembling the Coefficient Matrix')
    Atic = tic;
    FDM=assembleCoeffMat(NX,NY,dx,dy);
    toc(Atic)
end
% Initialize State Variables
[PSI] = InitializePSI(A,x,y,NX,NY);
OMEGA = InitializeOMEGA(PSI,NX,NY,dyy,d6);
[U,V] = InitializeVELOCITY(A,x,y,NX,NY);

t=0.0;

pIter = 0;
dIter = 0;
tstep = 0;
ii=2:NX-1;
jj=2:NY-1;
while t < tmax % start the time integration
    %% Compute dt
    [dt,alphaX,alphaY,alpha] = delta_t(U,V,dx,dy,dt_max,Re);
    pIter = pIter + 1;
    dIter = dIter + 1;
    tstep = tstep + 1;
    
    %% Omega Calc
    %OMEGA(jj,ii) = OMEGA(jj,ii) + dt*(-((PSI(jj+1,ii)-PSI(jj-1,ii)).*(OMEGA(jj,ii+1)-OMEGA(jj,ii-1))) + ((PSI(jj,ii+1)-PSI(jj,ii-1)).*(OMEGA(jj+1,ii)-OMEGA(jj-1,ii)))./(4*h^2) + Visc*(OMEGA(jj,ii+1)+OMEGA(jj,ii-1)+OMEGA(jj+1,ii)+OMEGA(jj-1,ii)-4.0*OMEGA(jj,ii))./(h*h)); % vorticity
    OMEGA(jj,ii) = OMEGA(jj,ii).*(1-alpha./(d(jj,ii).^2)) + ...
        OMEGA(jj,ii+1).*(-(dt/(2*dx)).*U(jj,ii+1).*dip1 + alphaX)./(d(jj,ii).^2) + ...
        OMEGA(jj,ii-1).*( (dt/(2*dx)).*U(jj,ii-1).*dim1 + alphaX)./(d(jj,ii).^2) + ...
        OMEGA(jj+1,ii).*(-(dt/(2*dy)).*V(jj+1,ii).*djp1 + alphaY)./(d(jj,ii).^2) + ...
        OMEGA(jj-1,ii).*( (dt/(2*dy)).*V(jj-1,ii).*djm1 + alphaY)./(d(jj,ii).^2);
    
    %% PsiCalc
    
    if strcmpi(method,'jacobi')
        %% Vectorized PsiCalc %%
        PsiTol = 1;
        kPsi = 0;
        while PsiTol > tol
            kPsi = kPsi+1;
            PSI0 = PSI;
            PSI2 = KappaA.*(dxx.*OMEGA(jj,ii).*D2(jj,ii) + PSI0(jj,ii+1) + ...
                PSI0(jj,ii-1) + Kappa2.*(PSI0(jj+1,ii) + PSI0(jj-1,ii)));
            
            PSI(jj,ii) = PSI2;
            PsiTol = max(max(abs(PSI-PSI0)));
        end
    elseif strcmpi(method,'direct')
        b = computeRHS(NX,NY,x,y,dx,dy,D2,OMEGA,PSI,A,IBL);
        X = FDM\b;
        PSI(jj,ii) = reshape(X',NY-2,NX-2);
        PSI = PSI_BC(PSI,x,y,NX,NY,IBL,A);
    elseif strcmpi(method,'linfactor')
        if tstep == 1
            disp('LINFACTOR Method Chosen'); disp('Performing Initial Matrix Factorization');
            FAtic = tic;
            F = linfactor(FDM);
            toc(FAtic)
        end
        b = computeRHS(NX,NY,x,y,dx,dy,D2,OMEGA,PSI,A,IBL);
        X = linfactor(F,b);
        PSI(jj,ii) = reshape(X',NY-2,NX-2);
        PSI = PSI_BC(PSI,x,y,NX,NY,IBL,A);
    elseif strcmpi(method,'inverseExplicit')
        if tstep == 1
            disp('Inverse Explicit Method Chosen'); disp('Performing Initial Matrix Inversion');
            mItic = tic;
            invA = full(A\speye(size(FDM)));
            toc(mItic)
        end
        b = computeRHS(NX,NY,x,y,dx,dy,D2,OMEGA,PSI,A,IBL);
        X = invA * b;
        PSI(jj,ii) = reshape(X,NY-2,NX-2);
        PSI = PSI_BC(PSI,x,y,NX,NY,IBL,A);
    elseif strcmpi(method,'iterative')
        if tstep == 1
            disp('Iterative Method Chosen'); disp('Computing Preconditioners');
            PAtic = tic;
            %             [Lp,Up] = ilu(FDM,struct('type','ilutp','droptol',1e-5));
            toc(PAtic)
        end
        b = computeRHS(NX,NY,x,y,dx,dy,D2,OMEGA,PSI,A,IBL);
        X0 = reshape(PSI(jj,ii),(NY-2)*(NX-2),1);
        [X,flag] = gmres(FDM,b,100,tol,1e5,[],[],X0);
        %         [X,flag] = pcg(FDM,b,tol,1e5,Lp,Up,X0);
        PSI(jj,ii) = reshape(X',NY-2,NX-2);
        PSI = PSI_BC(PSI,x,y,NX,NY,IBL,A);
    end
    %     X = M1*b;
    
    
    %% Boundary Conditions
    PSI = PSI_BC(PSI,x,y,NX,NY,IBL,A);
    OMEGA = OMEGA_BC(PSI,OMEGA,NX,NY,IBL,dy,d6);
    %     [U,V] = VELOCITY_BC(U,V,x,y,NX,NY,IBL,A);
    
    %% Velocity Calc
    U(jj,ii) = di(jj,ii) .* (PSI(jj+1,ii) - PSI(jj-1,ii)) ./ (2*dy);% ./ D3(jj,ii);
    V(jj,ii) = di(jj,ii) .* -(PSI(jj,ii+1) - PSI(jj,ii-1)) ./ (2*dx);% ./ d(jj,ii);
    [U,V] = VELOCITY_BC(U,V,x,y,NX,NY,IBL,A);
    %%
    t=t+dt; % print out t
    %% Plot
    if  t > tNext || tstep == 1 || t >= tmax %|| pIter == 1e3
        tNext = tNext + dtOUT;
        disp(['Timestep = ' num2str(tstep) ' Time = ' num2str(t) ' Time Elapsed = ' num2str(toc(tstart))])
        pIter = 0;
        [FIG,MOV] = VizFlowField(cCOORD,pCOORD,U,V,FIG,VIZ,MOV);
    end
    %% Data Output
    if RES.OutputData == true
        [pU,pV,VEL] = PARA_VEL(U,V,cCOORD,pCOORD);
        RES.VelProfile{1}.DATA(:,:,dIter) = VEL(RES.VelProfile{1}.jj,RES.VelProfile{1}.ii);
        RES.VelProfile{2}.DATA(:,:,dIter) = VEL(RES.VelProfile{2}.jj,RES.VelProfile{2}.ii);
        
        if dIter == dBuffer || tstep == 1 || t >= tmax && RES.OutputData == true
            dCycle = dCycle+1;
            dIter = 0;
            OutputData(cCOORD,pCOORD,RES,t,dCycle);
            RES.VelProfile{1}.DATA = zeros(length(RES.VelProfile{1}.jj),length(RES.VelProfile{1}.ii),dBuffer);
            RES.VelProfile{2}.DATA = zeros(length(RES.VelProfile{2}.jj),length(RES.VelProfile{2}.ii),dBuffer);
        end
    end
end
%%
if MOV.plotStateVar == true
    close(MOV.mov1);
end
if MOV.plotFlow == true
    close(MOV.mov2);
end
