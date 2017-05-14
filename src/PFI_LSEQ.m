%%
clear
close all
close('all');
clc
%% Geometry Definition
NX = 401;
NY = 1601;
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
ii=2:NX-1;
jj=2:NY-1;
d = gpuArray(sqrt(xx.^2 + yy.^2));
di = 1./d;
d2 = d.^2;
di_2dx = di(jj,ii) ./ (2*dx);
di_2dy = di(jj,ii) ./ (2*dy);


dip1 = gpuArray(sqrt(xx(jj,ii+1).^2+yy(jj,ii).^2));
dim1 = gpuArray(sqrt(xx(jj,ii-1).^2+yy(jj,ii).^2));
djp1 = gpuArray(sqrt(xx(jj,ii).^2+yy(jj+1,ii).^2));
djm1 = gpuArray(sqrt(xx(jj,ii).^2+yy(jj-1,ii).^2));

D2 = gpuArray(xx.^2+yy.^2);
D3(jj,ii)=sqrt(xx(jj,ii).^2+yy(jj,ii).^2);

d6 = x.^2+y(1)^2;

Kappa2=(dx/dy)^2;
KappaA = 1.0/(2.0*(1+Kappa2));

% Compute Mapping: Parabolic -> Cartesian
pCOORD.NX = NX;
pCOORD.NY = NY;
pCOORD.xMin = xMin;
pCOORD.xMax = xMax;
pCOORD.yMin = yMin;
pCOORD.yMax = yMax;
pCOORD.xx = xx;
pCOORD.yy = yy;
cCOORD = Para2Cart(pCOORD);

%% Flow Definition
Re = 700;
Visc = 1/Re;
Rc = Re*dx;
IBL = ceil(.05*NY);%40;
A = 2;
%% Solution Control
tol = 1e-5;
tmax = 10;
dt_max = 1/(10*Re);
method = 'krylov';

%% Results Options
RES.OutputData = true;
dBuffer = 1e2;
if RES.OutputData == true
    RES.VelProfile{1}.saveData = true;
    RES.VelProfile{1}.file = matfile('SurfaceVelocity.mat','Writable',true);
    RES.VelProfile{1}.ii = 1:NX;
    RES.VelProfile{1}.jj = 2;
    RES.VelProfile{1}.DATA = zeros(length(RES.VelProfile{1}.jj),length(RES.VelProfile{1}.ii),dBuffer);
    RES.VelProfile{1}.time = zeros(1,dBuffer);
    
    RES.VelProfile{2}.saveData = true;
    RES.VelProfile{2}.file = matfile('U-VelocityProfile.mat','Writable',true);
    RES.VelProfile{2}.ii = ceil(NX/2);
    RES.VelProfile{2}.jj = 1:2*IBL;
    RES.VelProfile{2}.DATA = zeros(length(RES.VelProfile{2}.jj),length(RES.VelProfile{2}.ii),dBuffer);
    RES.VelProfile{2}.time = zeros(1,dBuffer);
    
    RES.VelProfile{3}.saveData = true;
    RES.VelProfile{3}.file = matfile('V-VelocityProfile.mat','Writable',true);
    RES.VelProfile{3}.ii = ceil(NX/2);
    RES.VelProfile{3}.jj = 1:2*IBL;
    RES.VelProfile{3}.DATA = zeros(length(RES.VelProfile{3}.jj),length(RES.VelProfile{3}.ii),dBuffer);
    RES.VelProfile{3}.time = zeros(1,dBuffer);
end
dCycle = 0;
dtOUT = 1e-2;
outNext = 0;
tNext = 0;

%% Visualization Options
VIZ.plotStateVar = false;
VIZ.plotFlow = true;
dtPLOT = 0.1;
FIG = [];

%% Movie Options
MOV.plotStateVar = false;
MOV.plotFlow = false;

%% Image Options
IMG.plotStateVar = false;
IMG.plotFlow = true;

%% Logging Options
logName = 'PFI_LSEQ.log';
fLOG = fopen(logName,'w+');
fprintf(fLOG,'DIRECT NUMERICAL SIMULATION OF PARABOLA FLOW via LINEAR SYSTEM OF EQUATIONS\n');
fprintf(fLOG,['START TIME: ' datestr(clock,31) '\n']);
fprintf(fLOG,['REYNOLDS: ' num2str(Re) '\n']);
fprintf(fLOG,['CIRCULATION: ' num2str(A) '\n']);
fprintf(fLOG,['NX: ' num2str(NX) '\n']);
fprintf(fLOG,['NY: ' num2str(NY) '\n']);
fprintf(fLOG,['METHOD: ' method '\n']);
fprintf(fLOG,['MAX SIMULATION TIME: ' num2str(tmax) '\n']);
fclose(fLOG);
%% Initialization
tstart = tic;
if strcmpi(method,'jacobi')
    FDM = [];
    F = [];
else
    % Assemble the Coefficient Matrix
    disp('Assembling the Coefficient Matrix')
    Atic = tic;
    FDM = -assembleCoeffMat(NX,NY,dx,dy);
    toc(Atic)
    if strcmpi(method,'linfactor')
        disp('LINFACTOR Method Chosen'); disp('Performing Initial Matrix Factorization');
        FAtic = tic;
        F = linfactor(FDM);
        toc(FAtic)
    elseif strcmpi(method,'krylov')
        FDM = gpuArray(FDM);
        disp('Computing Krylov Preconditioners')
        Ptic = tic;
%         [M1,M2] = ilu(gather(FDM),struct('type','ilutp','droptol',1e-3));
%         M1 = gpuArray(Lp);
%         M2 = gpuArray(Up);
        M1 = ichol(gather(FDM),struct('michol','on'));
        M1 = gather(M1);
        toc(Ptic);
        F = struct('M1',[],'M2',[]);
    else
        F = [];
    end
end
% Initialize State Variables
[PSI] = InitializePSI(A,x,y,NX,NY);
OMEGA = InitializeOMEGA(PSI,NX,NY,dyy,d6); d6 = gpuArray(d6);
[U,V] = InitializeVELOCITY(A,x,y,NX,NY);

PSI = gpuArray(PSI);
OMEGA = gpuArray(OMEGA);
U = gpuArray(U);
V = gpuArray(V);

% x = gpuArray(x);
% y = gpuArray(y);

t=0.0;

dIter = 0;
tstep = 0;
pCycle = 0;
ii = 2:NX-1;
jj = 2:NY-1;
II = gpuArray(repmat(1:NX,NY,1));
JJ = gpuArray(repmat([1:NY]',1,NX));
while t < tmax % start the time integration
    %% Compute dt
    [dt,alphaX,alphaY,alpha] = delta_t(U,V,dx,dy,dt_max,Re);
    tstep = tstep + 1;
    
    %% Omega Calc
    OMEGA(jj,ii) = arrayfun(@OMEGA_CALC,OMEGA(jj,ii),OMEGA(jj,ii+1),OMEGA(jj,ii-1),OMEGA(jj+1,ii),OMEGA(jj-1,ii),U(jj,ii+1),U(jj,ii-1),V(jj+1,ii),V(jj-1,ii),d2(jj,ii),djp1,djm1,dip1,dim1,alpha,alphaX,alphaY,dx,dy,dt);
    % Boundary Conditions
    OMEGA = OMEGA_BC(PSI,OMEGA,NX,NY,IBL,dy,d6);
    
    %% PsiCalc
    PSI = PSI_CALC(FDM,F,OMEGA,PSI,D2,dx,dy,dxx,KappaA,Kappa2,NX,NY,ii,jj,II,JJ,tol,tstep,method);
    % Boundary Conditions
    PSI = PSI_BC(PSI,x,y,NX,NY,IBL,A);
    
    %% Velocity Calc
    [U(jj,ii),V(jj,ii)] = arrayfun(@VEL_CALC,PSI(jj+1,ii),PSI(jj-1,ii),PSI(jj,ii+1),PSI(jj,ii-1),di_2dx,di_2dy); %di(jj,ii),dx,dy);
    % Boundary Conditions
    [U,V] = VELOCITY_BC(U,V,x,y,NX,NY,IBL,A);
    
    %%
    t=t+dt; % print out t
    
    %% Data Output
    if RES.OutputData == true
        [pU,pV,VEL] = PARA_VEL(U,V,cCOORD,pCOORD);
        VEL = gather(VEL);
        pU = gather(pU);
        pV = gather(pV);
        if t > outNext || tstep == 1 || t >= tmax
            outNext = outNext + dtOUT;
            dIter = dIter + 1;
            RES.VelProfile{1}.DATA(:,:,dIter) = pU(RES.VelProfile{1}.jj,RES.VelProfile{1}.ii);
            RES.VelProfile{1}.time(dIter) = t;
            RES.VelProfile{2}.DATA(:,:,dIter) = pU(RES.VelProfile{2}.jj,RES.VelProfile{2}.ii);
            RES.VelProfile{2}.time(dIter) = t;
            RES.VelProfile{3}.DATA(:,:,dIter) = pV(RES.VelProfile{3}.jj,RES.VelProfile{3}.ii);
            RES.VelProfile{3}.time(dIter) = t;
        end
        
        if (dIter == dBuffer || t >= tmax) && RES.OutputData == true
            dCycle = dCycle+1;
            dIter = 0;
            OutputData(cCOORD,pCOORD,RES,t,dCycle);
            RES.VelProfile{1}.DATA = zeros(length(RES.VelProfile{1}.jj),length(RES.VelProfile{1}.ii),dBuffer);
            RES.VelProfile{1}.time = zeros(1,dBuffer);
            RES.VelProfile{2}.DATA = zeros(length(RES.VelProfile{2}.jj),length(RES.VelProfile{2}.ii),dBuffer);
            RES.VelProfile{2}.time = zeros(1,dBuffer);
            RES.VelProfile{3}.DATA = zeros(length(RES.VelProfile{3}.jj),length(RES.VelProfile{3}.ii),dBuffer);
            RES.VelProfile{3}.time = zeros(1,dBuffer);
        end
    end
    
    %% Plot
    if  t > tNext || tstep == 1 || t >= tmax
        tNext = tNext + dtPLOT;
        pCycle = pCycle + 1;
        tElapsed = toc(tstart);
        disp(['Timestep = ' num2str(tstep) ' Time = ' num2str(t) ' Time Elapsed = ' num2str(tElapsed)])
        fLOG = fopen(logName,'a');
        fprintf(fLOG, ['Timestep = ' num2str(tstep) ' Time = ' num2str(t) ' Time Elapsed = ' num2str(tElapsed) '\r']);
        fclose(fLOG);
        [FIG,MOV] = VizFlowField(cCOORD,pCOORD,PSI,OMEGA,U,V,FIG,VIZ,MOV,IMG,pCycle);
    end
end
%%
if MOV.plotStateVar == true
    close(MOV.mov1);
end
if MOV.plotFlow == true
    close(MOV.mov2);
end
%%
STATE.OMEGA = gather(OMEGA);
STATE.PSI = gather(PSI);
STATE.U = gather(U);
STATE.V = gather(V);
%
FLOW.Re = Re;
FLOW.A = A;
FLOW.IBL = IBL;
FLOW.t = t;
%
LSEQ.F = F;
LSEQ.FDM = FDM;
save('PFI_Final.mat','STATE','FLOW','pCOORD','cCOORD','LSEQ','-v7.3')