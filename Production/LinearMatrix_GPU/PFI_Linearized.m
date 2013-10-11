clear
close all
clc
%% Inputs
Re = 4200;
A = 1.0;
dt = 0.0002380;
Nx = 1000;
My = 1000;
IBL = 100;
xmin = -20;
xmax = 20;
ymin = 1;
ymax = 11;
tol = 1e-5;
end_time = 100000;
log_update = 100;
output_write = 1000;

filename = 'PFI_Tests';
WriteVars = {'Omega' ',' 'Psi' ',' 'U' ',' 'V'};
formats = {'ASCII'};
%% Generate Grid
x = linspace(xmin,xmax,Nx);
dx = x(2) - x(1);
dx2 = 2 * dx;
dxx = dx^2;

y = linspace(ymin,ymax,My);
dy = y(2) - y(1);
dy2 = 2 * dy;
dyy = dy^2;


% !!!!!!!!!! COMMENT HERE !!!!!
Kappa2=(dx/dy)^2;
KappaA = 1.0/(2.0*(1+Kappa2));
Rc = Re*dx;

Cx = dt/dx;
Cx2 = .5*Cx;
Cy = dt/dy;
Cy2 = .5*Cy;

C = max(Cx,Cy);
alphaX = dt/(dxx*Re);
alphaY = dt/(dyy*Re);
alpha = 2*alphaX + 2*alphaY;


% !!!!!!!!!! COMMENT HERE !!!!!
xx = repmat(x,My,1);
yy = repmat(y',1,Nx);

% !!!!!!!!!! COMMENT HERE !!!!!
% ii=2:Nx-1;
% jj=2:My-1;
% d(1:My-2,1:Nx-2)=sqrt(xx(jj,ii).^2+yy(jj,ii).^2);
% d2(1:My-2,1:Nx-2)=xx(2:My-1,2:Nx-1).^2+yy(2:My-1,2:Nx-1).^2;
% d3(jj,ii)=sqrt(xx(jj,ii).^2+yy(jj,ii).^2);
% d6=x.^2+y(1)^2;

d = sqrt(xx.^2 + yy.^2);
d2 = xx.^2 + yy.^2;
d3 = sqrt(xx.^2 + yy.^2);
d6 = x.^2 + y(1)^2;


% !!!!!!!!!! COMMENT HERE !!!!!
ii = 1:Nx;
jj = 1:My;
iimat = repmat(ii,My,1);
jjmat = repmat(jj',1,Nx);
% make temporary xx2 and yy2 with ghost points
xx2 = repmat([x(1)-dx x x(end)+dx],My+2,1);
yy2 = repmat([y(1)-dy y y(end)+dy]',1,Nx+2);
dip1 = sqrt(xx2(jj+1,ii+2).^2 + yy2(jj+1,ii+1).^2);
dim1 = sqrt(xx2(jj+1,ii).^2 + yy2(jj+1,ii+1).^2);
djp1 = sqrt(xx2(jj+1,ii+1).^2 + yy2(jj+2,ii+1).^2);
djm1 = sqrt(xx2(jj+1,ii+1).^2 + yy2(jj,ii+1).^2);
clear xx2 yy2
% ii = 2:Nx-1;
% jj = 2:My-1;
% dip1(1:My-2,1:Nx-2) = sqrt(xx(jj,ii+1).^2+yy(jj,ii).^2);
% dim1(1:My-2,1:Nx-2) = sqrt(xx(jj,ii-1).^2+yy(jj,ii).^2);
% djp1(1:My-2,1:Nx-2) = sqrt(xx(jj,ii).^2+yy(jj+1,ii).^2);
% djm1(1:My-2,1:Nx-2) = sqrt(xx(jj,ii).^2+yy(jj-1,ii).^2);

%% INITIALIZE FIELD VARIABLES
V=zeros(My,Nx);
U=zeros(My,Nx);
Psi=zeros(My,Nx);
Omega=zeros(My,Nx);

for jj = 1:Nx
    for ii = 1:My
        V(ii,jj) = -(y(ii)-1)/sqrt(x(jj)^2+y(ii)^2);
    end
end

for jj = 1:Nx
    for ii = 2:My
        U(ii,jj) = (x(jj)+A)/sqrt(x(jj)^2+y(ii)^2);
    end
end

for jj = 1:Nx
    for ii = 1:My
        Psi(ii,jj) = (x(jj)+A)*(y(ii)-1);
    end
end

for ii = 1:Nx
    Omega(1,ii) = (7.0*Psi(1,ii)-8.0*Psi(2,ii)+Psi(3,ii)) / (2.0*dyy)/d6(1,ii);
end
%% GPU ARRAYS
Omega = gpuArray(Omega);
Psi = gpuArray(Psi);
U = gpuArray(U);
V = gpuArray(V);
xx = gpuArray(xx);
yy = gpuArray(yy);
d = gpuArray(d);
dip1 = gpuArray(dip1);
dim1 = gpuArray(dim1);
djp1 = gpuArray(djp1);
djm1 = gpuArray(djm1);


%% GENERATE SPARSE FINITE DIFFERENCE MATRIX
disp('Generating FD Coefficient Matrix...')
[AMat] = FDCM(Nx,My,dx,dy);
disp('Done')

%% GENERATE PRECONDITIONERS IF APPLICABLE
[M1,M2] = ilu(AMat);

%% INTERIOR LOCATION -- MATRIX->VECTOR
Z = false(Nx*My,1);
kk = 0;
for ii = 1:Nx
    for jj = 1:My
        kk = kk+1;
        if ii ~= 1 && ii ~=Nx && jj ~= 1 && jj ~= My
           Z(kk) = true;
        else 
           Z(kk) = false;
        end
    end
end
%% INITIALIZE COUNTER VARIABLES
timestep = 0;
output_counter = 0;
log_counter = 0;
% iter = zeros(end_time,1);
%% CFD BEGIN
while timestep < end_time
    % Step Counters
    timestep = timestep + 1;
    output_counter = output_counter + 1;
    log_counter = log_counter + 1;
    if log_counter == log_update
        log_counter = 0;
        disp(['Timestep = ' num2str(timestep)])
    end
    if timestep == 1
        PFI_FileWrite(filename,timestep,formats,My,gather(Omega),gather(Psi),gather(U),gather(V))
    end
    
    
    %%% COMPUTE FIELD VARIABLES
    % VORTICITY CALCULATION
    Omega(2:My-1,2:Nx-1) = arrayfun(@OmegaCalc,Omega(2:My-1,2:Nx-1),Omega(2:My-1,3:Nx),Omega(2:My-1,1:Nx-2),Omega(3:My,2:Nx-1),Omega(1:My-2,2:Nx-1),U(2:My-1,3:Nx),U(2:My-1,1:Nx-2),V(3:My,2:Nx-1),V(1:My-2,2:Nx-1),Nx,My,alpha,alphaX,alphaY,Cx2,Cy2,d(2:My-1,2:Nx-1),dip1(2:My-1,2:Nx-1),dim1(2:My-1,2:Nx-1),djp1(2:My-1,2:Nx-1),djm1(2:My-1,2:Nx-1));
    % STREAM-FUNCTION CALCULATION
    b = arrayfun(@gen_b,Psi(1:My-2,2:Nx-1),Psi(3:My,2:Nx-1),Psi(2:My-1,1:Nx-2),Psi(2:My-1,3:Nx),Omega(2:My-1,2:Nx-1),xx(2:My-1,2:Nx-1),yy(2:My-1,2:Nx-1),iimat(2:My-1,2:Nx-1),jjmat(2:My-1,2:Nx-1),dx,dy,Nx,My);
    b = reshape(gather(b),(Nx-2)*(My-2),1);
    [Psi(Z),solve_flag,relres,iter(timestep),resvec] = PsiCalc(AMat,gather(Psi(Z)),b,M1,M2,tol);
    % APPLY BOUNDARY CONDITIONS
    [Omega,Psi,U,V] = BoundaryCalc(Omega,Psi,U,V,x,y,dyy,d6,A,IBL,Nx,My);
    % VELOCITY CALCULATION
    [U(2:My-1,2:Nx-1),V(2:My-1,2:Nx-1)] = arrayfun(@VelocityCalc,Psi(3:My,2:Nx-1),Psi(1:My-2,2:Nx-1),Psi(2:My-1,3:Nx),Psi(2:My-1,1:Nx-2),d3(2:My-1,2:Nx-1),dx2,dy2);
    
    
    if output_counter == output_write
        output_counter = 0;
        PFI_FileWrite(filename,timestep,formats,My,gather(Omega),gather(Psi),gather(U),gather(V))
    end
end







