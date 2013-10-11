clear
close all
clc
%% Inputs
Re = 700;
A = 1.0;
dt = 0.0007;
Nx = 2000;
My = 2000;
IBL = 100;
xmin = -20;
xmax = 20;
ymin = 1;
ymax = 11;
tol = 1e-5;
end_time = 10;
log_update = 1;
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

iimat = repmat([2:Nx-1],My-2,1);
jjmat = repmat([2:My-1]',1,Nx-2);
% %% LINEARIZE THE VARIABLES
% % !!!!!!!!!! COMMENT HERE !!!!!
% Omega = reshape(Omega,Nx*My,1);
% Psi = reshape(Psi,Nx*My,1);
% U = reshape(U,Nx*My,1);
% V = reshape(V,Nx*My,1);
% 
% % !!!!!!!!!! COMMENT HERE !!!!!
% xx = reshape(xx,Nx*My,1);
% % !!!!!!!!!! COMMENT HERE !!!!!
% d = reshape(d,size(d,1)*size(d,2),1);
% d2 = reshape(d2,size(d2,1)*size(d2,2),1);
% d3 = reshape(d3,size(d3,1)*size(d3,2),1);
% d6 = reshape(d6,size(d6,1)*size(d6,2),1);
% 
% % !!!!!!!!!! COMMENT HERE !!!!!
% dip1 = reshape(dip1,size(dip1,1)*size(dip1,2),1);
% dim1 = reshape(dim1,size(dim1,1)*size(dim1,2),1);
% djp1 = reshape(djp1,size(djp1,1)*size(djp1,2),1);
% djm1 = reshape(djm1,size(djm1,1)*size(djm1,2),1);

%% GENERATE SPARSE FINITE DIFFERENCE MATRIX
[AMat] = FDCM(Nx,My,dx,dy);

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
iter = zeros(end_time,1);
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
    
    Omega0 = Omega;
    Psi0 = Psi;
    
    %%% COMPUTE FIELD VARIABLES
    % VORTICITY CALCULATION
    Omega(2:My-1,2:Nx-1) = OmegaCalc(Omega0,U,V,Nx,My,alpha,alphaX,alphaY,Cx2,Cy2,d,dip1,dim1,djp1,djm1);
    % STREAM-FUNCTION CALCULATION
    b = gen_b(Psi,Omega,d,dx,dy,Nx,My);
    [Psi(Z),solve_flag,relres,iter(timestep),resvec] = PsiCalc(AMat,Psi0(Z),b,M1,M2,tol);
    % APPLY BOUNDARY CONDITIONS
    [Omega,Psi,U,V] = BoundaryCalc(Omega,Psi,U,V,x,y,dyy,d6,A,IBL,Nx,My);
    % VELOCITY CALCULATION
    [U(2:My-1,2:Nx-1),V(2:My-1,2:Nx-1)] = VelocityCalc(Psi,d3,dx2,dy2,Nx,My);
    
end







