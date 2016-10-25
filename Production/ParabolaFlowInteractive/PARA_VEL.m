function [pU,pV,VEL] = PARA_VEL(U,V,cCOORD,pCOORD)
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
% pU=zeros(NY,NX);
% pV=zeros(NY,NX);

% pU(2:NY,1:NX)=(U(2:NY,1:NX).*dety(2:NY,1:NX)-V(2:NY,1:NX).*dmuy(2:NY,1:NX)).*sqrt(mu(2:NY,1:NX).^2+eta(2:NY,1:NX).^2);  %Vx
% pV(2:NY,1:NX)=-(U(2:NY,1:NX).*detx(2:NY,1:NX)-V(2:NY,1:NX).*dmux(2:NY,1:NX)).*sqrt(mu(2:NY,1:NX).^2+eta(2:NY,1:NX).^2); %Vy

pU = (U.*dety - V.*dmuy) .* sqrt(mu.^2 + eta.^2);  %Vx
pV = -(U.*detx - V.*dmux) .* sqrt(mu.^2 + eta.^2); %Vy
pU(1,:) = 0;
pV(1,:) = 0;

%   Creates V-vctr from Vx and Vy
VEL=zeros(size(pU,1),NX);
VEL(1:size(pU,1),1:round(NX/2))=sqrt(pU(1:size(pU,1),1:round(NX/2)).^2+pV(1:size(pU,1),1:round(NX/2)).^2);
VEL(1:size(pU,1),round(NX/2):NX)=sqrt(pU(1:size(pU,1),round(NX/2):NX).^2+pV(1:size(pU,1),round(NX/2):NX).^2).*sign(pU(1:size(pU,1),round(NX/2):NX));

