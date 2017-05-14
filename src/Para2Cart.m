function cCOORD = Para2Cart(pCOORD)
NX = pCOORD.NX;
NY = pCOORD.NY;
xMin = pCOORD.xMin;
xMax = pCOORD.xMax;
yMin = pCOORD.yMin;
yMax = pCOORD.yMax;

xskew = 1;
yskew = 1;

pX=zeros(1,NX);
pY=zeros(1,NY);
for ii=1:NX
    if (ii<=abs(xMin/(xMax-xMin))*(NX-1))
        pX(ii)=(-1/2)*abs(xMax-xMin)*abs((ii-((NX-1)/2)-1)/((NX-1)/2))^xskew;
    else
        pX(ii)=(1/2)*abs(xMax-xMin)*abs((ii-((NX-1)/2)-1)/((NX-1)/2))^xskew;
    end
end

for ii=1:NY
    pY(ii)=yMin+abs(yMax-yMin)*((ii-1)/(NY-1))^yskew;
end


[mu,eta]=meshgrid(pX,pY);
pX=1/2*(mu.^2-eta.^2);
pY=mu.*eta;
bx=1/2*(mu(1,:).^2-1);
by=mu(1,:);

r=sqrt(pX.^2+pY.^2);
dmux=zeros(NY,NX);
dmuy=zeros(NY,NX);
detx=zeros(NY,NX);
dety=zeros(NY,NX);

dmux(2:NY,1:round((NX/2))-1)=1/2.*(pX(2:NY,1:round((NX/2))-1)./r(2:NY,1:round((NX/2))-1)+1)./(mu(2:NY,1:round((NX/2))-1));
dmux(2:NY,round((NX/2))+1:NX)=1/2.*(pX(2:NY,round((NX/2))+1:NX)./r(2:NY,round((NX/2))+1:NX)+1)./(mu(2:NY,round((NX/2))+1:NX));
dmuy(2:NY,1:round((NX/2))-1)=1/2.*(pY(2:NY,1:round((NX/2))-1)./r(2:NY,1:round((NX/2))-1))./mu(2:NY,1:round((NX/2))-1) ;
dmuy(2:NY,round((NX/2))+1:NX)=1/2.*(pY(2:NY,round((NX/2))+1:NX)./r(2:NY,round((NX/2))+1:NX))./mu(2:NY,round((NX/2))+1:NX) ;
detx(2:NY,1:round((NX/2))-1)=1/2.*(pX(2:NY,1:round((NX/2))-1)./r(2:NY,1:round((NX/2))-1)-1)./(eta(2:NY,1:round((NX/2))-1));
detx(2:NY,round((NX/2))+1:NX)=1/2.*(pX(2:NY,round((NX/2))+1:NX)./r(2:NY,round((NX/2))+1:NX)-1)./(eta(2:NY,round((NX/2))+1:NX));
dety(2:NY,1:round((NX/2))-1)=1/2.*(pY(2:NY,1:round((NX/2))-1)./r(2:NY,1:round((NX/2))-1))./eta(2:NY,1:round((NX/2))-1) ;
dety(2:NY,round((NX/2))+1:NX)=1/2.*(pY(2:NY,round((NX/2))+1:NX)./r(2:NY,round((NX/2))+1:NX))./eta(2:NY,round((NX/2))+1:NX) ;

% average mu=0 line:
dmux(1:NY,round(NX/2))=(dmux(1:NY,round((NX/2))-1)+dmux(1:NY,round((NX/2))+1))/2;
dmuy(1:NY,round(NX/2))=(dmuy(1:NY,round((NX/2))-1)+dmuy(1:NY,round((NX/2))+1))/2;
detx(1:NY,round(NX/2))=(detx(1:NY,round((NX/2))-1)+detx(1:NY,round((NX/2))+1))/2;
dety(1:NY,round(NX/2))=(dety(1:NY,round((NX/2))-1)+dety(1:NY,round((NX/2))+1))/2;

% Package Output
cCOORD.bx = bx;
cCOORD.by = by;
cCOORD.pX = pX;
cCOORD.pY = pY;
cCOORD.mu = gpuArray(mu);
cCOORD.eta = gpuArray(eta);
cCOORD.dmux = gpuArray(dmux);
cCOORD.dmuy = gpuArray(dmuy);
cCOORD.detx = gpuArray(detx);
cCOORD.dety = gpuArray(dety);
