function [U,V] = VELOCITY_BC(U,V,x,y,NX,NY,IBL,A)

% Lower
jj = 1;
U(jj,1:NX) = 0;
V(jj,1:NX) = 0;

% Upper
jj = NY;
U(jj,1:NX) = (x+A)./sqrt(x.^2+y(jj)^2);
V(jj,1:NX) =-(y(jj)-1)./sqrt(x.^2+y(jj)^2);

% Left
ii = 1;
U(2:IBL,ii) = U(2:IBL,ii+1);
V(2:IBL,ii) = V(2:IBL,ii+1);
U(IBL+1:NY-1,ii) = (x(ii)+A)./sqrt(x(ii)^2+y(IBL+1:NY-1).^2);
V(IBL+1:NY-1,ii) =-(y(IBL+1:NY-1)-1)./sqrt(x(ii)^2+y(IBL+1:NY-1).^2);

% Right
ii = NX;
U(2:IBL,ii) = U(2:IBL,ii-1);
V(2:IBL,ii) = V(2:IBL,ii-1);
U(IBL+1:NY-1,ii) = (x(ii)+A)./sqrt(x(ii)^2+y(IBL+1:NY-1).^2);
V(IBL+1:NY-1,ii) =-(y(IBL+1:NY-1)-1)./sqrt(x(ii)^2+y(IBL+1:NY-1).^2);