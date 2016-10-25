function [U,V] = InitializeVELOCITY(A,x,y,NX,NY)

U = zeros(NY,NX);
for ii=1:NX
    for jj=2:NY
        U(jj,ii)=(x(ii)+A)/sqrt(x(ii)^2+y(jj)^2);
    end
end

V = zeros(NY,NX);
for ii=1:NX
    for jj=2:NY
        V(jj,ii)=-(y(jj)-1)/sqrt(x(ii)^2+y(jj)^2);
    end
end
