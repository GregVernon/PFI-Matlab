function [U,V] = VelocityCalc(Psi,U,V,dx2,dy2,x,y,Nx,My)
% (jj,ii)       -->  jj + Nx*(ii-1)
% (jj,ii+1)     -->  jj + Nx*ii
% (jj,ii-1)     -->  jj + Nx*(ii-2)
% (jj+1,ii)     -->  (jj+1) + Nx*(ii-1)
% (jj-1,ii)     -->  (jj-1) + Nx*(ii-1)

for ii = 2:Nx-1
    for jj = 2:My-1
        U(jj + Nx*(ii-1)) = (Psi((jj+1) + Nx*(ii-1)) - Psi((jj-1) + Nx*(ii-1))) / (dy2*sqrt(x(ii)^2 + y(jj)^2));
        V(jj + Nx*(ii-1)) = -(Psi(jj + Nx*ii) - Psi(jj + Nx*(ii-2))) / (dx2*sqrt(x(ii)^2 + y(jj)^2));
    end
end
end