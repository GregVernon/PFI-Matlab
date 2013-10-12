function [Omega,Psi,U,V] = BoundaryCalc(Omega,Psi,U,V,Nx,My,A,IBL,x,y,dyy,d,d6)
% (jj,ii)       -->  jj + Nx*(ii-1)
% (jj,ii+1)     -->  jj + Nx*ii
% (jj,ii-1)     -->  jj + Nx*(ii-2)
% (jj+1,ii)     -->  (jj+1) + Nx*(ii-1)
% (jj-1,ii)     -->  (jj-1) + Nx*(ii-1)
%% NO SLIP BOUNDARY CONDITION
jj = 1;
for ii = 1:Nx
    Psi(jj + Nx*(ii-1)) = 0;
    Omega(jj + Nx*(ii-1)) = (7.0*Psi(jj + Nx*(ii-1))-8.0*Psi(jj+1 + Nx*(ii-1))+Psi(jj+2 + Nx*(ii-1)))/(2.0*dyy)/d6(ii);
    U(jj + Nx*(ii-1)) = 0;
    V(jj + Nx*(ii-1)) = 0;
end

%% UPPER BOUNDARY CONDITION
jj=My;
for ii = 1:Nx
    Omega(jj + Nx*(ii-1)) = 0;
    Psi(jj + Nx*(ii-1)) = (x(ii)+A)*(y(jj)-1);
    U(jj + Nx*(ii-1)) = (x(ii)+A)/d(jj + Nx*(ii-1));
    V(jj + Nx*(ii-1)) = -(y(jj)-1)/d(jj + Nx*(ii-1));
end

%% SIDE BOUNDARY CONDITIONS
for jj = 2:My
    if jj <= IBL
        %%%% INSIDE BOUNDARY LAYER
        ii=1;
        Omega(jj + Nx*(ii-1))= Omega(jj + Nx*ii);
        Psi(jj + Nx*(ii-1)) = Psi(jj + Nx*ii) ;
        U(jj + Nx*(ii-1)) = U(jj + Nx*ii);
        V(jj + Nx*(ii-1)) = V(jj + Nx*ii);
        ii=Nx;
        Omega(jj + Nx*(ii-1)) = Omega(jj + Nx*(ii-2));
        Psi(jj + Nx*(ii-1)) = Psi(jj + Nx*(ii-2));
        U(jj + Nx*(ii-1)) = U(jj + Nx*(ii-2));
        V(jj + Nx*(ii-1)) = V(jj + Nx*(ii-2));
    else
        %%%% OUTSIDE BOUNDARY LAYER
        ii=1;
        Omega(jj + Nx*(ii-1)) = 0.0;
        Psi(jj + Nx*(ii-1)) = (x(ii)+A)*(y(jj)-1);
        U(jj + Nx*(ii-1)) = (x(ii)+A)/d(jj + Nx*(ii-1));
        V(jj + Nx*(ii-1)) = -(y(jj)-1)/d(jj + Nx*(ii-1));
        ii=Nx;
        Omega(jj + Nx*(ii-1)) = 0.0;
        Psi(jj + Nx*(ii-1)) = (x(ii)+A)*(y(jj)-1);
        U(jj + Nx*(ii-1)) = (x(ii)+A)/d(jj + Nx*(ii-1));
        V(jj + Nx*(ii-1)) = -(y(jj)-1)/d(jj + Nx*(ii-1));
    end
end
end