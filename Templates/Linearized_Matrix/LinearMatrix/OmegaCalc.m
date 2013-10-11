function Omega2=OmegaCalc(Omega0,U,V,Nx,My,alpha,alphaX,alphaY,Cx2,Cy2,d,dip1,dim1,djp1,djm1)

    Omega2 = Omega0(2:My-1,2:Nx-1).*(1-alpha./(d(2:My-1,2:Nx-1).^2)) + ...
    Omega0(2:My-1,3:Nx).*(-Cx2.*U(2:My-1,3:Nx).*dip1(2:My-1,2:Nx-1) + alphaX)./(d(2:My-1,2:Nx-1).^2) + ...
    Omega0(2:My-1,1:Nx-2).*( Cx2.*U(2:My-1,1:Nx-2).*dim1(2:My-1,2:Nx-1) + alphaX)./(d(2:My-1,2:Nx-1).^2) + ...
    Omega0(3:My,2:Nx-1).*(-Cy2.*V(3:My,2:Nx-1).*djp1(2:My-1,2:Nx-1) + alphaY)./(d(2:My-1,2:Nx-1).^2) + ...
    Omega0(1:My-2,2:Nx-1).*( Cy2.*V(1:My-2,2:Nx-1).*djm1(2:My-1,2:Nx-1) + alphaY)./(d(2:My-1,2:Nx-1).^2);
end
