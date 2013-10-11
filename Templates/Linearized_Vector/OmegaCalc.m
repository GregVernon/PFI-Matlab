function Omega = OmegaCalc(Omega0,Omega,U,V,Nx,My,alpha,alphaX,alphaY,Cx2,Cy2,d,dip1,dim1,djp1,djm1)

% ii=2:Nx-1;
% jj=2:My-1;

% (jj,ii)       -->  jj + Nx*(ii-1)
% (jj,ii+1)     -->  jj + Nx*ii
% (jj,ii-1)     -->  jj + Nx*(ii-2)
% (jj+1,ii)     -->  (jj+1) + Nx*(ii-1)
% (jj-1,ii)     -->  (jj-1) + Nx*(ii-1)
for ii = 2:Nx-1;
    for jj = 2:My-1;
        
        Omega(jj + Nx*(ii-1)) = Omega0(jj + Nx*(ii-1)).*(1-alpha./(d(jj + Nx*(ii-1)).^2)) + ...
            Omega0(jj + Nx*ii).*(-Cx2.*U(jj + Nx*ii).*dip1(jj + Nx*(ii-1)) + alphaX)./(d(jj + Nx*(ii-1)).^2) + ...
            Omega0(jj + Nx*(ii-2)).*( Cx2.*U(jj + Nx*(ii-2)).*dim1(jj + Nx*(ii-1)) + alphaX)./(d(jj + Nx*(ii-1)).^2) + ...
            Omega0((jj+1) + Nx*(ii-1)).*(-Cy2.*V((jj+1) + Nx*(ii-1)).*djp1(jj + Nx*(ii-1)) + alphaY)./(d(jj + Nx*(ii-1)).^2) + ...
            Omega0((jj-1) + Nx*(ii-1)).*( Cy2.*V((jj-1) + Nx*(ii-1)).*djm1(jj + Nx*(ii-1)) + alphaY)./(d(jj + Nx*(ii-1)).^2);
        
        
    end
end
end