function Omega2=OMEGACALC(Omega0,Omega0ip1,Omega0im1,Omega0jp1,Omega0jm1...
    ,uip1,uim1,vjp1,vjm1,Cx2,Cy2,d,dip1,dim1,djp1,djm1,alpha,alphaX,alphaY)

    Omega2 = Omega0.*(1-alpha./(d.*d)) + ...
    Omega0ip1.*(-Cx2.*uip1.*dip1 + alphaX)./(d.*d) + ...
    Omega0im1.*( Cx2.*uim1.*dim1 + alphaX)./(d.*d) + ...
    Omega0jp1.*(-Cy2.*vjp1.*djp1 + alphaY)./(d.*d) + ...
    Omega0jm1.*( Cy2.*vjm1.*djm1 + alphaY)./(d.*d);
end