function [Psi2]=PSI_JACOBI(Omega,Psi0ip1,Psi0im1,Psi0jp1,Psi0jm1,d2,dxx,KappaA,Kappa2)

Psi2 = KappaA.*(dxx.*Omega.*d2 + Psi0ip1 + ...
    Psi0im1 + Kappa2.*(Psi0jp1 + Psi0jm1));

end