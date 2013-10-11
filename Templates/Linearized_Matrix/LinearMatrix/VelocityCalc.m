function [U,V]=VelocityCalc(Psi,d3,dx2,dy2,Nx,My)
    U =(Psi(3:My,2:Nx-1) - Psi(1:My-2,2:Nx-1))./(dy2)./d3(2:My-1,2:Nx-1);
    V = -(Psi(2:My-1,3:Nx) - Psi(2:My-1,1:Nx-2))./(dx2)./d3(2:My-1,2:Nx-1);
end