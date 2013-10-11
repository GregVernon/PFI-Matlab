function [U,V]=VelocityCalc(Psijp1,Psijm1,Psiip1,Psiim1,d3,dx2,dy2)
    U =(Psijp1 - Psijm1)./(dy2)./d3;
    V = -(Psiip1- Psiim1)./(dx2)./d3;
end