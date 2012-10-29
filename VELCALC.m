function [u,v]=VELCALC(Psiip1,Psiim1,Psijp1,Psijm1,d3,dx2,dy2)
    u =(Psiip1-Psiim1)./(dy2)./d3;
    v = -(Psijp1-Psijm1)./(dx2)./d3;
end