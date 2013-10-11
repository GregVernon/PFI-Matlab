function [b] = gen_b(Psijm1,Psijp1,Psiim1,Psiip1,Omega,xx,yy,ii,jj,dx,dy,Nx,My)
    
    if ii == 2 && jj == 2 % Top Left 
        b = -dx.^2.*Psijm1 - dy.^2.*Psiim1 - (dx.^2.*dy.^2).*(xx.^2 + yy.^2).*Omega;
    elseif ii == 2 && jj == My-1 % Bottom Left
        b = -dx.^2.*Psijp1 - dy.^2.*Psiim1 - (dx.^2.*dy.^2).*(xx.^2 + yy.^2).*Omega;
    elseif ii == Nx-1 && jj == 2 % Top Right
        b = -dx.^2.*Psijm1 - dy.^2.*Psiip1 - (dx.^2.*dy.^2).*(xx.^2 + yy.^2).*Omega;
    elseif ii == Nx-1 && jj == My-1 % Bottom Right
        b = -dx.^2.*Psijp1 - dy.^2.*Psiip1 - (dx.^2.*dy.^2).*(xx.^2 + yy.^2).*Omega;
    elseif ii == 2 % Left Side
        b = -dy.^2.*Psiim1 - (dx.^2.*dy.^2).*(xx.^2 + yy.^2).*Omega;
    elseif ii == Nx-1 % Right Side
        b = -dy.^2.*Psiip1 - (dx.^2.*dy.^2).*(xx.^2 + yy.^2).*Omega;
    elseif jj == 2 % Top Side
        b = -dx.^2.*Psijm1 - (dx.^2.*dy.^2).*(xx.^2 + yy.^2).*Omega;
    elseif jj == My-1 % Bottom Side
        b = -dx.^2.*Psijp1 - (dx.^2.*dy.^2).*(xx.^2 + yy.^2).*Omega;
    else % Interior Point w/ no BC
        b = -(dx.^2.*dy.^2).*(xx.^2 + yy.^2).*Omega;
    end


end