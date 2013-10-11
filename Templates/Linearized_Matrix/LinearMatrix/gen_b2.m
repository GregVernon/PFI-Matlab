function [b] = gen_b2(Psi,Omega,d,dx,dy,Nx,My)
% (jj,ii)       -->  jj + Nx*(ii-1)
% (jj,ii+1)     -->  jj + Nx*ii
% (jj,ii-1)     -->  jj + Nx*(ii-2)
% (jj+1,ii)     -->  (jj+1) + Nx*(ii-1)
% (jj-1,ii)     -->  (jj-1) + Nx*(ii-1)
b = zeros((Nx-2)*(My-2),1);
kk = 0;
for ii = 2:Nx-1
    for jj = 2:My-1
        kk = kk+1;
        if ii == 2 && jj == 2 % Top Left
            b(kk) = -dx^2*Psi((jj-1) + Nx*(ii-1)) - dy^2*Psi(jj + Nx*(ii-2)) - (dx^2*dy^2)*d(jj + Nx*(ii-1))*Omega(jj + Nx*(ii-1));
        elseif ii == 2 && jj == My-1 % Bottom Left
            b(kk) = -dx^2*Psi((jj+1) + Nx*(ii-1)) - dy^2*Psi(jj + Nx*(ii-2)) - (dx^2*dy^2)*d(jj + Nx*(ii-1))*Omega(jj + Nx*(ii-1));
        elseif ii == Nx-1 && jj == 2 % Top Right
            b(kk) = -dx^2*Psi((jj-1) + Nx*(ii-1)) - dy^2*Psi(jj + Nx*ii) - (dx^2*dy^2)*d(jj + Nx*(ii-1))*Omega(jj + Nx*(ii-1));
        elseif ii == Nx-1 && jj == My-1 % Bottom Right
            b(kk) = -dx^2*Psi((jj+1) + Nx*(ii-1)) - dy^2*Psi(jj + Nx*ii) - (dx^2*dy^2)*d(jj + Nx*(ii-1))*Omega(jj + Nx*(ii-1));
        elseif ii == 2 % Left Side
            b(kk) = -dy^2*Psi(jj + Nx*(ii-2)) - (dx^2*dy^2)*d(jj + Nx*(ii-1))*Omega(jj + Nx*(ii-1));
        elseif ii == Nx-1 % Right Side
            b(kk) = -dy^2*Psi(jj + Nx*ii) - (dx^2*dy^2)*d(jj + Nx*(ii-1))*Omega(jj + Nx*(ii-1));
        elseif jj == 2 % Top Side
            b(kk) = -dx^2*Psi((jj-1) + Nx*(ii-1)) - (dx^2*dy^2)*d(jj + Nx*(ii-1))*Omega(jj + Nx*(ii-1));
        elseif jj == My-1 % Bottom Side
            b(kk) = -dx^2*Psi((jj+1) + Nx*(ii-1)) - (dx^2*dy^2)*d(jj + Nx*(ii-1))*Omega(jj + Nx*(ii-1));
        else % Interior Point w/ no BC
            b(kk) = -(dx^2*dy^2)*d(jj + Nx*(ii-1))*Omega(jj + Nx*(ii-1));
        end
    end
end

end