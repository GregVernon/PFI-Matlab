function [b] = gen_b2(Psi,Omega,d,dx,dy,Nx,My)
% (jj,ii)       -->  jj + Nx*(ii-1)
% (jj,ii+1)     -->  jj + Nx*ii
% (jj,ii-1)     -->  jj + Nx*(ii-2)
% (jj+1,ii)     -->  (jj+1) + Nx*(ii-1)
% (jj-1,ii)     -->  (jj-1) + Nx*(ii-1)
b = zeros((Nx-2)*(My-2),1);
kk = reshape([1:(Nx-2)*(My-2)]',My-2,Nx-2);
TL1 = kk(1,1);
BL1 = kk(My-2,1);
TR1 = kk(1,Nx-2);
BR1 = kk(My-2,Nx-2);
LS1 = kk(2:My-3,1);
RS1 = kk(2:My-3,Nx-2);
TS1 = kk(1,2:Nx-3);
BS1 = kk(My-2,2:Nx-3);
IP1 = kk(2:My-3,2:Nx-3);

MASK = reshape([1:Nx*My]',My,Nx);
TL2 = MASK(2,2);
BL2 = MASK(My-1,2);
TR2 = MASK(2,Nx-1);
BR2 = MASK(My-1,Nx-1);
LS2 = MASK(3:My-2,2);
RS2 = MASK(3:My-2,Nx-1);
TS2 = MASK(2,3:Nx-2);
BS2 = MASK(My-1,3:Nx-2);
IP2 = MASK(3:My-2,3:Nx-2);
%% Top Left -- ii == 2 && jj == 2
b(TL1) = -dx^2*Psi(1,2) - dy^2*Psi(2,1) - (dx^2 * dy^2)*d(TL2).^2*Omega(TL2);
%% Bottom Left -- ii == 2 && jj == My-1
b(BL1) = -dx^2*Psi(My,2) - dy^2*Psi(My-1,1) - (dx^2 * dy^2)*d(BL2).^2*Omega(BL2);
%% Top Right -- ii == Nx-1 && jj == 2
b(TR1) = -dx^2*Psi(1,Nx-1) - dy^2*Psi(2,Nx) - (dx^2 * dy^2)*d(TR2).^2*Omega(TR2);
%% Bottom Right -- ii == Nx-1 && jj == My-1
b(BR1) = -dx^2*Psi(My,Nx-1) - dy^2*Psi(My-1,Nx) -  (dx^2 * dy^2)*d(BR2).^2*Omega(BR2);
%% Left Side -- ii == 2
b(LS1) = -dy^2*Psi(3:My-2,1) -  (dx^2 * dy^2)*d(LS2).^2.*Omega(LS2);
%% Right Side -- ii == Nx-1
b(RS1) = -dy^2*Psi(3:My-2,Nx) - (dx^2 * dy^2)*d(RS2).^2.*Omega(RS2);
%% Top Side -- jj == 2
b(TS1) = -dx^2*Psi(1,3:Nx-2) - (dx^2 * dy^2)*d(TS2).^2.*Omega(TS2);
%% Bottom Side -- jj == My-1
b(BS1) = -dx^2*Psi(My,3:Nx-2) - (dx^2 * dy^2)*d(BS2).^2.*Omega(BS2);
%% Interior Point
b(IP1) = -(dx^2 * dy^2) * d(IP2).^2.*Omega(IP2);


end