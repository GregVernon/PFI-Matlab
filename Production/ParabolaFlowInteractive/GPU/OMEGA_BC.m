function OMEGA = OMEGA_BC(PSI,OMEGA,NX,NY,IBL,dy,d6)

% Lower
jj = 1;
OMEGA(jj,1:NX)=(7.0*PSI(1,1:NX)-8.0*PSI(2,1:NX)+PSI(3,1:NX))...
                    /(2.0*dy^2)./d6;
% Upper
jj = NY;
OMEGA(jj,1:NX) = 0;

% Left
ii = 1;
OMEGA(2:IBL,ii)= OMEGA(2:IBL,ii+1);
OMEGA(IBL+1:NY-1,ii) = 0.0;

% Right
ii = NX;
OMEGA(2:IBL,ii) = OMEGA(2:IBL,ii-1);
OMEGA(IBL+1:NY-1,ii) = 0.0;
