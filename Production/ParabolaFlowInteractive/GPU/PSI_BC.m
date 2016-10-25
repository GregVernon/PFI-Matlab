function PSI = PSI_BC(PSI,x,y,NX,NY,IBL,A)

% Lower
jj = 1;
PSI(jj,1:NX) = 0;

% Upper
jj = NY;
PSI(jj,1:NX) = (x+A)*(y(jj)-1);

% Left
ii = 1;
PSI(2:IBL,ii) = PSI(2:IBL,ii+1) ;
PSI(IBL+1:NY-1,ii) = (x(ii)+A)*(y(IBL+1:NY-1)-1);

% Right
ii = NX;
PSI(2:IBL,ii) = PSI(2:IBL,ii-1) ;
PSI(IBL+1:NY-1,ii) = (x(ii)+A)*(y(IBL+1:NY-1)-1);