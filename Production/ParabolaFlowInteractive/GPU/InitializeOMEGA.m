function OMEGA = InitializeOMEGA(PSI,NX,NY,dyy,d6)

OMEGA = zeros(NY,NX);
for ii=1:NX
    OMEGA(1,ii) = (7.0*PSI(1,ii)-8.0*PSI(2,ii)+PSI(3,ii))/(2.0*dyy)/d6(1,ii);
end