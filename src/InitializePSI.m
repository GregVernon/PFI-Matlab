function PSI = InitializePSI(A,x,y,NX,NY)

PSI = zeros(NY,NX);
for ii=1:NX
    for jj=1:NY
        PSI(jj,ii)=(x(ii)+A)*(y(jj)-1);
    end
end