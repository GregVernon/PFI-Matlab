function x = G_PCG(A,b,x0,Tol)
tStart = tic;
r = b-A*x0;
s = r;
x = x0;
kPsi = 0;
At = A';
while norm(r) > Tol
    kPsi = kPsi+1;
    a = r'*r/((At*s)'*s);
    x = x + a*s;
    
    r1 = r - (At*a)'*s;
    beta = r1'*r1/(r'*r);
    s = r1 + beta*s;
    r = r1;
end
tEnd = toc(tStart);