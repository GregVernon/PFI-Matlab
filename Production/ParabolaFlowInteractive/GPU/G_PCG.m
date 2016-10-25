function x = G_PCG(A,b,x0,Tol)
tStart = tic;
r = b-A*x0;
s = r;
x = x0;
kPsi = 0;
while norm(r) > Tol
    kPsi = kPsi+1;
%     a = r'*r/(s'*A*s);
    a = r'*r/((A'*s)'*s);
    x = x + a*s;
%     r1 = r - a*A*s;
    r1 = r - (A'*a)'*s;
    beta = r1'*r1/(r'*r);
    s = r1 + beta*s;
    r = r1;
    %         clc
    % %         partPsi = reshape(x,My,Nx);
    %         partresid = reshape(r,My,Nx);
    %         figure(gcf)
    % %         surf(partPsi,'LineStyle','none')
    %         imagesc(partresid)
    %         drawnow
    %         loglog(kPsi,norm(r),'k')
    %         drawnow
    %         hold on
    %         disp(num2str(kPsi))
    %         disp(num2str(norm(r)))
end
tEnd = toc(tStart);