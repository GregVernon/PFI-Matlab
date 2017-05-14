function PSI = PSI_CALC(FDM,F,OMEGA,PSI,D2,dx,dy,dxx,KappaA,Kappa2,NX,NY,ii,jj,II,JJ,tol,tstep,method)
if strcmpi(method,'jacobi')
    PsiTol = 1;
    kPsi = 0;
    while PsiTol > tol
        kPsi = kPsi+1;
        PSI0 = PSI;
        PSI(jj,ii) = arrayfun(@PSI_JACOBI,OMEGA(jj,ii),PSI(jj,ii+1),PSI(jj,ii-1),PSI(jj+1,ii),PSI(jj-1,ii),D2(jj,ii),dxx,KappaA,Kappa2);
        PsiTol = max(max(abs(PSI-PSI0)));
    end
elseif strcmpi(method,'direct')
    b = zeros(NY,NX,'gpuArray');
    b(jj,ii) = -arrayfun(@computeRHS,NX,NY,dx,dy,D2(jj,ii),OMEGA(jj,ii),PSI(jj+1,ii),PSI(jj-1,ii),PSI(jj,ii+1),PSI(jj,ii-1),II(jj,ii),JJ(jj,ii),b(jj,ii));
    b = reshape(b(2:NY-1,2:NX-1),(NX-2)*(NY-2),1);
    X = gather(FDM)\gather(b);
    PSI(jj,ii) = gpuArray(reshape(X',NY-2,NX-2));
    %     PSI = PSI_BC(PSI,x,y,NX,NY,IBL,A);
elseif strcmpi(method,'linfactor')
    if tstep == 1
        disp('LINFACTOR Method Chosen'); disp('Performing Initial Matrix Factorization');
        %         FAtic = tic;
        %         F = linfactor(FDM);
        %         toc(FAtic)
    end
    b = zeros(NY,NX,'gpuArray');
    b(jj,ii) = -arrayfun(@computeRHS,NX,NY,dx,dy,D2(jj,ii),OMEGA(jj,ii),PSI(jj+1,ii),PSI(jj-1,ii),PSI(jj,ii+1),PSI(jj,ii-1),II(jj,ii),JJ(jj,ii),b(jj,ii));
    b = reshape(b(2:NY-1,2:NX-1),(NX-2)*(NY-2),1);
    
    X = linfactor(F,gather(b));
    X = gpuArray(X);
    PSI(jj,ii) = reshape(X',NY-2,NX-2);
    %     PSI = PSI_BC(PSI,x,y,NX,NY,IBL,A);
elseif strcmpi(method,'inverseExplicit')
    if tstep == 1
        disp('Inverse Explicit Method Chosen'); disp('Performing Initial Matrix Inversion');
        mItic = tic;
        invA = full(A\speye(size(FDM)));
        toc(mItic)
    end
    b = zeros(NY,NX,'gpuArray');
    b(jj,ii) = -arrayfun(@computeRHS,NX,NY,dx,dy,D2(jj,ii),OMEGA(jj,ii),PSI(jj+1,ii),PSI(jj-1,ii),PSI(jj,ii+1),PSI(jj,ii-1),II(jj,ii),JJ(jj,ii),b(jj,ii));
    b = reshape(b(2:NY-1,2:NX-1),(NX-2)*(NY-2),1);
    X = invA * b;
    PSI(jj,ii) = reshape(X,NY-2,NX-2);
elseif strcmpi(method,'krylov')
    if tstep == 1
        disp('Krylov Method Chosen');
    end
    b = zeros(NY,NX,'gpuArray');
    b(jj,ii) = -arrayfun(@computeRHS,NX,NY,dx,dy,D2(jj,ii),OMEGA(jj,ii),PSI(jj+1,ii),PSI(jj-1,ii),PSI(jj,ii+1),PSI(jj,ii-1),II(jj,ii),JJ(jj,ii),b(jj,ii));
    b = reshape(b(2:NY-1,2:NX-1),(NX-2)*(NY-2),1);
    X0 = reshape(PSI(jj,ii),(NY-2)*(NX-2),1);
    
%    [X,flag,relres,iter] = pcg(FDM,b,tol,1e5,F.M1,F.M2,X0);
%    [X,flag,relres,iter] = gmres(FDM,b,100,tol,1e5,[],[],X0);
%     X = G_PCG(FDM,b,X0,tol);
    X = conjgrad(FDM,b,X0,tol);    
    PSI(jj,ii) = reshape(X',NY-2,NX-2);
end