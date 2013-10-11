function [Psi,solve_flag,relres,iter,resvec] = PsiCalc(AMat,Psi0,b,M1,M2,tol)
[Psi,solve_flag,relres,iter,resvec] = tfqmr(AMat,b,tol,1e4,M1,M2,Psi0);
end