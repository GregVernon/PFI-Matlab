function [U,V] = VEL_CALC(PSI_jp1,PSI_jm1,PSI_ip1,PSI_im1,di_2dx,di_2dy) %di,dx,dy)

% U = di .* (PSI_jp1 - PSI_jm1) ./ (2*dy);
% V = di .* -(PSI_ip1 - PSI_im1) ./ (2*dx);

U = di_2dy .* (PSI_jp1 - PSI_jm1);
V = di_2dx .* -(PSI_ip1 - PSI_im1);