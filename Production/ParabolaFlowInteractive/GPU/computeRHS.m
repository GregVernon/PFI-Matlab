function b = computeRHS(NX,NY,dx,dy,D2,OMEGA,PSI_jp1,PSI_jm1,PSI_ip1,PSI_im1,ii,jj,b)
%% Generate b for a domain from 2:N-1
%% Bottom Left Node
if jj == 2 && ii == 2
    % b(jj,ii) = (dx^2 * D2(jj,ii) * OMEGA(jj,ii) + PSI(jj,ii-1) + (dx^2/dy^2)*PSI(jj-1,ii));
    b = -((dx^2 * dy^2) * D2 .* OMEGA) + (dy^2*PSI_im1 + dx^2*PSI_jm1);
elseif jj == 2 && ii == NX-1;
    % Bottom Right Node
    b = -((dx^2 * dy^2) * D2 .* OMEGA) + (dy^2*PSI_ip1 - dx^2*PSI_jm1);
elseif jj == NY-1 && ii == 2;
    % Top Left Node
    b = -((dx^2 * dy^2) * D2 .* OMEGA) - (dx^2*PSI_jp1 + dy^2*PSI_im1);
elseif jj == NY-1 && ii == NX-1
    % Top Right Node
    b = -((dx^2 * dy^2) * D2 .* OMEGA) - (dx^2*PSI_jp1 + dy^2*PSI_ip1);
elseif jj > 2 && jj < NY-1 && ii == 2
    % Left Node
    b = -((dx^2 * dy^2) * D2 .* OMEGA) - (dy^2*PSI_im1);
elseif jj > 2 && jj < NY-1 && ii == NX-1
    % Right Node
    b = -((dx^2 * dy^2) * D2 .* OMEGA) - (dy^2*PSI_ip1);
elseif jj == 2 && ii > 2 && ii < NX-1
    % Bottom Node
    b = -((dx^2 * dy^2) * D2 .* OMEGA) - (dx^2*PSI_jm1);
elseif jj == NY-1 && ii > 2 && ii < NX-1
    % Top Node
    b = -((dx^2 * dy^2) * D2 .* OMEGA) - (dx^2*PSI_jp1);
elseif jj >= 3 && jj <= NY-2 && ii >= 3 && ii <= NX-2;
    % Interior Node
    b = -((dx^2 * dy^2) * D2 .* OMEGA);
end
% Reshape into a vector
% b = reshape(b(2:NY-1,2:NX-1),(NX-2)*(NY-2),1);