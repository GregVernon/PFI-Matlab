function b = computeRHS_2(NX,NY,dx,dy,D2,OMEGA,PSI,b)
%% Generate b for a domain from 2:N-1
%% Bottom Left Node
% b(jj,ii) = (dx^2 * D2(jj,ii) * OMEGA(jj,ii) + PSI(jj,ii-1) + (dx^2/dy^2)*PSI(jj-1,ii));
% if jj == 2 && ii == 2
jj = 2;
ii = 2;
b(jj,ii) = -((dx^2 * dy^2) * D2(jj,ii) .* OMEGA(jj,ii)) + (dy^2*PSI(jj,ii-1) + dx^2*PSI(jj-1,ii));
%% Bottom Right Node
% elseif jj == 2 && ii == NX-1
jj = 2;
ii = NX-1;
b(jj,ii) = -((dx^2 * dy^2) * D2(jj,ii) .* OMEGA(jj,ii)) + (dy^2*PSI(jj,ii+1) - dx^2*PSI(jj-1,ii));
%% Top Left Node
% elseif jj == NY-1 && ii == 2
jj = NY-1;
ii = 2;
b(jj,ii) = -((dx^2 * dy^2) * D2(jj,ii) .* OMEGA(jj,ii)) - (dx^2*PSI(jj+1,ii) + dy^2*PSI(jj,ii-1));
%% Top Right Node
% elseif jj == NY-1 && ii == NX-1
jj = NY-1;
ii = NX-1;
b(jj,ii) = -((dx^2 * dy^2) * D2(jj,ii) .* OMEGA(jj,ii)) - (dx^2*PSI(jj+1,ii) + dy^2*PSI(jj,ii+1));
%% Left Node
% elseif jj > 2 && jj < NY-1 && ii == 2
jj = 3:NY-2;
ii = 2;
b(jj,ii) = -((dx^2 * dy^2) * D2(jj,ii) .* OMEGA(jj,ii)) - (dy^2*PSI(jj,ii-1));
%% Right Node
% elseif jj > 2 && jj < NY-1 && ii == NX-1
jj = 3:NY-2;
ii = NX-1;
b(jj,ii) = -((dx^2 * dy^2) * D2(jj,ii) .* OMEGA(jj,ii)) - (dy^2*PSI(jj,ii+1));
%% Bottom Node
% elseif jj == 2 && ii > 2 && ii < NX-1
jj = 2;
ii = 3:NX-2;
b(jj,ii) = -((dx^2 * dy^2) * D2(jj,ii) .* OMEGA(jj,ii)) - (dx^2*PSI(jj-1,ii));
%% Top Node
%elseif jj == NY-1 && ii > 2 && ii < NX-1
jj = NY-1;
ii = 3:NX-2;
b(jj,ii) = -((dx^2 * dy^2) * D2(jj,ii) .* OMEGA(jj,ii)) - (dx^2*PSI(jj+1,ii));

%% Interior Node
% elseif jj >= 3 && jj <= NY-2 && ii >= 3 && ii <= NX-2
jj = 3:NY-2;
ii = 3:NX-2;
b(jj,ii) = -((dx^2 * dy^2) * D2(jj,ii) .* OMEGA(jj,ii));

end
