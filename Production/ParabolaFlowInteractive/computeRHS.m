function b = computeRHS(NX,NY,x,y,dx,dy,D2,OMEGA,PSI,A,IBL)
%% Generate b for a domain from 1:N
%% Unvectorized
% b = zeros(N^2,1);
% node = 0;
% for jj = 1:N
%     for ii = 1:N
%         node = node+1;
%         if ii == 1 || ii == N || jj == 1 || jj == N
%             % External Node
%             b(node) = Omega(jj,ii);
%         elseif ii == 2 && jj == 2
%             % Bottom Left Node
%             b(node) = -(dx^2 * D2(jj,ii) * Omega(jj,ii) - Psi(jj,ii-1) - Psi(jj-1,ii))/4;
%         elseif ii == N-1 && jj == 2
%             % Bottom Right Node
%             b(node) = -(dx^2 * D2(jj,ii) * Omega(jj,ii) - Psi(jj,ii+1) - Psi(jj-1,ii))/4;
%         elseif ii == 2 && jj == N-1
%             % Top Left Node
%             b(node) = -(dx^2 * D2(jj,ii) * Omega(jj,ii) - Psi(jj,ii-1) - Psi(jj+1,ii))/4;
%         elseif ii == N-1 && jj == N-1
%             % Top Right Node
%             b(node) = -(dx^2 * D2(jj,ii) * Omega(jj,ii) - Psi(jj,ii+1) - Psi(jj+1,ii))/4;
%         elseif ii == 2
%             % Left Node
%             b(node) = -(dx^2 * D2(jj,ii) * Omega(jj,ii) - Psi(jj,ii-1))/4;
%         elseif ii == N-1
%             % Right Node
%             b(node) = -(dx^2 * D2(jj,ii) * Omega(jj,ii) - Psi(jj,ii+1))/4;
%         elseif jj == 2
%             % Bottom Node
%             b(node) = -(dx^2 * D2(jj,ii) * Omega(jj,ii) - Psi(jj-1,ii))/4;
%         elseif jj == N-1
%             % Top Node
%             b(node) = -(dx^2 * D2(jj,ii) * Omega(jj,ii) - Psi(jj+1,ii))/4;
%         else
%             % Interior Node
%             b(node) = -(dx^2 * D2(jj,ii) * Omega(jj,ii))/4;
%         end
%     end
% end

%% Vectorized
b = zeros(NY,NX);
%% External Nodes
% % Left
% ii = 1;
% b(2:IBL,ii) = 0; %OMEGA(:,ii);
% b(IBL+1:NY-1,ii) = (x(ii)+A)*(y(IBL+1:NY-1)-1);
% 
% % Right
% ii = NX;
% b(2:IBL,ii) = 0; %OMEGA(:,ii);
% b(IBL+1:NY-1,ii) = (x(ii)+A)*(y(IBL+1:NY-1)-1);
% 
% % Bottom
% jj = 1;
% b(jj,:) = 0; %OMEGA(jj,:);
% 
% % Top
% jj = NY;
% b(jj,:) = (x+A) * (y(jj) - 1); %OMEGA(jj,:);

%% Bottom Left Node
jj = 2;
ii = 2;
% b(jj,ii) = (dx^2 * D2(jj,ii) * OMEGA(jj,ii) + PSI(jj,ii-1) + (dx^2/dy^2)*PSI(jj-1,ii));
b(jj,ii) = -((dx^2 * dy^2) * D2(jj,ii) .* OMEGA(jj,ii)) - (dy^2*PSI(jj,ii-1)) - (dx^2*PSI(jj-1,ii));
% Bottom Right Node
jj = 2;
ii = NX-1;
% b(jj,ii) = (dx^2 * D2(jj,ii) * OMEGA(jj,ii) + PSI(jj,ii+1) + (dx^2/dy^2)*PSI(jj-1,ii));
b(jj,ii) = -((dx^2 * dy^2) * D2(jj,ii) .* OMEGA(jj,ii)) - (dy^2*PSI(jj,ii+1) - (dx^2*PSI(jj-1,ii)));
% Top Left Node
jj = NY-1;
ii = 2;
% b(jj,ii) = (dx^2 * D2(jj,ii) * OMEGA(jj,ii) + PSI(jj,ii-1) + (dx^2/dy^2)*PSI(jj+1,ii));
b(jj,ii) = -((dx^2 * dy^2) * D2(jj,ii) .* OMEGA(jj,ii)) - (dy^2*PSI(jj,ii-1) - (dx^2*PSI(jj+1,ii)));
% Top Right Node
jj = NY-1;
ii = NX-1;
% b(jj,ii) = (dx^2 * D2(jj,ii) * OMEGA(jj,ii) + PSI(jj,ii+1) + (dx^2/dy^2)*PSI(jj+1,ii));
b(jj,ii) = -((dx^2 * dy^2) * D2(jj,ii) .* OMEGA(jj,ii)) - (dy^2*PSI(jj,ii+1) - (dx^2*PSI(jj+1,ii)));
% Left Node
jj = 2:NY-1;
ii = 2;
% b(jj,ii) = (dx^2 * D2(jj,ii) .* OMEGA(jj,ii) + PSI(jj,ii-1));
b(jj,ii) = -((dx^2 * dy^2) * D2(jj,ii) .* OMEGA(jj,ii)) - (dy^2*PSI(jj,ii-1));
% Right Node
jj = 2:NY-1;
ii = NX-1;
% b(jj,ii) = (dx^2 * D2(jj,ii) .* OMEGA(jj,ii) + PSI(jj,ii+1));
b(jj,ii) = -((dx^2 * dy^2) * D2(jj,ii) .* OMEGA(jj,ii)) - (dy^2*PSI(jj,ii+1));
% Bottom Node
jj = 2;
ii = 2:NX-1;
% b(jj,ii) = (dx^2 * D2(jj,ii) .* OMEGA(jj,ii) + (dx^2/dy^2)*PSI(jj-1,ii));
b(jj,ii) = -((dx^2 * dy^2) * D2(jj,ii) .* OMEGA(jj,ii)) - (dx^2*PSI(jj-1,ii));
% Top Node
jj = NY-1;
ii = 2:NX-1;
% b(jj,ii) = (dx^2 * D2(jj,ii) .* OMEGA(jj,ii) + (dx^2/dy^2)*PSI(jj+1,ii));
b(jj,ii) = -((dx^2 * dy^2) * D2(jj,ii) .* OMEGA(jj,ii)) - (dx^2*PSI(jj+1,ii));
% Interior Node
jj = 3:NY-2;
ii = 3:NX-2;
% b(jj,ii) = (dx^2 * D2(jj,ii) .* OMEGA(jj,ii));
b(jj,ii) = -((dx^2 * dy^2) * D2(jj,ii) .* OMEGA(jj,ii));

% Reshape into a vector
b = reshape(b(2:NY-1,2:NX-1),(NX-2)*(NY-2),1);