function A = assembleCoeffMat(NX, NY, dx, dy)
%% Generate A for a domain from 1:N
StencilTemplate = spalloc(NY-2,NX-2,5);
sRow = cell(NY-2,NX-2);
sCol = cell(NY-2,NX-2);
sVal = cell(NY-2,NX-2);
% syms C L R B T
C = -2*(dy^2 + dx^2); %-4;
L = dy^2; %1;
R = dy^2; %1;
B = dx^2; %1;
T = dx^2; %1;
nNodes = (NX-2)*(NY-2);
parfor node = 1:nNodes
    Stencil = StencilTemplate;
    [jj,ii] = ind2sub([NY-2,NX-2],node);
    if ii == 1 && jj == 1
        % Bottom Left Node
        % Slower
        Stencil(jj,ii) = C;     %-4;
        Stencil(jj,ii+1) = R;   %1;
        Stencil(jj+1,ii) = T;   %1;
        
        % Faster
        %             srow = [jj jj jj+1];
        %             scol = [ii ii+1 ii];
        %             sval = [C R T];
        %             Stencil = sparse(srow,scol,sval);
    elseif ii == NX-2 && jj == 1
        % Bottom Right Node
        % Slower
        Stencil(jj,ii) = C;     %-4;
        Stencil(jj,ii-1) = L;   %1;
        Stencil(jj+1,ii) = T;   %1;
        
        % Faster
        %             srow = [jj jj jj+1];
        %             scol = [ii ii-1 ii];
        %             sval = [C L T];
        %             Stencil = sparse(srow,scol,sval);
    elseif ii == 1 && jj == NY-2
        % Top Left Node
        % Slower
        Stencil(jj,ii) = C;     %-4;
        Stencil(jj-1,ii) = B;   %1;
        Stencil(jj,ii+1) = R;   %1;
        
        % Faster
        %             srow = [jj jj-1 jj];
        %             scol = [ii ii ii+1];
        %             sval = [C B R];
        %             Stencil = sparse(srow,scol,sval);
    elseif ii == NX-2 && jj == NY-2
        % Top Right Node
        % Slower
        Stencil(jj,ii) = C;     %-4;
        Stencil(jj,ii-1) = L;   %1;
        Stencil(jj-1,ii) = B;   %1;
        
        % Faster
        %             srow = [jj jj jj-1];
        %             scol = [ii ii-1 ii];
        %             sval = [C L B];
        %             Stencil = sparse(srow,scol,sval);
    elseif ii == 1
        %             Left Node
        % Slower
        Stencil(jj,ii) = C;     %-4;
        Stencil(jj,ii+1) = R;   %1;
        Stencil(jj+1,ii) = T;   %1;
        Stencil(jj-1,ii) = B;   %1;
        
        % Faster
        %             srow = [jj jj jj+1 jj-1];
        %             scol = [ii ii+1 ii ii];
        %             sval = [C R T B];
        %             Stencil = sparse(srow,scol,sval);
    elseif ii == NX-2
        % Right Node
        % Slower
        Stencil(jj,ii) = C;     %-4;
        Stencil(jj,ii-1) = L;   %1;
        Stencil(jj+1,ii) = T;   %1;
        Stencil(jj-1,ii) = B;   %1;
        
        % Faster
        %             srow = [jj jj jj+1 jj-1];
        %             scol = [ii ii-1 ii ii];
        %             sval = [C L T B];
        %             Stencil = sparse(srow,scol,sval);
    elseif jj == 1
        % Bottom Node
        % Slower
        Stencil(jj,ii) = C;     %-4;
        Stencil(jj+1,ii) = T;   %1;
        Stencil(jj,ii-1) = L;   %1;
        Stencil(jj,ii+1) = R;   %1;
        
        % Faster
        %             srow = [jj jj+1 jj jj];
        %             scol = [ii ii ii-1 ii+1];
        %             sval = [C T L R];
        %             Stencil = sparse(srow,scol,sval);
    elseif jj == NY-2
        % Top Node
        % Slower
        Stencil(jj,ii) = C;     %-4;
        Stencil(jj-1,ii) = B;   %1;
        Stencil(jj,ii+1) = R;   %1;
        Stencil(jj,ii-1) = L;   %1;
        
        % Faster
        %             srow = [jj jj-1 jj jj];
        %             scol = [ii ii ii+1 ii-1];
        %             sval = [C B R L];
        %             Stencil = sparse(srow,scol,sval);
    else
        % Interior Node
        %Slower
        Stencil(jj,ii) = C;     %-4;
        Stencil(jj,ii+1) = R;   %1;
        Stencil(jj,ii-1) = L;   %1;
        Stencil(jj+1,ii) = T;   %1;
        Stencil(jj-1,ii) = B;   %1;
        
        %Faster
        %             srow = [jj jj jj jj+1 jj-1];
        %             scol = [ii ii+1 ii-1 ii ii];
        %             sval = [C R L T B];
        %             Stencil = sparse(srow,scol,sval);
    end
    %flipud(Stencil);  % Just for viz purposes
    [indx] = find(Stencil);
    val = Stencil(indx);
    
    sRow{node} = ones(length(val),1)*node;
    sCol{node} = indx;
    sVal{node} = val;
end

sRow = cell2mat(reshape(sRow,(NX-2)*(NY-2),1));
sCol = cell2mat(reshape(sCol,(NX-2)*(NY-2),1));
sVal = cell2mat(reshape(sVal,(NX-2)*(NY-2),1));
A = sparse(sRow,sCol,sVal);
% spy(A)