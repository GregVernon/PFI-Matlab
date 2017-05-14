function A = assembleCoeffMat(NX, NY, dx, dy)
%% Generate A for a domain from 1:N
% STENCIL = spalloc(N^2,N^2,5*N^2);
StencilTemplate = spalloc(NY-2,NX-2,5);
% StencilTemplate = sym(zeros(N));
sRow = zeros(5*(NX-2)*(NY-2),1);
sCol = zeros(5*(NX-2)*(NY-2),1);
sVal = zeros(5*(NX-2)*(NY-2),1);
% syms C L R B T
C = -2*(dy^2 + dx^2); %-4;
L = dy^2;%1;
R = dy^2;%1;
B = dx^2;%1;
T = dx^2;%1;
kk = 0;
node = 0;
for ii = 1:NX-2
    for jj = 1:NY-2
        Stencil = StencilTemplate;
        node = node+1;
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
        % Unvectorized
        %         for iVal = 1:length(val)
        %             kk = kk+1;
        %             sRow(kk) = node;
        %             sCol(kk) = indx(iVal);
        %             sVal(kk) = val(iVal);
        %         end
        
        % Vectorized
        KK = kk+1:kk+length(val);
        sRow(KK) = node;
        sCol(KK) = indx;
        sVal(KK) = val;
        kk = kk+length(val);
        %tmp = transpose(Stencil);
        %tmp1 = reshape(tmp,1,N^2);
        %STENCIL(kk,:) = reshape(transpose(Stencil),1,N^2);
    end
end
sRow(kk+1:end) = [];
sCol(kk+1:end) = [];
sVal(kk+1:end) = [];
A = sparse(sRow,sCol,sVal);
% spy(A)