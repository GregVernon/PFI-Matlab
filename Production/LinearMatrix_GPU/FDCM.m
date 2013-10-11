function [AMat] = FDCM(Nx,My,dx,dy)

col = zeros(5*Nx^2,1);
row = zeros(5*Nx^2,1);
Val = zeros(5*Nx^2,1);

kk = 0;
kk2 = 0;
for jj = 2:My-1
    for ii = 2:Nx-1
        kk = kk+1;
        kk2 = kk2+1;
        
%         if ~rem(kk2,2^13)
%             clc
%             disp('Constructing Coefficient Matrix:')
%             disp(['Evaluating Node ' num2str(kk2) ' of ' num2str(Nx^2)])
%         end
        
        row(kk) = kk2;
        col(kk) = kk2;
        Val(kk) = -2*(dx^2 + dy^2);
        if ii == 2 && jj == 2           % Node A
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 + 1;
            Val(kk) = dx^2;
            
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 + (Nx-2);
            Val(kk) = dy^2;
            
        elseif ii == 2 && jj == My-1     % Node M
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 + 1;
            Val(kk) = dx^2;
            
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 - (Nx-2);
            Val(kk) = dy^2;
            
        elseif ii == Nx-1 && jj == 2     % Node D
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 - 1;
            Val(kk) = dx^2;
            
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 + (Nx-2);
            Val(kk) = dy^2;
            
        elseif ii == Nx-1 && jj == My-1   % Node P
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 - 1;
            Val(kk) = dx^2;
            
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 -(Nx-2);
            Val(kk) = dy^2;
            
        elseif ii == 2                  % Nodes E & I
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 + 1;
            Val(kk) = dx^2;
            
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 + (Nx-2);
            Val(kk) = dy^2;
            
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 - (Nx-2);
            Val(kk) = dy^2;
            
        elseif ii == Nx-1                % Nodes H & L
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 - 1;
            Val(kk) = dx^2;
            
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 - (Nx-2);
            Val(kk) = dy^2;
            
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 + (Nx-2);
            Val(kk) = dy^2;
            
        elseif jj == 2                  % Nodes B & C
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 + 1;
            Val(kk) = dx^2;
            
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 - 1;
            Val(kk) = dx^2;
            
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 + (Nx-2);
            Val(kk) = dy^2;
            
        elseif jj == My-1                % Nodes N & O
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 + 1;
            Val(kk) = dx^2;
            
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 - 1;
            Val(kk) = dx^2;
            
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 - (Nx-2);
            Val(kk) = dy^2;
            
        else                            % Nodes F & G & J & K
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 + 1;
            Val(kk) = dx^2;
            
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 - 1;
            Val(kk) = dx^2;
            
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 + (Nx-2);
            Val(kk) = dy^2;
            
            kk = kk+1;
            row(kk) = kk2;
            col(kk) = kk2 - (Nx-2);
            Val(kk) = dy^2;
            
        end
    end
end
trimcol = col==0;
col(trimcol) = [];
row(trimcol) = [];
Val(trimcol) =[];
AMat = sparse(row,col,Val);
% figure
% spy(A)
