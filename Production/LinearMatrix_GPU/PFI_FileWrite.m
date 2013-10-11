function PFI_FileWrite(filename,timestep,formats,My,Omega,Psi,U,V)

for fout = 1:length(formats)
    if strcmpi(formats{fout},'ASCII')
        fid = fopen(strcat(filename,'.pfi'),'a');
        fprintf(fid,'%s%i\n','TIMESTEP: ',timestep);
        
        %         for nvar = 1:length(varargin{1})
        
        %             if strcmpi(varargin{1}(nvar),'Omega')
        fprintf(fid,'%s\n','Omega:');
        for jj = 1:My
            fprintf(fid,'%E\t',Omega(jj,:));
            fprintf(fid,'\n');
        end
        %             end
        
        %             if strcmpi(varargin{1}(nvar),'Psi')
        fprintf(fid,'%s\n','Psi:');
        for jj = 1:My
            fprintf(fid,'%E\t',Psi(jj,:));
            fprintf(fid,'\n');
        end
        %             end
        
        %             if strcmpi(varargin{1}(nvar),'U')
        fprintf(fid,'%s\n','U:');
        for jj = 1:My
            fprintf(fid,'%E\t',U(jj,:));
            fprintf(fid,'\n');
        end
        %             end
        
        %             if strcmpi(varargin{1}(nvar),'V')
        fprintf(fid,'%s\n','V:');
        for jj = 1:My
            fprintf(fid,'%E\t',V(jj,:));
            fprintf(fid,'\n');
        end
        %             end
        %         end
        fclose('all');
    elseif strcmpi(formats{fout},'MAT')
        OutVars = genvarname(strcat('Timestep_',num2str(timestep,'%i')));
        evalc([OutVars '= struct(''Omega'',Omega,''Psi'',Psi,''U'',U,''V'',V)''']);
        save([filename '.mat'],OutVars,'-append','-mat')
    end
end