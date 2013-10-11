%Output file name
tic
filename=('TestFile.vtk');
Omega = gather(Omega);
Psi = gather(Psi);
U = gather(U);
V = gather(V);
xx = gather(xx);
yy = gather(yy);
%load the MATLAB 3D flow example
% load wind
% tic

nr_of_elements=numel(U);
fid = fopen(filename, 'w'); 

%ASCII file header
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'VTK from Matlab\n');
fprintf(fid, 'BINARY\n\n');
fprintf(fid, 'DATASET STRUCTURED_GRID\n');
fprintf(fid, ['DIMENSIONS ' num2str(size(U,1)) ' ' num2str(size(U,2)) ' ' num2str(size(U,3)) '\n']);
fprintf(fid, ['POINTS ' num2str(nr_of_elements) ' float\n']);
fclose(fid);

%append binary x,y,z data
fid = fopen(filename, 'a'); 
fwrite(fid, [reshape(xx,1,nr_of_elements);  reshape(yy,1,nr_of_elements); zeros(1,nr_of_elements)],'float','b');

%append another ASCII sub header
fprintf(fid, ['\nPOINT_DATA ' num2str(nr_of_elements) '\n']);
fprintf(fid, 'VECTORS velocity_vectors float\n');

%append binary u,v,w data
fwrite(fid, [reshape(U,1,nr_of_elements);  reshape(V,1,nr_of_elements); zeros(1,nr_of_elements)],'float','b');

%append another binary u,v,w data set
% fprintf(fid, '\nVECTORS another_vector_set float\n'); %ASCII header
% fwrite(fid, [reshape(U*10,1,nr_of_elements);  reshape(V*2,1,nr_of_elements); reshape(zeros(My*Nx),1,nr_of_elements)],'float','b'); %binary data

%append some scalar data
fprintf(fid, '\nSCALARS VelocityMagnitude float\n'); %ASCII header
fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
fwrite (fid, reshape(sqrt(U.^2+V.^2),1,nr_of_elements),'float','b'); %binary data

%append some scalar data
fprintf(fid, '\nSCALARS Omega float\n'); %ASCII header
fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
fwrite (fid, reshape(Omega,1,nr_of_elements),'float','b'); %binary data

%append some scalar data
fprintf(fid, '\nSCALARS Psi float\n'); %ASCII header
fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
fwrite (fid, reshape(Psi,1,nr_of_elements),'float','b'); %binary data
fclose(fid);
toc