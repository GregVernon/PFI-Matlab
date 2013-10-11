The purpose of this repository is to contain code to be used as an algorithm template
for transcription of code into another language (Python - C - Fortran)

This code is written with matrices reshaped into vector form:

% (jj,ii)       -->  jj + Nx*(ii-1)
% (jj,ii+1)     -->  jj + Nx*ii
% (jj,ii-1)     -->  jj + Nx*(ii-2)
% (jj+1,ii)     -->  (jj+1) + Nx*(ii-1)
% (jj-1,ii)     -->  (jj-1) + Nx*(ii-1)

where jj is the steps related to y and ii to x. 

Remember that Matlab uses:  (Row,Column) format.