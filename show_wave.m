clc;  clear; close('all');
fid=fopen('C:\Users\PC\Documents\Visual Studio 2015\Projects\MPI_test\test.bin');
%domain_size=fread(fid,2,'double');
res=fread(fid,[198,198],'double');
contour(res);
