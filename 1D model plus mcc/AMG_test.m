clear;
clc
nz=4404;
EPS0 = 8.854e-12;
dz=2.3509e-5;
A=blktridiag(2,-1,-1,nz);
B=full(A);
B(1,:)=0*B(1,:);
B(1,1)=1;
B(nz,:)=0*B(nz,:);
B(nz,nz)=1;
den(1)=0;
den(nz)=0;
den=den/EPS0*(dz)^2;
% [phi,k]=G_S(B,den,phi,1e-3,2000);

load den.mat 
%K=A;
K=B;
%F=b;
F=den;
w=zeros(size(B,1),1);
level=3;
node_order=[-1 1 1];
relax_it=2; 
relax_para=1;  
post_smoothing=1; 
max_iter=20000; 
tol=1e-06; 
pc_type=2;
connection_threshold=0.25;
 [w,error,iter,flag]=AMG(K,F, w, level,  relax_it, relax_para, ...
                         post_smoothing, max_iter, tol,  pc_type, connection_threshold);
