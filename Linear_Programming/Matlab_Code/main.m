clc;clear clc;

A = [1 0 1 0 0; 0 1 0 1 0; 1 1 0 0 1];
% A为约束方程组系数矩阵
b = [4;6;8];
%b为约束方程组常数项
c = [-2;-5;0;0;0];
u = [2;3;2;3;3];
[x,N] = affscale(A,b,c,u)

%{
 A=[1 0 1 0 0; 0 1 0 1 0; 1 1 0 0 1];
 b=[4;6;8];
 c=[-2;-5;0;0;0];
 v=[3;4;5];
 Binv=eye(3);
[x,v,Binv] = revsimp(A,b,c,v,Binv)
%}
%{
 A=[0 5 1 0 0; 6 2 0 1 0;1 1 0 0 1];
 b=[15;24;5];
 c=[-2;-1;0;0;0]';
 [ X, z,BIndex ] = simplex( A, b, c );
%}
%{
 A=[0 5 1 0 0; 6 2 0 1 0;1 1 0 0 1];
 b=[15;24;5];
 c=[-2;-1;0;0;0]';
 [X1,z1,v1,X2, z2,v2 ] = twophase( A, b, c )
%}