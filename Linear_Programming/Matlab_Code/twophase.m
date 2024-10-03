function [ x,v ] = twophase( A, b, c )
% 两阶段法
n=length(c);
m=length(b);
%选择进基向量
v=n*ones(m,1);
for i=1:m
    v(i)=v(i)+i;
end
%求解加入人工变量的第一阶段
[x,v]=simplex([A eye(m)],b,[zeros(n,1);ones(m,1)],v);
%去掉人工变量 计算第二阶段
Binv=inv(A(:,v));
A=Binv*A;
b=Binv*b;
[x,v]=simplex(A,b,c,v);
end