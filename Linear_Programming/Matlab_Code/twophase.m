function [ x,v ] = twophase( A, b, c )
% ���׶η�
n=length(c);
m=length(b);
%ѡ���������
v=n*ones(m,1);
for i=1:m
    v(i)=v(i)+i;
end
%�������˹������ĵ�һ�׶�
[x,v]=simplex([A eye(m)],b,[zeros(n,1);ones(m,1)],v);
%ȥ���˹����� ����ڶ��׶�
Binv=inv(A(:,v));
A=Binv*A;
b=Binv*b;
[x,v]=simplex(A,b,c,v);
end