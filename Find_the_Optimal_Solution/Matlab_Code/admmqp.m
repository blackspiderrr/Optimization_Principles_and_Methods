function [z]=admmqp(H,f,A,u,l)
n = size(H,1);
z = ones(n,1);
s = ones(n,1);
v = ones(n,1);
maxiter = 1000;
rho = 1.0;
iM = inv(H+rho*A'*A);
k=0;

    while k<maxiter
         k=k+1;
         z=-iM*(f+A'*(rho*(v-s)));
         Az=A*z;
        s=min(max(Az+v,l),u);
        v=v+Az-s;
    end





end