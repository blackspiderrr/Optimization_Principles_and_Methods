function [z]=admmfast(Q,f,x,G,W,S,L)

k=1;
maxiter = 100;
iMG = inv(Q)*G';
iMc = inv(Q)*f*x;
GL=1/L*G;
bL=1/L*(W+S*x);

y0 = zeros(80,1);
y = zeros(80,1);
while k<maxiter
beta=max((k-1)/(k+2),0);
w=y+beta*(y-y0);
z=-(iMG*w+iMc);
s=GL*z-bL;
y0=y;
y=w+s;
k=k+1;
end
z= -inv(Q)*(f*x + G'*y);
end