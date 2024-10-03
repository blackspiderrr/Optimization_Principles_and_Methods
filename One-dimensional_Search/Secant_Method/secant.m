function [x,v] = secant(g,x1,x2,e)
%使用割线法求解求解g(x)=0
%确认参数--割线法迭代--更新变量   
if nargin < 4
    e = 1e-5;
    if nargin < 3
        x1 = 0;
        x2 = 1;
        if nargin < 1
            disp('至少需要提供函数g');
        end
    end
end
while abs(x2-x1) >= abs(x1)*e
    x0 = x1;
    x1 = x2;
    g0 = g(x0);
    g1 = g(x1);
    x2 = x1 - g1*(x1-x0)/(g1-g0);
end 
x = x2;
v = g(x);
fprintf('x=%g g(x)=%g\n', x, v);
end