function [x,N] = affscale(A,b,c,u)
% 仿射尺度法
% u:初始点
xnew=u;
n=length(c);
% 设置迭代次数
max_iter=8;
% 设置合适的步长
alpha = 0.99;
for k = 1:max_iter
    xcurr=xnew;
    D = diag(xcurr);
    Abar = A*D;
    % 构造正交投影算子
    Pbar = eye(n) - Abar'*inv(Abar*Abar')*Abar;
    d = -D*Pbar*D*c;
    % 确定合适的步长
    if d ~= zeros(n,1)
        nonzd = find(d<0);
        r = min(-xcurr(nonzd)./d(nonzd));
    else
        disp('Terminating: d = 0');
        break;
    end
    % 迭代序列
    xnew = xcurr + alpha*r*d; 
end
% 判断输出条件
if nargout >= 1
    x=xnew;
    if nargout == 2
        N=k;
    end
else
    disp('Final point =');
    disp(xnew');
    disp('Number of iterations =');
    disp(k);
end
end