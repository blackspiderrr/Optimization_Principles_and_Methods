function [x,v,Binv]=revsimp(A,b,c,v,Binv)
n=length(c);
m=length(b);
% 由A选取的列向量组成的基矩阵
cB=c(v(:));
% 计算当前的检验数
y0 = Binv*b; 
r = c'-cB'*Binv*A; 
options(5) = 1;%进基向量的选取规则
while ones(1,n)*(r' >= zeros(n,1)) ~= n %大于0的检验数小于变量维数n
    % 从小于0的检验数中选择最小检验数
    if options(5) == 0
        [r_q,q] = min(r);
    else
        q=1;
        while r(q) >= 0
            q=q+1;
        end
    end 
    % 计算yq
    yq = Binv*A(:,q);
    min_ratio = inf;
    p=0;
    % 如果不存在yi大于0，则算法停止问题有无界解；否则计算
    for i=1:m
        if yq(i)>0
            if y0(i)/yq(i) < min_ratio
                min_ratio = y0(i)/yq(i);
                p = i;
            end %if
        end %if
    end %for
    if p == 0
        disp('Problem unbounded');
        break;
    end %if
    %构造增广的修正型单纯形状表
    augrevtabl=pivot([Binv y0 yq],p,m+2);
    Binv=augrevtabl(:,1:m);
    y0=augrevtabl(:,m+1);
    v(p) = q;
    cB=c(v(:)); 
    r = c'-cB'*Binv*A; %row vector of relative cost coefficients
end
x=zeros(n,1);
x(v(:))=y0;
end