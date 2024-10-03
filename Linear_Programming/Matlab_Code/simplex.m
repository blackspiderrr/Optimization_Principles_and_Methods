function [x, v] = simplex( A, b, c ,v)
% 单纯形法的实现
%输入:系数矩阵、 约束函数的常数列向量、目标函数的系数向量、选取的基向量序号
%输出：最优解、基向量的序号向量
n=length(c);
m=length(b);
%构造增广矩阵规范性
cB=c(v(:));
r = c'-cB'*A;
cost = -cB'*b;
tabl=[A b;r cost];

%OPTIONS(5) specifies how the pivot element is selected;
% 0=choose the most negative relative cost coefficient;
% 1=use Bland’s rule.
options(5) = 0;
while ones(1,n)*(r' >= zeros(n,1)) ~= n
    % 选择合适进基向量的方式,选择出小于0的检验数
    if options(5) == 0
        [r_q,q] = min(r);
    else
        q=1;
        while r(q) >= 0
            q=q+1;
        end
    end
        min_ratio = inf;
        p=0;
        % 如果不存在y_iq大于0，则问题有无界解；否则计算p=argmin_i {y_i0/y_iq:y_iq>0}
        for i=1:m
            if tabl(i,q)>0
                if tabl(i,n+1)/tabl(i,q) < min_ratio %table 最后一列为y_0
                    min_ratio = tabl(i,n+1)/tabl(i,q);
                    p = i;
                end %if
            end %if
        end %for
        if p == 0
            disp('Problem unbounded');
            break;
        end %if
        %以元素(p,q)为枢轴元素进行枢轴变换，更新增广矩阵规范型
        tabl=pivot(tabl,p,q);
        v(p) = q;
        r = tabl(m+1,1:n);
end %while
    x=zeros(n,1);
    x(v(:))=tabl(1:m,n+1);
end
