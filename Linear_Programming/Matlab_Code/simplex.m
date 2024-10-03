function [x, v] = simplex( A, b, c ,v)
% �����η���ʵ��
%����:ϵ������ Լ�������ĳ�����������Ŀ�꺯����ϵ��������ѡȡ�Ļ��������
%��������Ž⡢���������������
n=length(c);
m=length(b);
%�����������淶��
cB=c(v(:));
r = c'-cB'*A;
cost = -cB'*b;
tabl=[A b;r cost];

%OPTIONS(5) specifies how the pivot element is selected;
% 0=choose the most negative relative cost coefficient;
% 1=use Bland��s rule.
options(5) = 0;
while ones(1,n)*(r' >= zeros(n,1)) ~= n
    % ѡ����ʽ��������ķ�ʽ,ѡ���С��0�ļ�����
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
        % ���������y_iq����0�����������޽�⣻�������p=argmin_i {y_i0/y_iq:y_iq>0}
        for i=1:m
            if tabl(i,q)>0
                if tabl(i,n+1)/tabl(i,q) < min_ratio %table ���һ��Ϊy_0
                    min_ratio = tabl(i,n+1)/tabl(i,q);
                    p = i;
                end %if
            end %if
        end %for
        if p == 0
            disp('Problem unbounded');
            break;
        end %if
        %��Ԫ��(p,q)Ϊ����Ԫ�ؽ�������任�������������淶��
        tabl=pivot(tabl,p,q);
        v(p) = q;
        r = tabl(m+1,1:n);
end %while
    x=zeros(n,1);
    x(v(:))=tabl(1:m,n+1);
end
