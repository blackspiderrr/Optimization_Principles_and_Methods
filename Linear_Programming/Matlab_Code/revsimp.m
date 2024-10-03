function [x,v,Binv]=revsimp(A,b,c,v,Binv)
n=length(c);
m=length(b);
% ��Aѡȡ����������ɵĻ�����
cB=c(v(:));
% ���㵱ǰ�ļ�����
y0 = Binv*b; 
r = c'-cB'*Binv*A; 
options(5) = 1;%����������ѡȡ����
while ones(1,n)*(r' >= zeros(n,1)) ~= n %����0�ļ�����С�ڱ���ά��n
    % ��С��0�ļ�������ѡ����С������
    if options(5) == 0
        [r_q,q] = min(r);
    else
        q=1;
        while r(q) >= 0
            q=q+1;
        end
    end 
    % ����yq
    yq = Binv*A(:,q);
    min_ratio = inf;
    p=0;
    % ���������yi����0�����㷨ֹͣ�������޽�⣻�������
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
    %��������������͵�����״��
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