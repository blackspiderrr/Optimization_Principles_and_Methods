function [x,N] = affscale(A,b,c,u)
% ����߶ȷ�
% u:��ʼ��
xnew=u;
n=length(c);
% ���õ�������
max_iter=8;
% ���ú��ʵĲ���
alpha = 0.99;
for k = 1:max_iter
    xcurr=xnew;
    D = diag(xcurr);
    Abar = A*D;
    % ��������ͶӰ����
    Pbar = eye(n) - Abar'*inv(Abar*Abar')*Abar;
    d = -D*Pbar*D*c;
    % ȷ�����ʵĲ���
    if d ~= zeros(n,1)
        nonzd = find(d<0);
        r = min(-xcurr(nonzd)./d(nonzd));
    else
        disp('Terminating: d = 0');
        break;
    end
    % ��������
    xnew = xcurr + alpha*r*d; 
end
% �ж��������
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