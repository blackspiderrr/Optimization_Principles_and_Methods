clear; clc; close all;
% ����ϵͳϵ������
A=[2 2 0; -2 -1 1;0 2 1.5]; B=[0; 0;1];
% ��ʼ״̬��-�����������һ���ص�Լ����Χ�ڣ��������޽�
x0=[2;1;1];
% Ԥ�ⲽ��
Np=10;
% �Ż�Ŀ���������Ȩ����
Q=eye(3); R=1;
% ת��Ϊ�ÿ�����ut��ʾ�ģ�����״̬�����Ƶ����̵ľ���
At=[]; Bt=[]; temp=[];
% ת����ļ�Ȩ����
Qt=[]; Rt=[];
%cost=[];
% ��Ȩ����ļ�����̣��Լ��Ƶ����̾���ĵ��ӹ���
for i=1:Np
        At=[At; A^i];
        Bt=[Bt zeros(size(Bt,1), size(B,2));
            A^(i-1)*B temp];
        temp=[A^(i-1)*B temp];
        
        Qt=[Qt zeros(size(Qt,1),size(Q,1));
            zeros(size(Q,1),size(Qt,1)) Q];
        Rt=[Rt zeros(size(Rt,1),size(R,1));
            zeros(size(R,1),size(Rt,1)) R];
  
end
% ������ut��������
%lb=-1*ones(Np,1);
%ub=5*ones(Np,1);
% ������ut�ĳ�ʼֵ
%u0=zeros(Np,1);
% ת������Ż�Ŀ�꺯������ѭ���Ż�������H��ı��ʽΪ�Ż�Ŀ�����һ��
H=2*(Bt'*Qt*Bt + Rt);
% ת������Ż��еĲ���ʽԼ�����ϵ�����󣬺���ѭ���е�biΪ����ʽ�ұ�
Ai=[Bt;-Bt;eye(10);-eye(10)];
%Ai=[Bt;eye(10)];
% ����u������ÿһ�����õĿ�����
d=[20;10;20; 20;10;20; 20;10;20 ;20;10;20 ;20;10;20 ;20;10;20 ;20;10;20 ;20;10;20 ;20;10;20 ;20;10;20];
u=[];
x=x0;
xk=x0;
[a,b]=eig(Ai*inv(H)*Ai');
lamada = max(diag(b))+1;
h=1;
S=zeros(80,1);
for k=1:50
    % һ��׼�����������ж����Ż�
    bi=[d-At*xk; d+At*xk;10*ones(10,1);10*ones(10,1)];
    %ub=[d-At*xk;10*ones(10,1)];
    %lb=[d-At*xk;10*ones(10,1)];
    [ut]=admmfast(H,(2*xk'*At'*Qt*Bt)',h,Ai,bi,S,lamada);
    %[ut]=admmqp(H,(2*xk'*At'*Qt*Bt)',Ai,);
    % ��ʾ������Ƿ�����
    %fprintf('%d\n', exitflag)

    % �����Ż��õ��Ŀ������ĵ�һ��Ԫ����Ϊʵ�����õĿ����������뵽ԭϵͳ�еõ���һ��ʱ�̵�״̬��
    u(k) = ut(1);
    x(:, k+1) = A*x(:, k) + B*u(k);
    xk = x(:, k+1);
    % ���Ż���ʼֵ�����޸ģ�����Ԥ��ֵ�ĺ����Ϊ��һ���ĳ�ʼֵ
   % u0 = [ut(2:Np); ut(Np)];
    
end


figure();
plot(x');
legend('x_1','x_2','x_2');

figure();
plot(u);
legend('u');
