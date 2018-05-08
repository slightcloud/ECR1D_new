function[x,k]=G_S(A,b,x0,eps,M)
%��˹���¶��������󷽳���Ľ⣨����ʽ��⣩
%AΪ�������ϵ������bΪ��������Ҷ���
%xΪ���Է�����Ľ��ˣ�x0Ϊ������ֵ
%epsΪ����ޣ�MΪ������������
if nargin==3
    eps= 1.0e-6;%Ĭ�Ͼ���
    M = 10000;%��������ʱĬ�Ϻ���������
elseif nargin ==4
    M = 10000;%������Ĭ��ֵ
elseif nargin<3
    error('��������');
    return
end
[n,m]=size(A);
nb=length(b);
%�������������е�ά�������ʱ��ֹͣ���㣬�����������Ϣ
if n~=m
    error('����A�����������������!');
    return;
end
%�����������Ҷ����ά����ƥ��ʱ��ֹͣ���㣬�����������Ϣ
if n~=nb
    error('����A�����������b�ĳ������!');
    return;
end 
L =zeros(n,n);
U =zeros(n,n);
D =zeros(n,n);
for i=2:n
    for j=1:i-1
        L(i,j)=-A(i,j);
    end
end
for i=1:n-1
    for j=i+1:n
        U(i,j)=-A(i,j);
    end
end
for i=1:n
    D(i,i)=A(i,i);
end
B=inv(D-L)*U;        %BΪ��������
g=inv(D-L)*b;        %gΪ�Ҷ���
pr=max(abs(eig(B))); %����������װ뾶
if pr>=1
    error('���������װ뾶����1������������');
    return;
end
k=0;
tol=1;
while tol>=eps 
    x = B*x0+g;
    k = k+1;         %��������
    tol = norm(x-x0);%ǰ������������������
    x0 = x;
    if(k>=M)
        disp('Warning: ��������̫�࣬���ܲ�������');
        return;
    end
end
