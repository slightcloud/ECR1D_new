function phi_r=G_S_SOR(phi,rho,tol,w,N,dz,nz)
%*******************This function calculate phi**********************
global EPS0
R_freq=100;

for j=1:N     %��������
    
    phi(1)=phi(2);
    for i=2:nz-1  %�м���ĵ���
        phi_tem=0.5*(rho(i)/EPS0*dz*dz+phi(i-1)+phi(i+1));   %�������޲�ַ����õ�
        phi(i)=phi(i)+w*(phi_tem-phi(i));
        
    end
    
    %             phi_tem(1)=phi(2);     %��ڴ�Ϊ��һ��߽�����
    %             phi_tem(nz)=phi(nz-1); %���ڴ�ҲΪ��һ��߽�����
    
    if mod(j,R_freq)==0
        r_i=phi(2)-phi(1);
        r_sum=r_i*r_i;
        for i=2:nz-1
            r_i=(phi(i+1)-2*phi(i)+phi(i-1))/(dz*dz)+rho(i)/EPS0;
            r_sum=r_sum+r_i*r_i;
        end
        R=sqrt(r_sum)/nz;
%         R=norm(phi_tem-phi); %������ε����Ĳв�
        if (R<=tol)          %��������
%             phi=phi_tem;
            fprintf('Congratulations! The solver of the electric potential is converged at %d step\n',j);
            break;
        end
    end
%     phi=phi_tem;
end
if j>=N                  %���������ⲻ����������ֹ���������
    fprintf('The result of the electric potential is not convergent, the program is stopped!\n');
%     break;
end
phi_r=phi;
