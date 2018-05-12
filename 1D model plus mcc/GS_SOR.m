%****************************************************************
%���Ӻ����ø�˹-���¶��ӳ��ɳڵ����ķ���������%
%****************************************************************

function phi = GS_SOR(rho,nz,tol,w,choice)
%***************************************************************
%rho�������ÿ������ϵĵ���ܶȣ�
%nz�������ϵͳ�ĸ��������
%tol������ǵ���������ݲ
%w��������ɳ�����
%choice��ȡֵΪ1��2��
%          choice=1����ʾ�����߽�ĵ���Ϊ0
%          choice=2����ʾ�����߽�ĵ��ƶ�z��һ�׵���Ϊ0
%***************************************************************
global dz EPS0    %����ȫ�ֱ���
solver_it = 20000;  %��������
phi = zeros(nz,1);  %�����Ƹ���ֵ
%**********************choice = 1*******************************
if choice == 1
    for i = 1:solver_it
        g = 0.5*((rho(2:nz-1,1)/EPS0)*dz*dz+phi(1:nz-2,1)+phi(3:nz,1));   %ֻ���2��nz-1�ĵ���ֵ�������ǲ��ɷ���
        phi(2:nz-1,1) = phi(2:nz-1,1)+w*(g-phi(2:nz-1,1));                %�����ɳڵ���
        if mod(i,25) == 0
            res = rho(2:nz-1,1)/EPS0+(phi(1:nz-2,1)-2*phi(2:nz-1,1)+phi(3:nz,1))/(dz*dz); %��������Ĳв�
            sum_res = sum(res.^2);          %������вв�ĺ�
            judge = sqrt(sum_res/nz);       %���ƽ���в�
            if judge < tol                  %��ƽ���в����趨�ݲ���жԱȣ�������в�Ҫ����ֹͣѭ��
                break;
            end
        end
    end
    if i == solver_it
        fprintf('The interation is not converged!!!\n');
    end
end

%**********************choice = 2*******************************
if choice == 2
    for i = 1:solver_it
        g = 0.5*((rho(2:nz-1,1)/EPS0)*dz*dz+phi(1:nz-2,1)+phi(3:nz,1));   %ֻ���2��nz-1�ĵ���ֵ�������ǲ��ɷ���
        phi(2:nz-1,1) = phi(2:nz-1,1)+w*(g-phi(2:nz-1,1));                %�����ɳڵ���
        phi(1,1) = phi(2,1);                                              %��߽�ĵ���
        phi(nz,1) = phi(nz-1,1);                                          %�ұ߽�ĵ���
        if mod(i,25) == 0
            res = rho(2:nz-1,1)/EPS0+(phi(1:nz-2,1)-2*phi(2:nz-1,1)+phi(3:nz,1))/(dz*dz); %��������Ĳв�
            sum_res = sum(res.^2);          %������вв�ĺ�
            judge = sqrt(sum_res/nz);       %���ƽ���в�
            if judge < tol                  %��ƽ���в����趨�ݲ���жԱȣ�������в�Ҫ����ֹͣѭ��
                break;
            end
        end
    end
    if i == solver_it
        fprintf('The interation is not converged!!!\n');
    end
end
        

