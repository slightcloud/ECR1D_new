%****************************************************************
%该子函数用高斯-赛德尔加超松弛迭代的方法求解电势%
%****************************************************************

function phi = GS_SOR(rho,nz,tol,w,choice)
%***************************************************************
%rho代表的是每个格点上的电荷密度；
%nz代表的是系统的格点总数；
%tol代表的是迭代结果的容差；
%w代表的是松弛因子
%choice的取值为1和2，
%          choice=1：表示两个边界的电势为0
%          choice=2：表示两个边界的电势对z的一阶导数为0
%***************************************************************
global dz EPS0    %声明全局变量
solver_it = 20000;  %迭代步数
phi = zeros(nz,1);  %给电势赋初值
%**********************choice = 1*******************************
if choice == 1
    for i = 1:solver_it
        g = 0.5*((rho(2:nz-1,1)/EPS0)*dz*dz+phi(1:nz-2,1)+phi(3:nz,1));   %只求解2到nz-1的电势值，依据是泊松方程
        phi(2:nz-1,1) = phi(2:nz-1,1)+w*(g-phi(2:nz-1,1));                %做超松弛迭代
        if mod(i,25) == 0
            res = rho(2:nz-1,1)/EPS0+(phi(1:nz-2,1)-2*phi(2:nz-1,1)+phi(3:nz,1))/(dz*dz); %求出迭代的残差
            sum_res = sum(res.^2);          %求出所有残差的和
            judge = sqrt(sum_res/nz);       %求出平均残差
            if judge < tol                  %将平均残差与设定容差进行对比，若满足残差要求则停止循环
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
        g = 0.5*((rho(2:nz-1,1)/EPS0)*dz*dz+phi(1:nz-2,1)+phi(3:nz,1));   %只求解2到nz-1的电势值，依据是泊松方程
        phi(2:nz-1,1) = phi(2:nz-1,1)+w*(g-phi(2:nz-1,1));                %做超松弛迭代
        phi(1,1) = phi(2,1);                                              %左边界的电势
        phi(nz,1) = phi(nz-1,1);                                          %右边界的电势
        if mod(i,25) == 0
            res = rho(2:nz-1,1)/EPS0+(phi(1:nz-2,1)-2*phi(2:nz-1,1)+phi(3:nz,1))/(dz*dz); %求出迭代的残差
            sum_res = sum(res.^2);          %求出所有残差的和
            judge = sqrt(sum_res/nz);       %求出平均残差
            if judge < tol                  %将平均残差与设定容差进行对比，若满足残差要求则停止循环
                break;
            end
        end
    end
    if i == solver_it
        fprintf('The interation is not converged!!!\n');
    end
end
        

