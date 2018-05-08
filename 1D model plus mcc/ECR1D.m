%1 2/2D ECR plasma simulation using PIC-MCC
clear variables
clc

global q_e m_e EPS0

%定义基本参量
q_e=1.602e-19;              %电子电量，单位：C
AMU = 1.661e-27;            %质子质量，单位：kg
m_e=9.1e-31;                %电子质量，单位：kg
m_Li=7*AMU;                 %Li-7质量，单位：kg
EPS0 = 8.854e-12;           %真空中的介电常数，单位：F/m
MU0=4*pi()*1e-7;            %真空中的磁导率，单位：H/m
k=1.38e-23;                 %玻尔兹曼常数，单位：J/K
c=3e8;                      %光速，单位：m/s
%c=3e2;                     %将光速改小，以观察粒子对电磁波的影响（结果显示这样改有问题）
Te_average=10;              %平均电子温度，单位：eV
n_e_average=5e17;           %平均电子密度，单位：m^-3
n0=1e18;                    %中性原子密度，单位：m^-3
T=273.15+300;               %系统温度为300℃，单位：K
lamda=sqrt(EPS0*Te_average/(n0*q_e));    %德拜长度，单位：m
z_min=0;                                 %z方向的最小值，单位：m
L_real=0.1035;                           %离子源z方向的最大值，单位：m
dz=lamda;                                %设置空间间隔
dt_e=dz/(2*c);                             %设置时间间隔，为保证模拟精确度，保证dt<=dz/c，一般设为dz/(2*c)单位：s
dt_i=100*dt_e;                           %离子的运动时间设置为电子运动时间的100倍
nz=round(L_real/dz)+1;                   %格点数
z_max=dz*(nz-1);                         %模拟区域z方向的最大值
spwt=1e8;                                %每个宏粒子包含的真实粒子数
vth_e=sqrt(2*q_e*0.5/m_e);               %电子的热速度，假设电子的热能为0.5 eV，单位：m/s
N_e=100000;                                %初始电子数目为N_e
N_Li1=0;                                 %Li+的数目
N_Li2=0;                                 %Li2+的数目
N_Li3=0;                                 %Li3+的数目
N_Li_ex=0;                               %激发态Li原子的数目
tol=1e-6;                               %电势求解的精度值
sor=1.4;                                 %超松驰迭代的因子
max_part=1000000;                           %粒子容器的最大值
%B_ex_z=0.0875;                          %外部永磁铁产生的磁感应强度，单位：T
E0=1.5e8;                                %微波电场强度的幅值，单位：V/m（不确定）
B0=0.0001;                               %微波磁感应强度的幅值，单位:T（不确定）
w=2.45;                                  %微波的频率,单位:GHz
step_num=1000;                          %总共运行的时间步数
N_phi=5e5;                               %电势求解的迭代步数
R_i_max=7.8034e-7;                       %各个价态总电离率的最大值
E_i_Li1=5.4;                             %Li0->Li1的电离能，单位：eV
E_i_Li2=75.77;                           %Li1->Li2的电离能，单位：eV
E_i_Li3=122.664;                         %Li2->Li3的电离能，单位：eV
z=zeros(nz,1);                           %用于盛放每个格点的坐标
z(2:nz,1)=(1:nz-1)'*dz*1000;                  %每个格点的坐标值，单位：mm
B_ex_z=0.0001*(2595.7+32.4092*z-3.1038*z.^2-0.0130847*z.^3+0.00198809*z.^4-2.48996e-5*z.^5+8.9815e-8*z.^6);   %根据实验上测得的磁感应强度值进行拟合，单位：T
%plot(z,B_ex_z);
B_ir=[-7.97027	4.27075	-3.96662	1.78936	-0.400521	0.0346181
    -12.2121	9.29911	-6.78284	2.46302	-0.453047	0.0332396
    -34.5328	36.724	-21.6441	6.55211	-1.0141	0.0636916]; %fitting parameters of ionization rate


%预分配空间
phi=zeros(nz,1);               %电势
phi_tem=zeros(nz,1);           %电势求解的中间变量
B_mic=zeros(nz,3);             %微波的磁场
E_mic=zeros(nz,3);             %微波的电场
E_sta=zeros(nz,1);             %电子与离子自身产生的静电场
E_e=zeros(max_part,3);         %粒子位置处的电场强度
B_e=zeros(max_part,3);         %粒子位置处的磁感应强度
E_boun=zeros(step_num,3);      %用于设置吸收边界
B_boun=zeros(step_num,3);      %用于设置吸收边界
%den=zeros(nz,1);              %每个格点上的粒子数密度
%chg=zeros(nz,1);              %每个格点上的电荷比例
J=zeros(nz,3);                 %每个格点上的电流密度
pos_e=zeros(max_part,1);         %电子的位置信息
vel_e=zeros(max_part,3);         %电子的速度信息
pos_Li1=zeros(max_part,1);    %Li+的位置信息
vel_Li1=zeros(max_part,3);    %Li+的速度信息
pos_Li2=zeros(max_part,1);    %Li2+的位置信息
vel_Li2=zeros(max_part,3);    %Li2+的速度信息
pos_Li3=zeros(max_part,1);    %Li3+的位置信息
vel_Li3=zeros(max_part,3);    %Li3+的速度信息
pos_Li_ex=zeros(max_part,1);    %激发态的Li原子的位置信息
vel_Li_ex=zeros(max_part,3);    %激发态的Li原子的速度信息



%给定初始电子的位置和速度
pos_e(1:N_e)=rand(N_e,1)*z_max;                                                   %初始电子的位置在0到z_max之间

%vel(1:N_e/2,:)=vth_e*rand(N_e/2,3);
%vel(N_e/2+1:N_e,:)=-vth_e*rand(N_e/2,3);

%vel(1:N_e,:)=1e8;
vel_e(1:N_e,:)=vth_e*2*(rand(N_e,3)+rand(N_e,3)+rand(N_e,3)-1.5);                 %设置初始电子的速度在热速度的-3倍到3倍之间
vel_e(1:N_e,:)=UpdateVelocity(E_e(1:N_e,:),B_e(1:N_e,:),vel_e(1:N_e,:),-0.5*dt_e);    %将电子的速度向前移动半个时间步长

%开始主循环，即对时间的循环
for ts_i=1:step_num                 %运行的时间步数
    
    %chg=zeros(nz,4);                %用于存放每个格点分配到的电荷量，1~4列分别代表:电子，Li1+，Li2+，Li3+
    %den_vel_Li1=zeros(nz,3);        %用于存放每个格点上Li+的电流密度分量
    %den_vel_Li2=zeros(nz,3);        %用于存放每个格点上Li2+的电流密度分量
    %den_vel_Li3=zeros(nz,3);        %用于存放每个格点上Li3+的电流密度分量
    %count_Li=zeros(nz-1,3);            %用于计算每个单元格内的粒子数，1~3列分别代表Li+，Li2+和Li3+
    %count_g_Li=zeros(nz,3);            %用于计算每个格点上的粒子数密度，1~3列分别代表Li+，Li2+和Li3+
    
    %分配Li+的电荷和电流密度到临近格点上
    %for p=1:N_Li1
    
    %    fi=1+pos_Li1(p)/dz;                 %实际格点位置，为浮点数
    %    i=floor(fi);                        %对应整数格点位置
    %   hz=fi-i;                            %粒子与第i个格点之间的距离占距离步长的比例
    
    %         %将Li+电荷按照比例分配到临近的两个格点上
    %         chg(i,2)=chg(i,2)+(1-hz);      %分配到第i个格点上的电荷比例
    %         chg(i+1,2)=chg(i+1,2)+hz;      %分配到第i+1个格点上的电荷比例
    %
    %         count_Li(i,1)=count_Li(i,1)+1;       %对每个单元格内的Li+计数
    %
    %         %将速度按比例分配到临近的两个格点上
    %         den_vel_Li1(i,:)=den_vel_Li1(i,:)+spwt*q_e*vel_Li1(p,:)*(1-hz);   %分配到第i个格点上的速度
    %         den_vel_Li1(i+1,:)=den_vel_Li1(i,:)+spwt*q_e*vel_Li1(p,:)*hz;     %分配到第i+1个格点上的速度
    %     end
    %
    %     %分配Li2+的电荷和电流密度到临近格点上
    %     for p=1:N_Li2
    %
    %         fi=1+pos_Li2(p)/dz;                 %实际格点位置，为浮点数
    %         i=floor(fi);                        %对应整数格点位置
    %         hz=fi-i;                            %粒子与第i个格点之间的距离占距离步长的比例
    %
    %         %将Li2+电荷按照比例分配到临近的两个格点上
    %         chg(i,3)=chg(i,3)+(1-hz);      %分配到第i个格点上的电荷比例
    %         chg(i+1,3)=chg(i+1,3)+hz;      %分配到第i+1个格点上的电荷比例
    %
    %         count_Li(i,2)=count_Li(i,2)+1;       %对每个单元格内的Li2+计数
    %
    %         %将速度按比例分配到临近的两个格点上
    %         den_vel_Li2(i,:)=den_vel_Li2(i,:)+spwt*2*q_e*vel_Li2(p,:)*(1-hz);   %分配到第i个格点上的速度
    %         den_vel_Li2(i+1,:)=den_vel_Li2(i,:)+spwt*2*q_e*vel_Li2(p,:)*hz;     %分配到第i+1个格点上的速度
    %     end
    %
    %     %分配Li3+的电荷和电流密度到临近格点上
    %     for p=1:N_Li3
    %
    %         fi=1+pos_Li3(p)/dz;                 %实际格点位置，为浮点数
    %         i=floor(fi);                        %对应整数格点位置
    %         hz=fi-i;                            %粒子与第i个格点之间的距离占距离步长的比例
    %
    %         %将Li3+电荷按照比例分配到临近的两个格点上
    %         chg(i,4)=chg(i,4)+(1-hz);      %分配到第i个格点上的电荷比例
    %         chg(i+1,4)=chg(i+1,4)+hz;      %分配到第i+1个格点上的电荷比例
    %
    %         count_Li(i,3)=count_Li(i,3)+1;       %对每个单元格内的Li3+计数
    %
    %         %将速度按比例分配到临近的两个格点上
    %         den_vel_Li3(i,:)=den_vel_Li3(i,:)+spwt*3*q_e*vel_Li3(p,:)*(1-hz);   %分配到第i个格点上的速度
    %         den_vel_Li3(i+1,:)=den_vel_Li3(i,:)+spwt*3*q_e*vel_Li3(p,:)*hz;     %分配到第i+1个格点上的速度
    %     end
    %
    %     %计算每个格点上的粒子密度
    %     for count_i=2:nz-1
    %         count_g_Li(count_i,:)=0.5*(count_Li(count_i-1,:)+count_Li(count_i,:))/(2*dz);
    %     end
    %     count_g_Li(1,:)=count_Li(1,:)/dz;      %第一个格点的粒子密度
    %     count_g_Li(nz,:)=count_Li(nz-1,:)/dz;  %最后一个格点的粒子密度
    %
    %     J_i=1e5*(den_vel_Li1.*count_g_Li(:,1)+den_vel_Li2.*count_g_Li(:,2)+den_vel_Li3.*count_g_Li(:,3));  %离子在每个格点上的电流密度
    
    %开始电子的运动循环
    for ts_e=1:dt_i/dt_e                                                                 %离子运动一步，电子运动dt_i/dt_e步
        
        chg=zeros(nz,4);                %用于存放每个格点分配到的电荷量，1~4列分别代表:电子，Li1+，Li2+，Li3+
        den_vel_Li1=zeros(nz,3);        %用于存放每个格点上Li+的电流密度分量
        den_vel_Li2=zeros(nz,3);        %用于存放每个格点上Li2+的电流密度分量
        den_vel_Li3=zeros(nz,3);        %用于存放每个格点上Li3+的电流密度分量
        count_Li=zeros(nz-1,3);            %用于计算每个单元格内的粒子数，1~3列分别代表Li+，Li2+和Li3+
        count_g_Li=zeros(nz,3);            %用于计算每个格点上的粒子数密度，1~3列分别代表Li+，Li2+和Li3+
        distri_Li1=zeros(N_Li1,nz-1);     %用来盛放每个单元格内Li+的序数
        distri_Li2=zeros(N_Li2,nz-1);     %用来盛放每个单元格内Li2+的序数
        distri_Li3=zeros(N_Li3,nz-1);     %用来盛放每个单元格内Li3+的序数
        
        for p=1:N_Li1
            fi=1+pos_Li1(p)/dz;                 %实际格点位置，为浮点数
            i=floor(fi);                        %对应整数格点位置
            count_Li(i,1)=count_Li(i,1)+1;       %对每个单元格内的Li+计数
            distri_Li1(count_Li(i,1),i)=p;  %每个单元格内的Li+的序数
        end
        for p=1:N_Li2
            fi=1+pos_Li2(p)/dz;                 %实际格点位置，为浮点数
            i=floor(fi);                        %对应整数格点位置
            count_Li(i,2)=count_Li(i,2)+1;       %对每个单元格内的Li2+计数
            distri_Li2(count_Li(i,2),i)=p;  %每个单元格内的Li2+的序数
        end
        for p=1:N_Li3
            fi=1+pos_Li3(p)/dz;                 %实际格点位置，为浮点数
            i=floor(fi);                        %对应整数格点位置
            count_Li(i,3)=count_Li(i,3)+1;       %对每个单元格内的Li3+计数
            distri_Li3(count_Li(i,3),i)=p;  %每个单元格内的Li3+的序数
        end
        
        %首先进行MCC过程
        den_Li=count_Li/(100*dz);       %每个单元格内Li+，Li2+和Li3+的数密度，单位：cm-1
        nu_t=(n0+den_Li(:,1)+den_Li(:,2))*R_i_max; %每个单元格内总的电子碰撞频率，考虑要不要给锂离子乘以spwt???
        P_emax=1-exp(-nu_t*dt_e);                    %每个单元格内任意一个电子的最大碰撞概率
        count_e=zeros(nz-1,1);                     %用来对每个单元格内的电子进行计数
        distri_e=zeros(N_e,nz-1);                   %用来盛放每个单元格内的电子序数
        for p=1:N_e
            fi=1+pos_e(p)/dz;            %实际格点位置，为浮点数
            i=floor(fi);               %对应整数格点位置
            count_e(i)=count_e(i)+1;   %相应单元格内的计数器计数一次
            distri_e(count_e(i),i)=p;  %每个单元格内的电子序数
        end
        N_e_mcc=P_emax.*count_e;        %每个单元格内抽取的电子个数
        N_e_sam=zeros(ceil(max(N_e_mcc)),nz-1); %该数列用来存放每个单元随机抽出的电子的序号
        Ek_e=zeros(ceil(max(N_e_mcc)),nz-1);    %该矩阵用来存放抽样出来的电子的动能
        nu_Li1=zeros(ceil(max(N_e_mcc)),nz-1);  %该矩阵用来存放抽样出来的e+Li0->Li+的碰撞频率
        nu_Li2=zeros(ceil(max(N_e_mcc)),nz-1);  %该矩阵用来存放抽样出来的e+Li+->Li2+的碰撞频率
        nu_Li3=zeros(ceil(max(N_e_mcc)),nz-1);  %该矩阵用来存放抽样出来的e+Li2+->Li3+的碰撞频率
        for i=1:nz-1            %开始对每个单元格进行蒙特卡洛抽样
            N_e_sam(1:ceil(N_e_mcc(i)),i)=(randperm(count_e(i),ceil(N_e_mcc(i))))';  %randperm(n,k)函数用于从n个数中随机抽取不相同的k个数
            for j=1:N_e_mcc(i)
                Ek_e(j,i)=0.5*m_e*sum((vel_e(distri_e(N_e_sam(j,i),i),:)).^2)/q_e;  %算出每个抽取出来的电子的动能，并转换为eV单位
                nu_Li1(j,i)=n0*10^(B_ir(1,:)*[1;log10(Ek_e(j,i));(log10(Ek_e(j,i)))^2;(log10(Ek_e(j,i)))^3;(log10(Ek_e(j,i)))^4;(log10(Ek_e(j,i)))^5])/nu_t(i); %根据数据库中的数据拟合曲线求出抽样电子对应的产生Li+的碰撞频率
                nu_Li2(j,i)=den_Li(i,1)*10^(B_ir(2,:)*[1;log10(Ek_e(j,i));(log10(Ek_e(j,i)))^2;(log10(Ek_e(j,i)))^3;(log10(Ek_e(j,i)))^4;(log10(Ek_e(j,i)))^5])/nu_t(i); %产生Li2+的碰撞频率
                nu_Li3(j,i)=den_Li(i,2)*10^(B_ir(3,:)*[1;log10(Ek_e(j,i));(log10(Ek_e(j,i)))^2;(log10(Ek_e(j,i)))^3;(log10(Ek_e(j,i)))^4;(log10(Ek_e(j,i)))^5])/nu_t(i); %产生Li3+的碰撞频率
            end
        end
        nu_Li1_Li2=nu_Li1+nu_Li2;             %构造数列的第二项，第一项是nu_Li1
        nu_Li1_Li2_Li3=nu_Li1+nu_Li2+nu_Li3;  %构造数列的第三项
        for i=1:nz-1               %遍历每个单元格
            for j=1:N_e_mcc(i)     %遍历每个单元格内抽出的电子
                R=rand();          %R为（0，1）之间的随机数
                if R<=nu_Li1(j,i)&&Ek_e(j,i)>=E_i_Li1                                        %Li0->Li+事件发生
                    N_Li1=N_Li1+1;                                                           %产生了一个新的Li+
                    N_e=N_e+1;                                                               %产生了一个新的电子
                    pos_Li1(N_Li1)=pos_e(distri_e(N_e_sam(j,i),i));                              %新产生的Li+的位置与入射电子的相同
                    vel_Li1(N_Li1,:)=randraw('maxwell',k*T/m_Li,1)*random_unit_vector(3,1)'; %新产生的Li+的速度服从麦克斯韦分布
                    pos_e(N_e)=pos_Li1(N_Li1);                                               %新产生的电子和新产生的Li+在相同的位置
                    vel_e_temp=0.5*sqrt(2*q_e*(Ek_e(j,i)-E_i_Li1)/m_e)*random_unit_vector(3,2);  %认为电子对离子的能量没有贡献，电子动能减去电离能后平分给散射电子和新产生的电子，random_unit_vector函数用来产生各向同性的电子
                    vel_e(N_e,:)=vel_e_temp(:,1)';                                           %新生成的电子的速度
                    vel_e(distri_e(N_e_sam(j,i),i),:)=vel_e_temp(:,2)';                        %散射电子的速度
                end
                if R>nu_Li1(j,i)&&R<=nu_Li1_Li2(j,i)&&Ek_e(j,i)>=E_i_Li2                         %Li+->Li2+事件发生
                    N_Li2=N_Li2+1;                                                          %产生了一个新的Li2+
                    N_e=N_e+1;                                                              %产生了一个新的电子
                    O_Li1=distri_Li1(randperm(count_Li(i,1),1),i);                          %消亡的Li+对应的序号，随机在第i个单元格中选择
                    N_Li1=N_Li1-1;                                                          %Li+数目减少一个
                    pos_Li2(N_Li2)=pos_Li1(O_Li1);                                          %新产生的Li2+的位置与消亡的Li+的位置相同
                    vel_Li2(N_Li2,:)=vel_Li1(O_Li1,:);                                      %新产生的Li2+的速度与消亡掉的Li+的速度相同
                    pos_e(N_e)=pos_Li2(N_Li2);                                              %新产生的电子的位置和Li2+产生的位置相同
                    vel_e_temp=0.5*sqrt(2*q_e*(Ek_e(j,i)-E_i_Li2)/m_e)*random_unit_vector(3,2); %认为电子对离子的能量没有贡献，电子动能减去电离能后平分给散射电子和新产生的电子，random_unit_vector函数用来产生各向同性的电子
                    vel_e(N_e,:)=vel_e_temp(:,1)';                                          %新生成的电子的速度
                    vel_e(distri_e(N_e_sam(j,i),i),:)=vel_e_temp(:,2)';                       %散射电子的速度
                end
                if R>nu_Li1_Li2(j,i)&&R<=nu_Li1_Li2_Li3(j,i)&&Ek_e(j,i)>=E_i_Li3                 %Li2+->Li3+事件发生
                    N_Li3=N_Li3+1;                                                          %产生了一个新的Li3+
                    N_e=N_e+1;                                                              %产生了一个新的电子
                    O_Li2=distri_Li2(randperm(count_Li(i,2),1),i);                          %消亡的Li2+对应的序号，随机在第i个单元格中选择
                    N_Li2=N_Li2-1;                                                          %Li2+数目减少一个
                    pos_Li3(N_Li3)=pos_Li2(O_Li2);                                          %新产生的Li3+的位置与消亡的Li2+的位置相同
                    vel_Li3(N_Li3,:)=vel_Li2(O_Li2,:);                                      %新产生的Li3+的速度与消亡掉的Li2+的速度相同
                    pos_e(N_e)=pos_Li3(N_Li3);                                              %新产生的电子的位置和Li3+产生的位置相同
                    vel_e_temp=0.5*sqrt(2*q_e*(Ek_e(j,i)-E_i_Li3)/m_e)*random_unit_vector(3,2); %认为电子对离子的能量没有贡献，电子动能减去电离能后平分给散射电子和新产生的电子，random_unit_vector函数用来产生各向同性的电子
                    vel_e(N_e,:)=vel_e_temp(:,1)';                                          %新生成的电子的速度
                    vel_e(distri_e(N_e_sam(j,i),i),:)=vel_e_temp(:,2)';                       %散射电子的速度
                end
            end
        end
        
        %分配Li+的电荷和电流密度到临近格点上
        for p=1:N_Li1
            
            fi=1+pos_Li1(p)/dz;                 %实际格点位置，为浮点数
            i=floor(fi);                        %对应整数格点位置
            hz=fi-i;                            %粒子与第i个格点之间的距离占距离步长的比例
            
            %将Li+电荷按照比例分配到临近的两个格点上
            chg(i,2)=chg(i,2)+(1-hz);      %分配到第i个格点上的电荷比例
            chg(i+1,2)=chg(i+1,2)+hz;      %分配到第i+1个格点上的电荷比例
            
            count_Li(i,1)=count_Li(i,1)+1;       %对每个单元格内的Li+计数
            
            %将速度按比例分配到临近的两个格点上
            den_vel_Li1(i,:)=den_vel_Li1(i,:)+spwt*q_e*vel_Li1(p,:)*(1-hz);   %分配到第i个格点上的速度
            den_vel_Li1(i+1,:)=den_vel_Li1(i,:)+spwt*q_e*vel_Li1(p,:)*hz;     %分配到第i+1个格点上的速度
        end
        
        %分配Li2+的电荷和电流密度到临近格点上
        for p=1:N_Li2
            
            fi=1+pos_Li2(p)/dz;                 %实际格点位置，为浮点数
            i=floor(fi);                        %对应整数格点位置
            hz=fi-i;                            %粒子与第i个格点之间的距离占距离步长的比例
            
            %将Li2+电荷按照比例分配到临近的两个格点上
            chg(i,3)=chg(i,3)+(1-hz);      %分配到第i个格点上的电荷比例
            chg(i+1,3)=chg(i+1,3)+hz;      %分配到第i+1个格点上的电荷比例
            
            count_Li(i,2)=count_Li(i,2)+1;       %对每个单元格内的Li2+计数
            
            %将速度按比例分配到临近的两个格点上
            den_vel_Li2(i,:)=den_vel_Li2(i,:)+spwt*2*q_e*vel_Li2(p,:)*(1-hz);   %分配到第i个格点上的速度
            den_vel_Li2(i+1,:)=den_vel_Li2(i,:)+spwt*2*q_e*vel_Li2(p,:)*hz;     %分配到第i+1个格点上的速度
        end
        
        %分配Li3+的电荷和电流密度到临近格点上
        for p=1:N_Li3
            
            fi=1+pos_Li3(p)/dz;                 %实际格点位置，为浮点数
            i=floor(fi);                        %对应整数格点位置
            hz=fi-i;                            %粒子与第i个格点之间的距离占距离步长的比例
            
            %将Li3+电荷按照比例分配到临近的两个格点上
            chg(i,4)=chg(i,4)+(1-hz);      %分配到第i个格点上的电荷比例
            chg(i+1,4)=chg(i+1,4)+hz;      %分配到第i+1个格点上的电荷比例
            
            count_Li(i,3)=count_Li(i,3)+1;       %对每个单元格内的Li3+计数
            
            %将速度按比例分配到临近的两个格点上
            den_vel_Li3(i,:)=den_vel_Li3(i,:)+spwt*3*q_e*vel_Li3(p,:)*(1-hz);   %分配到第i个格点上的速度
            den_vel_Li3(i+1,:)=den_vel_Li3(i,:)+spwt*3*q_e*vel_Li3(p,:)*hz;     %分配到第i+1个格点上的速度
        end
        
        %计算每个格点上的粒子密度
        for count_i=2:nz-1
            count_g_Li(count_i,:)=0.5*(count_Li(count_i-1,:)+count_Li(count_i,:))/(2*dz);
        end
        count_g_Li(1,:)=count_Li(1,:)/dz;      %第一个格点的粒子密度
        count_g_Li(nz,:)=count_Li(nz-1,:)/dz;  %最后一个格点的粒子密度
        
        J_i=1e5*(den_vel_Li1.*count_g_Li(:,1)+den_vel_Li2.*count_g_Li(:,2)+den_vel_Li3.*count_g_Li(:,3));  %离子在每个格点上的电流密度
        
        
        
        
        %计算电子的电荷密度和电流密度
        chg(:,1)=zeros(nz,1);                                                       %用于存放每个格点分配到的电子的电荷量
        %den=zeros(nz,1);                                                            %每个格点的电荷密度，单位：C/m
        den_vel_e=zeros(nz,3);                                                        %用于存放每个格点上的电流密度分量
        B_mic_old=B_mic;                                                            %B_mic_old用于和B_mic求平均，从而用于求解速度
        pos_e_half=zeros(N_e,1);                                                      %用于存放半时刻时的粒子位置
        
        for p_half=1:N_e                                                            %将所有粒子的位置向前移动到半时刻处
            pos_e_half(p_half,1)=pos_e(p_half,1)-vel_e(p_half,3)*0.5*dt_e;
            if pos_e_half(p_half,1)<0                                                 %两个if用来表示周期性边界条件
                pos_e_half(p_half,1)=pos_e_half(p_half,1)+z_max;
            end
            if pos_e_half(p_half,1)>z_max
                pos_e_half(p_half,1)=pos_e_half(p_half,1)-z_max;
            end
        end
        % pos_half=pos-vel(:,3)*0.5*dt;                                               %将电子的位置向前移动半个时间步长，以求得半时刻处的电流密度
        
        %首先将所有电子的电荷按照最近格点原则分配
        
        count=zeros(nz-1,1);            %用于计算每个单元格内的粒子数
        count_g=zeros(nz,1);            %用于计算每个格点上的粒子数密度
        
        %分配电子的电荷与电流密度到临近格点上
        for p=1:N_e
            fi=1+pos_e(p)/dz;            %实际格点位置，为浮点数
            i=floor(fi);               %对应整数格点位置
            hz=fi-i;                   %粒子与第i个格点之间的距离占距离步长的比例
            
            %将电子电荷按照比例分配到临近的两个格点上
            chg(i,1)=chg(i,1)+(1-hz);      %分配到第i个格点上的电荷比例
            chg(i+1,1)=chg(i+1,1)+hz;      %分配到第i+1个格点上的电荷比例
            
            fi_half=1+pos_e_half(p)/dz;  %半时刻处的格点位置
            i_half=floor(fi_half);     %半时刻处格点位置对应的整数格点位置
            hz_half=fi_half-i_half;    %半时刻处粒子与第i_half个格点之间的距离占距离步长的比例
            
            count(i_half)=count(i_half)+1;  %统计每个单元格内的粒子数目
            
            
            
            %将速度按比例分配到临近的两个格点上
            den_vel_e(i_half,:)=den_vel_e(i_half,:)+spwt*(-q_e)*vel_e(p,:)*(1-hz_half);  %分配到第i个格点上的速度
            den_vel_e(i_half+1,:)=den_vel_e(i_half+1,:)+spwt*(-q_e)*vel_e(p,:)*hz_half;  %分配到第i+1个格点上的速度
        end
        
        den=spwt*q_e*(-1*chg(:,1)+1*chg(:,2)+2*chg(:,3)+3*chg(:,4))/dz;       %每个格点上的电荷密度
        %J=spwt*q_e*(-1)*den_vel;  %每个格点上的电流密度
        %J=den_vel.*den;
        %J=1e9*den_vel;
        
        %计算每个格点上的电子密度
        for count_i=2:nz-1
            count_g(count_i)=0.5*(count(count_i-1)+count(count_i))/(2*dz);
        end
        count_g(1)=count(1)/dz;      %第一个格点的电子密度
        count_g(nz)=count(nz-1)/dz;  %最后一个格点的电子密度
        
        J=1e5*den_vel_e.*count_g+J_i;          %计算每个格点上的电流密度
        
        %应用边界条件，在边界上仅有一半的长度有贡献，因此边界上的密度应当×2
        den(1)=2*den(1);   %左边界
        den(nz)=2*den(nz); %右边界
        
        %开始求解泊松方程，得到格点上的电势，所用方法为高斯-赛德尔迭代法（G-S）
        % phi_tem=phi;
        
        
        %         phi=G_S_SOR(phi,den,tol,sor,N_phi,dz,nz);
        
        
        
        %         for j=1:N_phi     %迭代步数
        %
        %
        %             for i=2:nz-1  %中间格点的电势
        %                 phi_tem(i)=0.5*(den(i)/EPS0*dz*dz+phi(i-1)+phi(i+1));   %根据有限差分方法得到
        %             end
        %
        %             phi_tem(1)=phi(2);     %入口处为第一类边界条件
        %             phi_tem(nz)=phi(nz-1); %出口处也为第一类边界条件
        %
        %             if mod(j,10)==0
        %                 R=norm(phi_tem-phi); %求出本次迭代的残差
        %                 if (R<=tol)          %收敛条件
        %                     phi=phi_tem;
        %                     fprintf('Congratulations! The solver of the electric potential is converged at %d step\n',j);
        %                     break;
        %                 end
        %             end
        %             phi=phi_tem;
        %         end
        %         if j>=N_phi                  %如果电势求解不收敛，则终止程序的运行
        %             fprintf('The result of the electric potential is not convergent, the program is stopped!\n');
        %             break;
        %         end
        %
        %
        
        
        
        A=blktridiag(2,-1,-1,nz);
        B=full(A);
        B(1,:)=0*B(1,:);
        B(1,1)=1;
        B(nz,:)=0*B(nz,:);
        B(nz,nz)=1;
        den(1)=0;
        den(nz)=0;
        den=den/EPS0*(dz)^2;
        % [phi,k]=G_S(B,den,phi,1e-3,2000);
        K=B;
        %F=b;
        F=den;
        W=zeros(size(B,1),1);
        level=3;
        node_order=[-1 1 1];
        relax_it=2;
        relax_para=1;
        post_smoothing=1;
        max_iter=20000;
        % tol=1e-06;
        pc_type=2;
        connection_threshold=0.25;
        [phi,error,iter,flag]=AMG(K,F, W, level,  relax_it, relax_para, ...
            post_smoothing, max_iter, tol,  pc_type, connection_threshold);
        
        
        
        
        
        %求解每个格点上的静电场，由于是一维，因此只能得到z方向的电场强度
        E_sta(2:nz-1)=(phi(1:nz-2)-phi(3:nz))/2/dz;      %中间格点的静电场
        E_sta(1)=(phi(1)-phi(2))/dz;                     %入口处静电场
        E_sta(nz)=(phi(nz-1)-phi(nz))/dz;                %出口处静电场
        
        %计算微波电磁场的各个分量
        E_mic_xold=E0*cos(w*1e9*2*pi()*(ts_e-2)*dt_e);
        E_mic_yold=E0*sin(w*1e9*2*pi()*(ts_e-2)*dt_e);
        E_mic(1,1)=E0*cos(w*1e9*2*pi()*(ts_e-1)*dt_e);              %入口处x方向的电场强度E_x=E0cos(wt)
        % B_mic_xin_full=B0*sin(w*1e9*2*pi()*(ts-0.5)*dt);       %入口处x方向磁感应强度，时间为整
        B_mic_xin=B0*sin(w*1e9*2*pi()*(ts_e-0.5)*dt_e);             %入口处x方向的磁感应强度为B_x=B0sin(wt)，且时间上向后移动Δt/2
        E_mic(1,2)=E0*sin(w*1e9*2*pi()*(ts_e-1)*dt_e);              %入口处y方向的电场强度E_y=E0sin(wt)
        %B_mic_yin_full=B0*cos(w*1e9*2*pi()*(ts-0.5)*dt);       %入口处y方向磁感应强度，时间为整
        B_mic_yin=B0*cos(w*1e9*2*pi()*(ts_e-0.5)*dt_e);             %入口处y方向的磁感应强度B_y=B0cos(wt)，且时间上向后移动Δt/2
        %B_mic(1,1)=B_mic_xin+dt/(2*dz)*(E_mic(2,2)-E_mic(1,2));                   %将Bx在空间和时间上都向后移动Δz/2和Δt/2
        B_mic(1,1)=B_mic_xin+dz/2*(1/(c^2)*(E_mic(1,2)-E_mic_yold)/dt_e+MU0*J(1,2)); %入口处的Bx，将Bx在空间上向后移动Δz/2，向后差分
        B_mic(1,2)=B_mic_yin-dz/2*(1/(c^2)*(E_mic(1,1)-E_mic_xold)/dt_e+MU0*J(1,1)); %入口处的By，将By在空间上向后移动Δz/2，向后差分
        %B_mic(1,2)=B_mic_yin-dt/(2*dz)*(E_mic(2,1)-E_mic(1,1));                   %将By在空间和时间上都向后移动Δz/2和Δt/2
        
        % E_mic(nz,1)=0; %出口处x方向的电场强度0
        % E_mic(nz,2)=0; %出口处y方向的电场强度0
        
        B_mic_xout=B_mic(nz,1);
        B_mic_yout=B_mic(nz,2);
        
        %B_mic_xout=0;  %出口处x方向的磁感应强度0
        % B_mic_yout=0; %出口处y方向的磁感应强度0
        % B_mic(nz-1,1)=B_mic_xout-dz/2*MU0*J(nz,2); %出口处的Bx，向前差分
        % B_mic(nz-1,2)=B_mic_yout+dz/2*MU0*J(nz,1); %出口处的By，向前差分
        
        B_mic(2:nz-1,1)=B_mic(2:nz-1,1)+dt_e/dz*(E_mic(3:nz,2)-E_mic(2:nz-1,2));   %Bx
        B_mic(2:nz-1,2)=B_mic(2:nz-1,2)-dt_e/dz*(E_mic(3:nz,1)-E_mic(2:nz-1,1));   %By
        E_mic(2:nz-1,1)=E_mic(2:nz-1,1)-c^2*dt_e*(MU0*J(2:nz-1,1)+(B_mic(2:nz-1,2)-B_mic(1:nz-2,2))/dz);    %Ex
        E_mic(2:nz-1,2)=E_mic(2:nz-1,2)+c^2*dt_e*(-1*MU0*J(2:nz-1,2)+(B_mic(2:nz-1,1)-B_mic(1:nz-2,1))/dz); %Ey
        %E_mic(1:nz,3)=E_mic(1:nz,3)-c^2*MU0*dt*J(1:nz,3); %Ez
        
        E_boun(ts_e,:)=E_mic(nz-1,:);
        B_boun(ts_e,:)=B_mic(nz-1,:);
        if ts_e>2
            E_mic(nz,:)=E_boun(ts_e-2,:);
            B_mic(nz,:)=B_boun(ts_e-2,:);
        end
        
        %将格点上的静电场和电磁场全部分配到每个电子上，另外还要加上外部永磁体产生的磁场
        E_e=zeros(max_part,3);           %每个电子所在位置处的电场强度
        B_e=zeros(max_part,3);           %每个电子所在位置处的磁感应强度
        %B_e(:,3)=B_e(:,3)+B_ex_z;        %z方向的磁感应强度需要加上外部磁感应强度
        B_mic_ave=0.5*(B_mic_old+B_mic); %利用平均，将两个半时间的值平均到到整时间上
        p=1;                             %用于遍历所有粒子
        while p<=N_e
            fi=1+pos_e(p)/dz;  %实际格点位置，为浮点数
            i=floor(fi);     %对应整数格点位置
            hz=fi-i;         %粒子与第i个格点之间的距离占距离步长的比例
            E_e(p,3)=E_e(p,3)+E_sta(i)*(1-hz)+E_sta(i+1)*hz;     %分配粒子之间的静电场
            E_e(p,:)=E_e(p,:)+E_mic(i,:)*(1-hz)+E_mic(i+1,:)*hz; %分配微波产生的电场
            %B_e(p,3)=B_e(p,3)+B_ex_z(i)*(1-hz)+B_ex_z(i+1)*hz;   %将格点上的静磁场分配给每个粒子
            B_e(p,3)=1e-4*(2595.7+32.4092*(1000*pos_e(p))-3.1038*(1000*pos_e(p))^2-...   %插值得到粒子位置处的磁感应强度
                0.0130847*(1000*pos_e(p))^3+0.00198809*(1000*pos_e(p))^4-2.48996e-5*(1000*pos_e(p))^5+8.9815e-8*(1000*pos_e(p))^6);
            if fi<=1.5                                           %四个if用来确定粒子位置，从而为之分配微波产生的磁感应强度
                B_e(p,:)=B_e(p,:)+2*hz*B_mic_ave(1,:)+(1-2*hz)*[B_mic_xin,B_mic_yin,0];
            end
            if fi>=nz-0.5
                B_e(p,:)=B_e(p,:)+2*(1-hz)*B_mic_ave(nz-1,:)+(2*hz-1)*[B_mic_xout,B_mic_yout,0];
            end
            if fi>1.5 && fi<nz-0.5 && hz<=0.5
                B_e(p,:)=B_e(p,:)+(hz+0.5)*B_mic_ave(i,:)+(0.5-hz)*B_mic_ave(i+1,:);
            end
            if fi>1.5 && fi<nz-0.5 && hz>0.5
                B_e(p,:)=B_e(p,:)+(1.5-hz)*B_mic_ave(i,:)+(hz-0.5)*B_mic_ave(i+1,:);
            end
            
            %利用BORIS方法求出下一时刻的速度和位置
            vel_e(p,:)=UpdateVelocity(E_e(p,:),B_e(p,:),vel_e(p,:),dt_e); %更新速度
            pos_e(p,1)=pos_e(p,1)+vel_e(p,3)*dt_e;                          %更新位置
            
            %设置吸收边界条件，电子与壁相互作用后消亡掉
            if pos_e(p,1)<=0             %左边界
                pos_e(p,1)=pos_e(N_e,1); %将最后一个电子的位置替换到消亡位置处
                vel_e(p,:)=vel_e(N_e,:); %将最后一个电子的速度替换到消亡位置处
                N_e=N_e-1;               %电子数减少一个
            end
            if pos_e(p,1)>z_max          %右边界
                pos_e(p,1)=pos_e(N_e,1); %将最后一个电子的位置替换到消亡位置处
                vel_e(p,:)=vel_e(N_e,:); %将最后一个电子的速度替换到消亡位置处
                N_e=N_e-1;               %电子数减少一个
            end
            
            p=p+1;
        end
        
        fprintf('ts_e: %d, ts_i: %d,N_Li1: %d, N_Li2: %d, N_Li3: %d, P_emax: %f\n Average Energy: %f eV\n', ts_e,ts_i,N_Li1,N_Li2,N_Li3,max(P_emax),mean(Ek_e(1,:)));
        % plot(1:nz,E_mic(:,1));
        % pause(0.1);
    end
    
    %将格点上的静电场和电磁场全部分配到每个离子上，另外还要加上外部永磁体产生的磁场
    E_Li1=zeros(max_part,3);           %每个Li+所在位置处的电场强度
    E_Li2=zeros(max_part,3);          %每个Li2+所在位置处的电场强度
    E_Li3=zeros(max_part,3);          %每个Li3+所在位置处的电场强度
    B_Li1=zeros(max_part,3);           %每个Li+所在位置处的磁感应强度
    B_Li2=zeros(max_part,3);          %每个Li2+所在位置处的磁感应强度
    B_Li3=zeros(max_part,3);          %每个Li3+所在位置处的磁感应强度
    %B_e(:,3)=B_e(:,3)+B_ex_z;        %z方向的磁感应强度需要加上外部磁感应强度
    B_mic_ave=0.5*(B_mic_old+B_mic); %利用平均，将两个半时间的值平均到到整时间上
    %     N_i=[N_Li1,N_Li2,N_Li3];
    %     for i_Li=1:3
    
    p=1;                             %用于遍历所有Li+
    while p<=N_Li1
        fi=1+pos_e(p)/dz;  %实际格点位置，为浮点数
        i=floor(fi);     %对应整数格点位置
        hz=fi-i;         %粒子与第i个格点之间的距离占距离步长的比例
        %             i_str=num2str(i);
        %             p_str=num2str(p);
        %             hz_str=num2str(hz);
        %             eval(['E_Li',i_str,'(',p_str,',',num2str(3),')','=','E_Li',i_str,'(',p_str,',',num2str(3),')','+E_sta(',i_str,')*(',num2str(1-hz),'+E_sta(',num2str(i+1),')*',num2str(hz)]);     %分配粒子之间的静电场
        %eval(['E_Li',num2str(i),'(',num2str(p),',',':',')'])=eval(['E_Li',num2str(i),'(',num2str(p),',',':',')'])+E_mic(i,:)*(1-hz)+E_mic(i+1,:)*hz; %分配微波产生的电场
        %             eval(['E_Li',i_str,'(',p_str,',:',')','=','E_Li',i_str,'(',p_str,',:',')','+E_mic(',i_str,',:',')*(',num2str(1-hz),'+E_mic(',num2str(i+1),',:)','*',num2str(hz)]);
        E_Li1(p,3)=E_Li1(p,3)+E_sta(i)*(1-hz)+E_sta(i+1)*hz;     %分配粒子之间的静电场
        E_Li1(p,:)=E_Li1(p,:)+E_mic(i,:)*(1-hz)+E_mic(i+1,:)*hz; %分配微波产生的电场
        %B_e(p,3)=B_e(p,3)+B_ex_z(i)*(1-hz)+B_ex_z(i+1)*hz;   %将格点上的静磁场分配给每个粒子
        %             eval(['B_Li',i_str,'(',p_str,',3)=1e-4*(2595.7+32.4092*(1000*pos_Li',i_str,'(',p_str,')',')-3.1038*(1000*pos_Li',i_str,'(',p_str,'))^2-0.0130847*(1000*pos_Li'...
        %                 ,i_str,'(',p_str,'))^3+0.00198809*(1000*pos_Li',i_str,'(',p_str,'))^4-2.48996e-5*(1000*pos_Li',i_str,'(',p_str,'))^5+8.9815e-8*(1000*pos_Li',...
        %                 i_str,'(',p_str,'))^6)']);
        B_Li1(p,3)=1e-4*(2595.7+32.4092*(1000*pos_Li1(p))-3.1038*(1000*pos_Li1(p))^2-...   %插值得到粒子位置处的磁感应强度
            0.0130847*(1000*pos_Li1(p))^3+0.00198809*(1000*pos_Li1(p))^4-2.48996e-5*(1000*pos_Li1(p))^5+8.9815e-8*(1000*pos_Li1(p))^6);
        if fi<=1.5                                           %四个if用来确定粒子位置，从而为之分配微波产生的磁感应强度
            B_Li1(p,:)=B_Li1(p,:)+2*hz*B_mic_ave(1,:)+(1-2*hz)*[B_mic_xin,B_mic_yin,0];
            %                 eval(['B_Li',i_str,'(',p_str,':)=B_Li',i_str,'(',p_str,':)+2*',hz_str,'*B_mic_ave(1,:)+(1-2*',hz_str,')*[B_mic_xin,B_mic_yin,0]']);
        end
        if fi>=nz-0.5
            B_Li1(p,:)=B_Li1(p,:)+2*(1-hz)*B_mic_ave(nz-1,:)+(2*hz-1)*[B_mic_xout,B_mic_yout,0];
            %                 eval(['B_Li',i_str,'(',p_str,':)=B_Li',i_str,'(',p_str,':)+2*(1-',hz_str,')*B_mic_ave(',num2str(nz),'-1,:)+(2*',hz_str,'-1)*[B_mic_xout,B_mic_yout,0]']);
        end
        if fi>1.5 && fi<nz-0.5 && hz<=0.5
            B_Li1(p,:)=B_Li1(p,:)+(hz+0.5)*B_mic_ave(i,:)+(0.5-hz)*B_mic_ave(i+1,:);
        end
        if fi>1.5 && fi<nz-0.5 && hz>0.5
            B_Li1(p,:)=B_Li1(p,:)+(1.5-hz)*B_mic_ave(i,:)+(hz-0.5)*B_mic_ave(i+1,:);
        end
        
        %利用BORIS方法求出下一时刻的速度和位置
        vel_Li1(p,:)=UpdateVelocity(E_Li1(p,:),B_Li1(p,:),vel_Li1(p,:),dt_i); %更新速度
        pos_Li1(p,1)=pos_Li1(p,1)+vel_Li1(p,3)*dt_i;                          %更新位置
        
        %设置吸收边界条件，Li+与壁相互作用后消亡掉
        if pos_Li1(p,1)<0             %左边界
            pos_Li1(p,1)=pos_Li1(N_Li1,1); %将最后一个Li+的位置替换到消亡位置处
            vel_Li1(p,:)=vel_Li1(N_Li1,:); %将最后一个Li+的速度替换到消亡位置处
            N_Li1=N_Li1-1;               %Li+数减少一个
        end
        if pos_Li1(p,1)>z_max          %右边界
            pos_Li1(p,1)=pos_Li1(N_Li1,1); %将最后一个Li+的位置替换到消亡位置处
            vel_Li1(p,:)=vel_Li1(N_Li1,:); %将最后一个Li+的速度替换到消亡位置处
            N_Li1=N_Li1-1;               %Li+数减少一个
        end
        
        p=p+1;
    end
    %     end
    p=1;                             %用于遍历所有Li2+
    while p<=N_Li2
        fi=1+pos_e(p)/dz;  %实际格点位置，为浮点数
        i=floor(fi);     %对应整数格点位置
        hz=fi-i;         %粒子与第i个格点之间的距离占距离步长的比例
        
        E_Li2(p,3)=E_Li2(p,3)+E_sta(i)*(1-hz)+E_sta(i+1)*hz;     %分配粒子之间的静电场
        E_Li2(p,:)=E_Li2(p,:)+E_mic(i,:)*(1-hz)+E_mic(i+1,:)*hz; %分配微波产生的电场
        
        B_Li2(p,3)=1e-4*(2595.7+32.4092*(1000*pos_Li1(p))-3.1038*(1000*pos_Li1(p))^2-...   %插值得到粒子位置处的磁感应强度
            0.0130847*(1000*pos_Li1(p))^3+0.00198809*(1000*pos_Li1(p))^4-2.48996e-5*(1000*pos_Li1(p))^5+8.9815e-8*(1000*pos_Li1(p))^6);
        if fi<=1.5                                           %四个if用来确定粒子位置，从而为之分配微波产生的磁感应强度
            B_Li2(p,:)=B_Li2(p,:)+2*hz*B_mic_ave(1,:)+(1-2*hz)*[B_mic_xin,B_mic_yin,0];
            
        end
        if fi>=nz-0.5
            B_Li2(p,:)=B_Li2(p,:)+2*(1-hz)*B_mic_ave(nz-1,:)+(2*hz-1)*[B_mic_xout,B_mic_yout,0];
            
        end
        if fi>1.5 && fi<nz-0.5 && hz<=0.5
            B_Li2(p,:)=B_Li2(p,:)+(hz+0.5)*B_mic_ave(i,:)+(0.5-hz)*B_mic_ave(i+1,:);
        end
        if fi>1.5 && fi<nz-0.5 && hz>0.5
            B_Li2(p,:)=B_Li2(p,:)+(1.5-hz)*B_mic_ave(i,:)+(hz-0.5)*B_mic_ave(i+1,:);
        end
        
        %利用BORIS方法求出下一时刻的速度和位置
        vel_Li2(p,:)=UpdateVelocity(E_Li2(p,:),B_Li2(p,:),vel_Li2(p,:),dt_i); %更新速度
        pos_Li2(p,1)=pos_Li2(p,1)+vel_Li2(p,3)*dt_i;                          %更新位置
        
        %设置吸收边界条件，Li2+与壁相互作用后消亡掉
        if pos_Li2(p,1)<0             %左边界
            pos_Li2(p,1)=pos_Li2(N_Li2,1); %将最后一个Li2+的位置替换到消亡位置处
            vel_Li2(p,:)=vel_Li2(N_Li2,:); %将最后一个Li2+的速度替换到消亡位置处
            N_Li2=N_Li2-1;               %Li2+数减少一个
        end
        if pos_Li2(p,1)>z_max          %右边界
            pos_Li2(p,1)=pos_Li2(N_Li2,1); %将最后一个Li2+的位置替换到消亡位置处
            vel_Li2(p,:)=vel_Li2(N_Li2,:); %将最后一个Li2+的速度替换到消亡位置处
            N_Li2=N_Li2-1;               %Li2+数减少一个
        end
        
        p=p+1;
    end
    
    
    p=1;                             %用于遍历所有Li3+
    while p<=N_Li3
        fi=1+pos_e(p)/dz;  %实际格点位置，为浮点数
        i=floor(fi);     %对应整数格点位置
        hz=fi-i;         %粒子与第i个格点之间的距离占距离步长的比例
        
        E_Li3(p,3)=E_Li3(p,3)+E_sta(i)*(1-hz)+E_sta(i+1)*hz;     %分配粒子之间的静电场
        E_Li3(p,:)=E_Li3(p,:)+E_mic(i,:)*(1-hz)+E_mic(i+1,:)*hz; %分配微波产生的电场
        
        B_Li3(p,3)=1e-4*(2595.7+32.4092*(1000*pos_Li1(p))-3.1038*(1000*pos_Li1(p))^2-...   %插值得到粒子位置处的磁感应强度
            0.0130847*(1000*pos_Li1(p))^3+0.00198809*(1000*pos_Li1(p))^4-2.48996e-5*(1000*pos_Li1(p))^5+8.9815e-8*(1000*pos_Li1(p))^6);
        if fi<=1.5                                           %四个if用来确定粒子位置，从而为之分配微波产生的磁感应强度
            B_Li3(p,:)=B_Li3(p,:)+2*hz*B_mic_ave(1,:)+(1-2*hz)*[B_mic_xin,B_mic_yin,0];
            
        end
        if fi>=nz-0.5
            B_Li3(p,:)=B_Li3(p,:)+2*(1-hz)*B_mic_ave(nz-1,:)+(2*hz-1)*[B_mic_xout,B_mic_yout,0];
            
        end
        if fi>1.5 && fi<nz-0.5 && hz<=0.5
            B_Li3(p,:)=B_Li3(p,:)+(hz+0.5)*B_mic_ave(i,:)+(0.5-hz)*B_mic_ave(i+1,:);
        end
        if fi>1.5 && fi<nz-0.5 && hz>0.5
            B_Li3(p,:)=B_Li3(p,:)+(1.5-hz)*B_mic_ave(i,:)+(hz-0.5)*B_mic_ave(i+1,:);
        end
        
        %利用BORIS方法求出下一时刻的速度和位置
        vel_Li3(p,:)=UpdateVelocity(E_Li3(p,:),B_Li3(p,:),vel_Li3(p,:),dt_i); %更新速度
        pos_Li3(p,1)=pos_Li3(p,1)+vel_Li3(p,3)*dt_i;                          %更新位置
        
        %设置吸收边界条件，Li3+与壁相互作用后消亡掉
        if pos_Li3(p,1)<0             %左边界
            pos_Li3(p,1)=pos_Li3(N_Li3,1); %将最后一个Li3+的位置替换到消亡位置处
            vel_Li3(p,:)=vel_Li3(N_Li3,:); %将最后一个Li3+的速度替换到消亡位置处
            N_Li3=N_Li3-1;               %Li3+数减少一个
        end
        if pos_Li3(p,1)>z_max          %右边界
            pos_Li3(p,1)=pos_Li2(N_Li3,1); %将最后一个Li3+的位置替换到消亡位置处
            vel_Li2(p,:)=vel_Li2(N_Li3,:); %将最后一个Li3+的速度替换到消亡位置处
            N_Li3=N_Li3-1;               %Li3+数减少一个
        end
        
        p=p+1;
    end
end

filename=strcat(datestr(clock,'yy-mm-dd-HH-MM-SS'),'E_mic.txt');
fid=fopen(filename,'wt');       %将微波电场数据写到E_mic.txt中，三列分别是时Ex,Ey,Ez，行数代表的是nz个格点
[row,col]=size(E_mic);
for i=1:1:row
    for j=1:1:col
        if (j==col)
            fprintf(fid,'%g\n',E_mic(i,j));
        else
            fprintf(fid,'%g\t',E_mic(i,j));
        end
    end
end
fclose(fid);

filename=strcat(datestr(clock,'yy-mm-dd-HH-MM-SS'),'B_mic.txt');
fid=fopen(filename,'wt');       %将微波磁场数据写到B_mic.txt中，三列分别是时Bx,By,Bz，行数代表的是nz个格点
[row,col]=size(B_mic);
for i=1:1:row
    for j=1:1:col
        if (j==col)
            fprintf(fid,'%g\n',B_mic(i,j));
        else
            fprintf(fid,'%g\t',B_mic(i,j));
        end
    end
end
fclose(fid);

filename=strcat(datestr(clock,'yy-mm-dd-HH-MM-SS'),'B_ex_z.txt');
fid=fopen(filename,'wt');       %将外界永磁体产生的磁场写入B_ex_z.txt文件中，第一列为每个格点的z坐标，第二列为对应的外界磁感应强度
%[row,col]=[nz-1,2];
for i=1:1:nz
    for j=1:1:2
        if (j==2)
            fprintf(fid,'%g\n',B_ex_z(i));
        else
            fprintf(fid,'%f\t',(i-1)*dz);
        end
    end
end
fclose(fid);

filename=strcat(datestr(clock,'yy-mm-dd-HH-MM-SS'),'z-v phase.txt');
fid=fopen(filename,'wt');
for i=1:N_e
    for j=1:4
        if j==1
            fprintf(fid,'%f\t',pos_e(i));
        end
        if j>1&&j<4
            fprintf(fid,'%g\t',vel_e(i,j-1));
        end
        if j==4
            fprintf(fid,'%g\n',vel_e(i,j-1));
        end
    end
end
fclose(fid);


fprintf('finished!\n');

%输出微波电磁场分布的图像
figure(1);
subplot(2,3,1);
plot(1:nz,E_mic(:,1));
xlabel('grid');
ylabel('E_x');
subplot(2,3,2);
plot(1:nz,E_mic(:,2));
xlabel('grid');
ylabel('E_y');
subplot(2,3,3);
plot(1:nz,E_mic(:,3));
xlabel('grid');
ylabel('E_z');
subplot(2,3,4);
plot(1:nz,B_mic(:,1));
xlabel('grid');
ylabel('B_x');
subplot(2,3,5);
plot(1:nz,B_mic(:,2));
xlabel('grid');
ylabel('B_y');
subplot(2,3,6);
plot(1:nz,B_mic(:,3));
xlabel('grid');
ylabel('B_z');

%输出z-v图像
figure(2);
subplot(3,1,1);
scatter(pos_e(1:N_e),vel_e(1:N_e,1));
xlabel('z(m)');
ylabel('v_x(m/s)');
subplot(3,1,2);
scatter(pos_e(1:N_e),vel_e(1:N_e,2));
xlabel('z(m)');
ylabel('v_y(m/s)');
subplot(3,1,3);
scatter(pos_e(1:N_e),vel_e(1:N_e,3));
xlabel('z(m)');
ylabel('v_z(m/s)');


