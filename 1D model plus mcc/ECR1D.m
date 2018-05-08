%1 2/2D ECR plasma simulation using PIC-MCC
clear variables
clc

global q_e m_e EPS0

%�����������
q_e=1.602e-19;              %���ӵ�������λ��C
AMU = 1.661e-27;            %������������λ��kg
m_e=9.1e-31;                %������������λ��kg
m_Li=7*AMU;                 %Li-7��������λ��kg
EPS0 = 8.854e-12;           %����еĽ�糣������λ��F/m
MU0=4*pi()*1e-7;            %����еĴŵ��ʣ���λ��H/m
k=1.38e-23;                 %����������������λ��J/K
c=3e8;                      %���٣���λ��m/s
%c=3e2;                     %�����ٸ�С���Թ۲����ӶԵ�Ų���Ӱ�죨�����ʾ�����������⣩
Te_average=10;              %ƽ�������¶ȣ���λ��eV
n_e_average=5e17;           %ƽ�������ܶȣ���λ��m^-3
n0=1e18;                    %����ԭ���ܶȣ���λ��m^-3
T=273.15+300;               %ϵͳ�¶�Ϊ300�棬��λ��K
lamda=sqrt(EPS0*Te_average/(n0*q_e));    %�°ݳ��ȣ���λ��m
z_min=0;                                 %z�������Сֵ����λ��m
L_real=0.1035;                           %����Դz��������ֵ����λ��m
dz=lamda;                                %���ÿռ���
dt_e=dz/(2*c);                             %����ʱ������Ϊ��֤ģ�⾫ȷ�ȣ���֤dt<=dz/c��һ����Ϊdz/(2*c)��λ��s
dt_i=100*dt_e;                           %���ӵ��˶�ʱ������Ϊ�����˶�ʱ���100��
nz=round(L_real/dz)+1;                   %�����
z_max=dz*(nz-1);                         %ģ������z��������ֵ
spwt=1e8;                                %ÿ�������Ӱ�������ʵ������
vth_e=sqrt(2*q_e*0.5/m_e);               %���ӵ����ٶȣ�������ӵ�����Ϊ0.5 eV����λ��m/s
N_e=100000;                                %��ʼ������ĿΪN_e
N_Li1=0;                                 %Li+����Ŀ
N_Li2=0;                                 %Li2+����Ŀ
N_Li3=0;                                 %Li3+����Ŀ
N_Li_ex=0;                               %����̬Liԭ�ӵ���Ŀ
tol=1e-6;                               %�������ľ���ֵ
sor=1.4;                                 %���ɳ۵���������
max_part=1000000;                           %�������������ֵ
%B_ex_z=0.0875;                          %�ⲿ�����������ĴŸ�Ӧǿ�ȣ���λ��T
E0=1.5e8;                                %΢���糡ǿ�ȵķ�ֵ����λ��V/m����ȷ����
B0=0.0001;                               %΢���Ÿ�Ӧǿ�ȵķ�ֵ����λ:T����ȷ����
w=2.45;                                  %΢����Ƶ��,��λ:GHz
step_num=1000;                          %�ܹ����е�ʱ�䲽��
N_phi=5e5;                               %�������ĵ�������
R_i_max=7.8034e-7;                       %������̬�ܵ����ʵ����ֵ
E_i_Li1=5.4;                             %Li0->Li1�ĵ����ܣ���λ��eV
E_i_Li2=75.77;                           %Li1->Li2�ĵ����ܣ���λ��eV
E_i_Li3=122.664;                         %Li2->Li3�ĵ����ܣ���λ��eV
z=zeros(nz,1);                           %����ʢ��ÿ����������
z(2:nz,1)=(1:nz-1)'*dz*1000;                  %ÿ����������ֵ����λ��mm
B_ex_z=0.0001*(2595.7+32.4092*z-3.1038*z.^2-0.0130847*z.^3+0.00198809*z.^4-2.48996e-5*z.^5+8.9815e-8*z.^6);   %����ʵ���ϲ�õĴŸ�Ӧǿ��ֵ������ϣ���λ��T
%plot(z,B_ex_z);
B_ir=[-7.97027	4.27075	-3.96662	1.78936	-0.400521	0.0346181
    -12.2121	9.29911	-6.78284	2.46302	-0.453047	0.0332396
    -34.5328	36.724	-21.6441	6.55211	-1.0141	0.0636916]; %fitting parameters of ionization rate


%Ԥ����ռ�
phi=zeros(nz,1);               %����
phi_tem=zeros(nz,1);           %���������м����
B_mic=zeros(nz,3);             %΢���Ĵų�
E_mic=zeros(nz,3);             %΢���ĵ糡
E_sta=zeros(nz,1);             %������������������ľ��糡
E_e=zeros(max_part,3);         %����λ�ô��ĵ糡ǿ��
B_e=zeros(max_part,3);         %����λ�ô��ĴŸ�Ӧǿ��
E_boun=zeros(step_num,3);      %�����������ձ߽�
B_boun=zeros(step_num,3);      %�����������ձ߽�
%den=zeros(nz,1);              %ÿ������ϵ��������ܶ�
%chg=zeros(nz,1);              %ÿ������ϵĵ�ɱ���
J=zeros(nz,3);                 %ÿ������ϵĵ����ܶ�
pos_e=zeros(max_part,1);         %���ӵ�λ����Ϣ
vel_e=zeros(max_part,3);         %���ӵ��ٶ���Ϣ
pos_Li1=zeros(max_part,1);    %Li+��λ����Ϣ
vel_Li1=zeros(max_part,3);    %Li+���ٶ���Ϣ
pos_Li2=zeros(max_part,1);    %Li2+��λ����Ϣ
vel_Li2=zeros(max_part,3);    %Li2+���ٶ���Ϣ
pos_Li3=zeros(max_part,1);    %Li3+��λ����Ϣ
vel_Li3=zeros(max_part,3);    %Li3+���ٶ���Ϣ
pos_Li_ex=zeros(max_part,1);    %����̬��Liԭ�ӵ�λ����Ϣ
vel_Li_ex=zeros(max_part,3);    %����̬��Liԭ�ӵ��ٶ���Ϣ



%������ʼ���ӵ�λ�ú��ٶ�
pos_e(1:N_e)=rand(N_e,1)*z_max;                                                   %��ʼ���ӵ�λ����0��z_max֮��

%vel(1:N_e/2,:)=vth_e*rand(N_e/2,3);
%vel(N_e/2+1:N_e,:)=-vth_e*rand(N_e/2,3);

%vel(1:N_e,:)=1e8;
vel_e(1:N_e,:)=vth_e*2*(rand(N_e,3)+rand(N_e,3)+rand(N_e,3)-1.5);                 %���ó�ʼ���ӵ��ٶ������ٶȵ�-3����3��֮��
vel_e(1:N_e,:)=UpdateVelocity(E_e(1:N_e,:),B_e(1:N_e,:),vel_e(1:N_e,:),-0.5*dt_e);    %�����ӵ��ٶ���ǰ�ƶ����ʱ�䲽��

%��ʼ��ѭ��������ʱ���ѭ��
for ts_i=1:step_num                 %���е�ʱ�䲽��
    
    %chg=zeros(nz,4);                %���ڴ��ÿ�������䵽�ĵ������1~4�зֱ����:���ӣ�Li1+��Li2+��Li3+
    %den_vel_Li1=zeros(nz,3);        %���ڴ��ÿ�������Li+�ĵ����ܶȷ���
    %den_vel_Li2=zeros(nz,3);        %���ڴ��ÿ�������Li2+�ĵ����ܶȷ���
    %den_vel_Li3=zeros(nz,3);        %���ڴ��ÿ�������Li3+�ĵ����ܶȷ���
    %count_Li=zeros(nz-1,3);            %���ڼ���ÿ����Ԫ���ڵ���������1~3�зֱ����Li+��Li2+��Li3+
    %count_g_Li=zeros(nz,3);            %���ڼ���ÿ������ϵ��������ܶȣ�1~3�зֱ����Li+��Li2+��Li3+
    
    %����Li+�ĵ�ɺ͵����ܶȵ��ٽ������
    %for p=1:N_Li1
    
    %    fi=1+pos_Li1(p)/dz;                 %ʵ�ʸ��λ�ã�Ϊ������
    %    i=floor(fi);                        %��Ӧ�������λ��
    %   hz=fi-i;                            %�������i�����֮��ľ���ռ���벽���ı���
    
    %         %��Li+��ɰ��ձ������䵽�ٽ������������
    %         chg(i,2)=chg(i,2)+(1-hz);      %���䵽��i������ϵĵ�ɱ���
    %         chg(i+1,2)=chg(i+1,2)+hz;      %���䵽��i+1������ϵĵ�ɱ���
    %
    %         count_Li(i,1)=count_Li(i,1)+1;       %��ÿ����Ԫ���ڵ�Li+����
    %
    %         %���ٶȰ��������䵽�ٽ������������
    %         den_vel_Li1(i,:)=den_vel_Li1(i,:)+spwt*q_e*vel_Li1(p,:)*(1-hz);   %���䵽��i������ϵ��ٶ�
    %         den_vel_Li1(i+1,:)=den_vel_Li1(i,:)+spwt*q_e*vel_Li1(p,:)*hz;     %���䵽��i+1������ϵ��ٶ�
    %     end
    %
    %     %����Li2+�ĵ�ɺ͵����ܶȵ��ٽ������
    %     for p=1:N_Li2
    %
    %         fi=1+pos_Li2(p)/dz;                 %ʵ�ʸ��λ�ã�Ϊ������
    %         i=floor(fi);                        %��Ӧ�������λ��
    %         hz=fi-i;                            %�������i�����֮��ľ���ռ���벽���ı���
    %
    %         %��Li2+��ɰ��ձ������䵽�ٽ������������
    %         chg(i,3)=chg(i,3)+(1-hz);      %���䵽��i������ϵĵ�ɱ���
    %         chg(i+1,3)=chg(i+1,3)+hz;      %���䵽��i+1������ϵĵ�ɱ���
    %
    %         count_Li(i,2)=count_Li(i,2)+1;       %��ÿ����Ԫ���ڵ�Li2+����
    %
    %         %���ٶȰ��������䵽�ٽ������������
    %         den_vel_Li2(i,:)=den_vel_Li2(i,:)+spwt*2*q_e*vel_Li2(p,:)*(1-hz);   %���䵽��i������ϵ��ٶ�
    %         den_vel_Li2(i+1,:)=den_vel_Li2(i,:)+spwt*2*q_e*vel_Li2(p,:)*hz;     %���䵽��i+1������ϵ��ٶ�
    %     end
    %
    %     %����Li3+�ĵ�ɺ͵����ܶȵ��ٽ������
    %     for p=1:N_Li3
    %
    %         fi=1+pos_Li3(p)/dz;                 %ʵ�ʸ��λ�ã�Ϊ������
    %         i=floor(fi);                        %��Ӧ�������λ��
    %         hz=fi-i;                            %�������i�����֮��ľ���ռ���벽���ı���
    %
    %         %��Li3+��ɰ��ձ������䵽�ٽ������������
    %         chg(i,4)=chg(i,4)+(1-hz);      %���䵽��i������ϵĵ�ɱ���
    %         chg(i+1,4)=chg(i+1,4)+hz;      %���䵽��i+1������ϵĵ�ɱ���
    %
    %         count_Li(i,3)=count_Li(i,3)+1;       %��ÿ����Ԫ���ڵ�Li3+����
    %
    %         %���ٶȰ��������䵽�ٽ������������
    %         den_vel_Li3(i,:)=den_vel_Li3(i,:)+spwt*3*q_e*vel_Li3(p,:)*(1-hz);   %���䵽��i������ϵ��ٶ�
    %         den_vel_Li3(i+1,:)=den_vel_Li3(i,:)+spwt*3*q_e*vel_Li3(p,:)*hz;     %���䵽��i+1������ϵ��ٶ�
    %     end
    %
    %     %����ÿ������ϵ������ܶ�
    %     for count_i=2:nz-1
    %         count_g_Li(count_i,:)=0.5*(count_Li(count_i-1,:)+count_Li(count_i,:))/(2*dz);
    %     end
    %     count_g_Li(1,:)=count_Li(1,:)/dz;      %��һ�����������ܶ�
    %     count_g_Li(nz,:)=count_Li(nz-1,:)/dz;  %���һ�����������ܶ�
    %
    %     J_i=1e5*(den_vel_Li1.*count_g_Li(:,1)+den_vel_Li2.*count_g_Li(:,2)+den_vel_Li3.*count_g_Li(:,3));  %������ÿ������ϵĵ����ܶ�
    
    %��ʼ���ӵ��˶�ѭ��
    for ts_e=1:dt_i/dt_e                                                                 %�����˶�һ���������˶�dt_i/dt_e��
        
        chg=zeros(nz,4);                %���ڴ��ÿ�������䵽�ĵ������1~4�зֱ����:���ӣ�Li1+��Li2+��Li3+
        den_vel_Li1=zeros(nz,3);        %���ڴ��ÿ�������Li+�ĵ����ܶȷ���
        den_vel_Li2=zeros(nz,3);        %���ڴ��ÿ�������Li2+�ĵ����ܶȷ���
        den_vel_Li3=zeros(nz,3);        %���ڴ��ÿ�������Li3+�ĵ����ܶȷ���
        count_Li=zeros(nz-1,3);            %���ڼ���ÿ����Ԫ���ڵ���������1~3�зֱ����Li+��Li2+��Li3+
        count_g_Li=zeros(nz,3);            %���ڼ���ÿ������ϵ��������ܶȣ�1~3�зֱ����Li+��Li2+��Li3+
        distri_Li1=zeros(N_Li1,nz-1);     %����ʢ��ÿ����Ԫ����Li+������
        distri_Li2=zeros(N_Li2,nz-1);     %����ʢ��ÿ����Ԫ����Li2+������
        distri_Li3=zeros(N_Li3,nz-1);     %����ʢ��ÿ����Ԫ����Li3+������
        
        for p=1:N_Li1
            fi=1+pos_Li1(p)/dz;                 %ʵ�ʸ��λ�ã�Ϊ������
            i=floor(fi);                        %��Ӧ�������λ��
            count_Li(i,1)=count_Li(i,1)+1;       %��ÿ����Ԫ���ڵ�Li+����
            distri_Li1(count_Li(i,1),i)=p;  %ÿ����Ԫ���ڵ�Li+������
        end
        for p=1:N_Li2
            fi=1+pos_Li2(p)/dz;                 %ʵ�ʸ��λ�ã�Ϊ������
            i=floor(fi);                        %��Ӧ�������λ��
            count_Li(i,2)=count_Li(i,2)+1;       %��ÿ����Ԫ���ڵ�Li2+����
            distri_Li2(count_Li(i,2),i)=p;  %ÿ����Ԫ���ڵ�Li2+������
        end
        for p=1:N_Li3
            fi=1+pos_Li3(p)/dz;                 %ʵ�ʸ��λ�ã�Ϊ������
            i=floor(fi);                        %��Ӧ�������λ��
            count_Li(i,3)=count_Li(i,3)+1;       %��ÿ����Ԫ���ڵ�Li3+����
            distri_Li3(count_Li(i,3),i)=p;  %ÿ����Ԫ���ڵ�Li3+������
        end
        
        %���Ƚ���MCC����
        den_Li=count_Li/(100*dz);       %ÿ����Ԫ����Li+��Li2+��Li3+�����ܶȣ���λ��cm-1
        nu_t=(n0+den_Li(:,1)+den_Li(:,2))*R_i_max; %ÿ����Ԫ�����ܵĵ�����ײƵ�ʣ�����Ҫ��Ҫ������ӳ���spwt???
        P_emax=1-exp(-nu_t*dt_e);                    %ÿ����Ԫ��������һ�����ӵ������ײ����
        count_e=zeros(nz-1,1);                     %������ÿ����Ԫ���ڵĵ��ӽ��м���
        distri_e=zeros(N_e,nz-1);                   %����ʢ��ÿ����Ԫ���ڵĵ�������
        for p=1:N_e
            fi=1+pos_e(p)/dz;            %ʵ�ʸ��λ�ã�Ϊ������
            i=floor(fi);               %��Ӧ�������λ��
            count_e(i)=count_e(i)+1;   %��Ӧ��Ԫ���ڵļ���������һ��
            distri_e(count_e(i),i)=p;  %ÿ����Ԫ���ڵĵ�������
        end
        N_e_mcc=P_emax.*count_e;        %ÿ����Ԫ���ڳ�ȡ�ĵ��Ӹ���
        N_e_sam=zeros(ceil(max(N_e_mcc)),nz-1); %�������������ÿ����Ԫ�������ĵ��ӵ����
        Ek_e=zeros(ceil(max(N_e_mcc)),nz-1);    %�þ���������ų��������ĵ��ӵĶ���
        nu_Li1=zeros(ceil(max(N_e_mcc)),nz-1);  %�þ���������ų���������e+Li0->Li+����ײƵ��
        nu_Li2=zeros(ceil(max(N_e_mcc)),nz-1);  %�þ���������ų���������e+Li+->Li2+����ײƵ��
        nu_Li3=zeros(ceil(max(N_e_mcc)),nz-1);  %�þ���������ų���������e+Li2+->Li3+����ײƵ��
        for i=1:nz-1            %��ʼ��ÿ����Ԫ��������ؿ������
            N_e_sam(1:ceil(N_e_mcc(i)),i)=(randperm(count_e(i),ceil(N_e_mcc(i))))';  %randperm(n,k)�������ڴ�n�����������ȡ����ͬ��k����
            for j=1:N_e_mcc(i)
                Ek_e(j,i)=0.5*m_e*sum((vel_e(distri_e(N_e_sam(j,i),i),:)).^2)/q_e;  %���ÿ����ȡ�����ĵ��ӵĶ��ܣ���ת��ΪeV��λ
                nu_Li1(j,i)=n0*10^(B_ir(1,:)*[1;log10(Ek_e(j,i));(log10(Ek_e(j,i)))^2;(log10(Ek_e(j,i)))^3;(log10(Ek_e(j,i)))^4;(log10(Ek_e(j,i)))^5])/nu_t(i); %�������ݿ��е����������������������Ӷ�Ӧ�Ĳ���Li+����ײƵ��
                nu_Li2(j,i)=den_Li(i,1)*10^(B_ir(2,:)*[1;log10(Ek_e(j,i));(log10(Ek_e(j,i)))^2;(log10(Ek_e(j,i)))^3;(log10(Ek_e(j,i)))^4;(log10(Ek_e(j,i)))^5])/nu_t(i); %����Li2+����ײƵ��
                nu_Li3(j,i)=den_Li(i,2)*10^(B_ir(3,:)*[1;log10(Ek_e(j,i));(log10(Ek_e(j,i)))^2;(log10(Ek_e(j,i)))^3;(log10(Ek_e(j,i)))^4;(log10(Ek_e(j,i)))^5])/nu_t(i); %����Li3+����ײƵ��
            end
        end
        nu_Li1_Li2=nu_Li1+nu_Li2;             %�������еĵڶ����һ����nu_Li1
        nu_Li1_Li2_Li3=nu_Li1+nu_Li2+nu_Li3;  %�������еĵ�����
        for i=1:nz-1               %����ÿ����Ԫ��
            for j=1:N_e_mcc(i)     %����ÿ����Ԫ���ڳ���ĵ���
                R=rand();          %RΪ��0��1��֮��������
                if R<=nu_Li1(j,i)&&Ek_e(j,i)>=E_i_Li1                                        %Li0->Li+�¼�����
                    N_Li1=N_Li1+1;                                                           %������һ���µ�Li+
                    N_e=N_e+1;                                                               %������һ���µĵ���
                    pos_Li1(N_Li1)=pos_e(distri_e(N_e_sam(j,i),i));                              %�²�����Li+��λ����������ӵ���ͬ
                    vel_Li1(N_Li1,:)=randraw('maxwell',k*T/m_Li,1)*random_unit_vector(3,1)'; %�²�����Li+���ٶȷ������˹Τ�ֲ�
                    pos_e(N_e)=pos_Li1(N_Li1);                                               %�²����ĵ��Ӻ��²�����Li+����ͬ��λ��
                    vel_e_temp=0.5*sqrt(2*q_e*(Ek_e(j,i)-E_i_Li1)/m_e)*random_unit_vector(3,2);  %��Ϊ���Ӷ����ӵ�����û�й��ף����Ӷ��ܼ�ȥ�����ܺ�ƽ�ָ�ɢ����Ӻ��²����ĵ��ӣ�random_unit_vector����������������ͬ�Եĵ���
                    vel_e(N_e,:)=vel_e_temp(:,1)';                                           %�����ɵĵ��ӵ��ٶ�
                    vel_e(distri_e(N_e_sam(j,i),i),:)=vel_e_temp(:,2)';                        %ɢ����ӵ��ٶ�
                end
                if R>nu_Li1(j,i)&&R<=nu_Li1_Li2(j,i)&&Ek_e(j,i)>=E_i_Li2                         %Li+->Li2+�¼�����
                    N_Li2=N_Li2+1;                                                          %������һ���µ�Li2+
                    N_e=N_e+1;                                                              %������һ���µĵ���
                    O_Li1=distri_Li1(randperm(count_Li(i,1),1),i);                          %������Li+��Ӧ����ţ�����ڵ�i����Ԫ����ѡ��
                    N_Li1=N_Li1-1;                                                          %Li+��Ŀ����һ��
                    pos_Li2(N_Li2)=pos_Li1(O_Li1);                                          %�²�����Li2+��λ����������Li+��λ����ͬ
                    vel_Li2(N_Li2,:)=vel_Li1(O_Li1,:);                                      %�²�����Li2+���ٶ�����������Li+���ٶ���ͬ
                    pos_e(N_e)=pos_Li2(N_Li2);                                              %�²����ĵ��ӵ�λ�ú�Li2+������λ����ͬ
                    vel_e_temp=0.5*sqrt(2*q_e*(Ek_e(j,i)-E_i_Li2)/m_e)*random_unit_vector(3,2); %��Ϊ���Ӷ����ӵ�����û�й��ף����Ӷ��ܼ�ȥ�����ܺ�ƽ�ָ�ɢ����Ӻ��²����ĵ��ӣ�random_unit_vector����������������ͬ�Եĵ���
                    vel_e(N_e,:)=vel_e_temp(:,1)';                                          %�����ɵĵ��ӵ��ٶ�
                    vel_e(distri_e(N_e_sam(j,i),i),:)=vel_e_temp(:,2)';                       %ɢ����ӵ��ٶ�
                end
                if R>nu_Li1_Li2(j,i)&&R<=nu_Li1_Li2_Li3(j,i)&&Ek_e(j,i)>=E_i_Li3                 %Li2+->Li3+�¼�����
                    N_Li3=N_Li3+1;                                                          %������һ���µ�Li3+
                    N_e=N_e+1;                                                              %������һ���µĵ���
                    O_Li2=distri_Li2(randperm(count_Li(i,2),1),i);                          %������Li2+��Ӧ����ţ�����ڵ�i����Ԫ����ѡ��
                    N_Li2=N_Li2-1;                                                          %Li2+��Ŀ����һ��
                    pos_Li3(N_Li3)=pos_Li2(O_Li2);                                          %�²�����Li3+��λ����������Li2+��λ����ͬ
                    vel_Li3(N_Li3,:)=vel_Li2(O_Li2,:);                                      %�²�����Li3+���ٶ�����������Li2+���ٶ���ͬ
                    pos_e(N_e)=pos_Li3(N_Li3);                                              %�²����ĵ��ӵ�λ�ú�Li3+������λ����ͬ
                    vel_e_temp=0.5*sqrt(2*q_e*(Ek_e(j,i)-E_i_Li3)/m_e)*random_unit_vector(3,2); %��Ϊ���Ӷ����ӵ�����û�й��ף����Ӷ��ܼ�ȥ�����ܺ�ƽ�ָ�ɢ����Ӻ��²����ĵ��ӣ�random_unit_vector����������������ͬ�Եĵ���
                    vel_e(N_e,:)=vel_e_temp(:,1)';                                          %�����ɵĵ��ӵ��ٶ�
                    vel_e(distri_e(N_e_sam(j,i),i),:)=vel_e_temp(:,2)';                       %ɢ����ӵ��ٶ�
                end
            end
        end
        
        %����Li+�ĵ�ɺ͵����ܶȵ��ٽ������
        for p=1:N_Li1
            
            fi=1+pos_Li1(p)/dz;                 %ʵ�ʸ��λ�ã�Ϊ������
            i=floor(fi);                        %��Ӧ�������λ��
            hz=fi-i;                            %�������i�����֮��ľ���ռ���벽���ı���
            
            %��Li+��ɰ��ձ������䵽�ٽ������������
            chg(i,2)=chg(i,2)+(1-hz);      %���䵽��i������ϵĵ�ɱ���
            chg(i+1,2)=chg(i+1,2)+hz;      %���䵽��i+1������ϵĵ�ɱ���
            
            count_Li(i,1)=count_Li(i,1)+1;       %��ÿ����Ԫ���ڵ�Li+����
            
            %���ٶȰ��������䵽�ٽ������������
            den_vel_Li1(i,:)=den_vel_Li1(i,:)+spwt*q_e*vel_Li1(p,:)*(1-hz);   %���䵽��i������ϵ��ٶ�
            den_vel_Li1(i+1,:)=den_vel_Li1(i,:)+spwt*q_e*vel_Li1(p,:)*hz;     %���䵽��i+1������ϵ��ٶ�
        end
        
        %����Li2+�ĵ�ɺ͵����ܶȵ��ٽ������
        for p=1:N_Li2
            
            fi=1+pos_Li2(p)/dz;                 %ʵ�ʸ��λ�ã�Ϊ������
            i=floor(fi);                        %��Ӧ�������λ��
            hz=fi-i;                            %�������i�����֮��ľ���ռ���벽���ı���
            
            %��Li2+��ɰ��ձ������䵽�ٽ������������
            chg(i,3)=chg(i,3)+(1-hz);      %���䵽��i������ϵĵ�ɱ���
            chg(i+1,3)=chg(i+1,3)+hz;      %���䵽��i+1������ϵĵ�ɱ���
            
            count_Li(i,2)=count_Li(i,2)+1;       %��ÿ����Ԫ���ڵ�Li2+����
            
            %���ٶȰ��������䵽�ٽ������������
            den_vel_Li2(i,:)=den_vel_Li2(i,:)+spwt*2*q_e*vel_Li2(p,:)*(1-hz);   %���䵽��i������ϵ��ٶ�
            den_vel_Li2(i+1,:)=den_vel_Li2(i,:)+spwt*2*q_e*vel_Li2(p,:)*hz;     %���䵽��i+1������ϵ��ٶ�
        end
        
        %����Li3+�ĵ�ɺ͵����ܶȵ��ٽ������
        for p=1:N_Li3
            
            fi=1+pos_Li3(p)/dz;                 %ʵ�ʸ��λ�ã�Ϊ������
            i=floor(fi);                        %��Ӧ�������λ��
            hz=fi-i;                            %�������i�����֮��ľ���ռ���벽���ı���
            
            %��Li3+��ɰ��ձ������䵽�ٽ������������
            chg(i,4)=chg(i,4)+(1-hz);      %���䵽��i������ϵĵ�ɱ���
            chg(i+1,4)=chg(i+1,4)+hz;      %���䵽��i+1������ϵĵ�ɱ���
            
            count_Li(i,3)=count_Li(i,3)+1;       %��ÿ����Ԫ���ڵ�Li3+����
            
            %���ٶȰ��������䵽�ٽ������������
            den_vel_Li3(i,:)=den_vel_Li3(i,:)+spwt*3*q_e*vel_Li3(p,:)*(1-hz);   %���䵽��i������ϵ��ٶ�
            den_vel_Li3(i+1,:)=den_vel_Li3(i,:)+spwt*3*q_e*vel_Li3(p,:)*hz;     %���䵽��i+1������ϵ��ٶ�
        end
        
        %����ÿ������ϵ������ܶ�
        for count_i=2:nz-1
            count_g_Li(count_i,:)=0.5*(count_Li(count_i-1,:)+count_Li(count_i,:))/(2*dz);
        end
        count_g_Li(1,:)=count_Li(1,:)/dz;      %��һ�����������ܶ�
        count_g_Li(nz,:)=count_Li(nz-1,:)/dz;  %���һ�����������ܶ�
        
        J_i=1e5*(den_vel_Li1.*count_g_Li(:,1)+den_vel_Li2.*count_g_Li(:,2)+den_vel_Li3.*count_g_Li(:,3));  %������ÿ������ϵĵ����ܶ�
        
        
        
        
        %������ӵĵ���ܶȺ͵����ܶ�
        chg(:,1)=zeros(nz,1);                                                       %���ڴ��ÿ�������䵽�ĵ��ӵĵ����
        %den=zeros(nz,1);                                                            %ÿ�����ĵ���ܶȣ���λ��C/m
        den_vel_e=zeros(nz,3);                                                        %���ڴ��ÿ������ϵĵ����ܶȷ���
        B_mic_old=B_mic;                                                            %B_mic_old���ں�B_mic��ƽ�����Ӷ���������ٶ�
        pos_e_half=zeros(N_e,1);                                                      %���ڴ�Ű�ʱ��ʱ������λ��
        
        for p_half=1:N_e                                                            %���������ӵ�λ����ǰ�ƶ�����ʱ�̴�
            pos_e_half(p_half,1)=pos_e(p_half,1)-vel_e(p_half,3)*0.5*dt_e;
            if pos_e_half(p_half,1)<0                                                 %����if������ʾ�����Ա߽�����
                pos_e_half(p_half,1)=pos_e_half(p_half,1)+z_max;
            end
            if pos_e_half(p_half,1)>z_max
                pos_e_half(p_half,1)=pos_e_half(p_half,1)-z_max;
            end
        end
        % pos_half=pos-vel(:,3)*0.5*dt;                                               %�����ӵ�λ����ǰ�ƶ����ʱ�䲽��������ð�ʱ�̴��ĵ����ܶ�
        
        %���Ƚ����е��ӵĵ�ɰ���������ԭ�����
        
        count=zeros(nz-1,1);            %���ڼ���ÿ����Ԫ���ڵ�������
        count_g=zeros(nz,1);            %���ڼ���ÿ������ϵ��������ܶ�
        
        %������ӵĵ��������ܶȵ��ٽ������
        for p=1:N_e
            fi=1+pos_e(p)/dz;            %ʵ�ʸ��λ�ã�Ϊ������
            i=floor(fi);               %��Ӧ�������λ��
            hz=fi-i;                   %�������i�����֮��ľ���ռ���벽���ı���
            
            %�����ӵ�ɰ��ձ������䵽�ٽ������������
            chg(i,1)=chg(i,1)+(1-hz);      %���䵽��i������ϵĵ�ɱ���
            chg(i+1,1)=chg(i+1,1)+hz;      %���䵽��i+1������ϵĵ�ɱ���
            
            fi_half=1+pos_e_half(p)/dz;  %��ʱ�̴��ĸ��λ��
            i_half=floor(fi_half);     %��ʱ�̴����λ�ö�Ӧ���������λ��
            hz_half=fi_half-i_half;    %��ʱ�̴��������i_half�����֮��ľ���ռ���벽���ı���
            
            count(i_half)=count(i_half)+1;  %ͳ��ÿ����Ԫ���ڵ�������Ŀ
            
            
            
            %���ٶȰ��������䵽�ٽ������������
            den_vel_e(i_half,:)=den_vel_e(i_half,:)+spwt*(-q_e)*vel_e(p,:)*(1-hz_half);  %���䵽��i������ϵ��ٶ�
            den_vel_e(i_half+1,:)=den_vel_e(i_half+1,:)+spwt*(-q_e)*vel_e(p,:)*hz_half;  %���䵽��i+1������ϵ��ٶ�
        end
        
        den=spwt*q_e*(-1*chg(:,1)+1*chg(:,2)+2*chg(:,3)+3*chg(:,4))/dz;       %ÿ������ϵĵ���ܶ�
        %J=spwt*q_e*(-1)*den_vel;  %ÿ������ϵĵ����ܶ�
        %J=den_vel.*den;
        %J=1e9*den_vel;
        
        %����ÿ������ϵĵ����ܶ�
        for count_i=2:nz-1
            count_g(count_i)=0.5*(count(count_i-1)+count(count_i))/(2*dz);
        end
        count_g(1)=count(1)/dz;      %��һ�����ĵ����ܶ�
        count_g(nz)=count(nz-1)/dz;  %���һ�����ĵ����ܶ�
        
        J=1e5*den_vel_e.*count_g+J_i;          %����ÿ������ϵĵ����ܶ�
        
        %Ӧ�ñ߽��������ڱ߽��Ͻ���һ��ĳ����й��ף���˱߽��ϵ��ܶ�Ӧ����2
        den(1)=2*den(1);   %��߽�
        den(nz)=2*den(nz); %�ұ߽�
        
        %��ʼ��Ⲵ�ɷ��̣��õ�����ϵĵ��ƣ����÷���Ϊ��˹-���¶���������G-S��
        % phi_tem=phi;
        
        
        %         phi=G_S_SOR(phi,den,tol,sor,N_phi,dz,nz);
        
        
        
        %         for j=1:N_phi     %��������
        %
        %
        %             for i=2:nz-1  %�м���ĵ���
        %                 phi_tem(i)=0.5*(den(i)/EPS0*dz*dz+phi(i-1)+phi(i+1));   %�������޲�ַ����õ�
        %             end
        %
        %             phi_tem(1)=phi(2);     %��ڴ�Ϊ��һ��߽�����
        %             phi_tem(nz)=phi(nz-1); %���ڴ�ҲΪ��һ��߽�����
        %
        %             if mod(j,10)==0
        %                 R=norm(phi_tem-phi); %������ε����Ĳв�
        %                 if (R<=tol)          %��������
        %                     phi=phi_tem;
        %                     fprintf('Congratulations! The solver of the electric potential is converged at %d step\n',j);
        %                     break;
        %                 end
        %             end
        %             phi=phi_tem;
        %         end
        %         if j>=N_phi                  %���������ⲻ����������ֹ���������
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
        
        
        
        
        
        %���ÿ������ϵľ��糡��������һά�����ֻ�ܵõ�z����ĵ糡ǿ��
        E_sta(2:nz-1)=(phi(1:nz-2)-phi(3:nz))/2/dz;      %�м���ľ��糡
        E_sta(1)=(phi(1)-phi(2))/dz;                     %��ڴ����糡
        E_sta(nz)=(phi(nz-1)-phi(nz))/dz;                %���ڴ����糡
        
        %����΢����ų��ĸ�������
        E_mic_xold=E0*cos(w*1e9*2*pi()*(ts_e-2)*dt_e);
        E_mic_yold=E0*sin(w*1e9*2*pi()*(ts_e-2)*dt_e);
        E_mic(1,1)=E0*cos(w*1e9*2*pi()*(ts_e-1)*dt_e);              %��ڴ�x����ĵ糡ǿ��E_x=E0cos(wt)
        % B_mic_xin_full=B0*sin(w*1e9*2*pi()*(ts-0.5)*dt);       %��ڴ�x����Ÿ�Ӧǿ�ȣ�ʱ��Ϊ��
        B_mic_xin=B0*sin(w*1e9*2*pi()*(ts_e-0.5)*dt_e);             %��ڴ�x����ĴŸ�Ӧǿ��ΪB_x=B0sin(wt)����ʱ��������ƶ���t/2
        E_mic(1,2)=E0*sin(w*1e9*2*pi()*(ts_e-1)*dt_e);              %��ڴ�y����ĵ糡ǿ��E_y=E0sin(wt)
        %B_mic_yin_full=B0*cos(w*1e9*2*pi()*(ts-0.5)*dt);       %��ڴ�y����Ÿ�Ӧǿ�ȣ�ʱ��Ϊ��
        B_mic_yin=B0*cos(w*1e9*2*pi()*(ts_e-0.5)*dt_e);             %��ڴ�y����ĴŸ�Ӧǿ��B_y=B0cos(wt)����ʱ��������ƶ���t/2
        %B_mic(1,1)=B_mic_xin+dt/(2*dz)*(E_mic(2,2)-E_mic(1,2));                   %��Bx�ڿռ��ʱ���϶�����ƶ���z/2�ͦ�t/2
        B_mic(1,1)=B_mic_xin+dz/2*(1/(c^2)*(E_mic(1,2)-E_mic_yold)/dt_e+MU0*J(1,2)); %��ڴ���Bx����Bx�ڿռ�������ƶ���z/2�������
        B_mic(1,2)=B_mic_yin-dz/2*(1/(c^2)*(E_mic(1,1)-E_mic_xold)/dt_e+MU0*J(1,1)); %��ڴ���By����By�ڿռ�������ƶ���z/2�������
        %B_mic(1,2)=B_mic_yin-dt/(2*dz)*(E_mic(2,1)-E_mic(1,1));                   %��By�ڿռ��ʱ���϶�����ƶ���z/2�ͦ�t/2
        
        % E_mic(nz,1)=0; %���ڴ�x����ĵ糡ǿ��0
        % E_mic(nz,2)=0; %���ڴ�y����ĵ糡ǿ��0
        
        B_mic_xout=B_mic(nz,1);
        B_mic_yout=B_mic(nz,2);
        
        %B_mic_xout=0;  %���ڴ�x����ĴŸ�Ӧǿ��0
        % B_mic_yout=0; %���ڴ�y����ĴŸ�Ӧǿ��0
        % B_mic(nz-1,1)=B_mic_xout-dz/2*MU0*J(nz,2); %���ڴ���Bx����ǰ���
        % B_mic(nz-1,2)=B_mic_yout+dz/2*MU0*J(nz,1); %���ڴ���By����ǰ���
        
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
        
        %������ϵľ��糡�͵�ų�ȫ�����䵽ÿ�������ϣ����⻹Ҫ�����ⲿ����������Ĵų�
        E_e=zeros(max_part,3);           %ÿ����������λ�ô��ĵ糡ǿ��
        B_e=zeros(max_part,3);           %ÿ����������λ�ô��ĴŸ�Ӧǿ��
        %B_e(:,3)=B_e(:,3)+B_ex_z;        %z����ĴŸ�Ӧǿ����Ҫ�����ⲿ�Ÿ�Ӧǿ��
        B_mic_ave=0.5*(B_mic_old+B_mic); %����ƽ������������ʱ���ֵƽ��������ʱ����
        p=1;                             %���ڱ�����������
        while p<=N_e
            fi=1+pos_e(p)/dz;  %ʵ�ʸ��λ�ã�Ϊ������
            i=floor(fi);     %��Ӧ�������λ��
            hz=fi-i;         %�������i�����֮��ľ���ռ���벽���ı���
            E_e(p,3)=E_e(p,3)+E_sta(i)*(1-hz)+E_sta(i+1)*hz;     %��������֮��ľ��糡
            E_e(p,:)=E_e(p,:)+E_mic(i,:)*(1-hz)+E_mic(i+1,:)*hz; %����΢�������ĵ糡
            %B_e(p,3)=B_e(p,3)+B_ex_z(i)*(1-hz)+B_ex_z(i+1)*hz;   %������ϵľ��ų������ÿ������
            B_e(p,3)=1e-4*(2595.7+32.4092*(1000*pos_e(p))-3.1038*(1000*pos_e(p))^2-...   %��ֵ�õ�����λ�ô��ĴŸ�Ӧǿ��
                0.0130847*(1000*pos_e(p))^3+0.00198809*(1000*pos_e(p))^4-2.48996e-5*(1000*pos_e(p))^5+8.9815e-8*(1000*pos_e(p))^6);
            if fi<=1.5                                           %�ĸ�if����ȷ������λ�ã��Ӷ�Ϊ֮����΢�������ĴŸ�Ӧǿ��
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
            
            %����BORIS���������һʱ�̵��ٶȺ�λ��
            vel_e(p,:)=UpdateVelocity(E_e(p,:),B_e(p,:),vel_e(p,:),dt_e); %�����ٶ�
            pos_e(p,1)=pos_e(p,1)+vel_e(p,3)*dt_e;                          %����λ��
            
            %�������ձ߽���������������໥���ú�������
            if pos_e(p,1)<=0             %��߽�
                pos_e(p,1)=pos_e(N_e,1); %�����һ�����ӵ�λ���滻������λ�ô�
                vel_e(p,:)=vel_e(N_e,:); %�����һ�����ӵ��ٶ��滻������λ�ô�
                N_e=N_e-1;               %����������һ��
            end
            if pos_e(p,1)>z_max          %�ұ߽�
                pos_e(p,1)=pos_e(N_e,1); %�����һ�����ӵ�λ���滻������λ�ô�
                vel_e(p,:)=vel_e(N_e,:); %�����һ�����ӵ��ٶ��滻������λ�ô�
                N_e=N_e-1;               %����������һ��
            end
            
            p=p+1;
        end
        
        fprintf('ts_e: %d, ts_i: %d,N_Li1: %d, N_Li2: %d, N_Li3: %d, P_emax: %f\n Average Energy: %f eV\n', ts_e,ts_i,N_Li1,N_Li2,N_Li3,max(P_emax),mean(Ek_e(1,:)));
        % plot(1:nz,E_mic(:,1));
        % pause(0.1);
    end
    
    %������ϵľ��糡�͵�ų�ȫ�����䵽ÿ�������ϣ����⻹Ҫ�����ⲿ����������Ĵų�
    E_Li1=zeros(max_part,3);           %ÿ��Li+����λ�ô��ĵ糡ǿ��
    E_Li2=zeros(max_part,3);          %ÿ��Li2+����λ�ô��ĵ糡ǿ��
    E_Li3=zeros(max_part,3);          %ÿ��Li3+����λ�ô��ĵ糡ǿ��
    B_Li1=zeros(max_part,3);           %ÿ��Li+����λ�ô��ĴŸ�Ӧǿ��
    B_Li2=zeros(max_part,3);          %ÿ��Li2+����λ�ô��ĴŸ�Ӧǿ��
    B_Li3=zeros(max_part,3);          %ÿ��Li3+����λ�ô��ĴŸ�Ӧǿ��
    %B_e(:,3)=B_e(:,3)+B_ex_z;        %z����ĴŸ�Ӧǿ����Ҫ�����ⲿ�Ÿ�Ӧǿ��
    B_mic_ave=0.5*(B_mic_old+B_mic); %����ƽ������������ʱ���ֵƽ��������ʱ����
    %     N_i=[N_Li1,N_Li2,N_Li3];
    %     for i_Li=1:3
    
    p=1;                             %���ڱ�������Li+
    while p<=N_Li1
        fi=1+pos_e(p)/dz;  %ʵ�ʸ��λ�ã�Ϊ������
        i=floor(fi);     %��Ӧ�������λ��
        hz=fi-i;         %�������i�����֮��ľ���ռ���벽���ı���
        %             i_str=num2str(i);
        %             p_str=num2str(p);
        %             hz_str=num2str(hz);
        %             eval(['E_Li',i_str,'(',p_str,',',num2str(3),')','=','E_Li',i_str,'(',p_str,',',num2str(3),')','+E_sta(',i_str,')*(',num2str(1-hz),'+E_sta(',num2str(i+1),')*',num2str(hz)]);     %��������֮��ľ��糡
        %eval(['E_Li',num2str(i),'(',num2str(p),',',':',')'])=eval(['E_Li',num2str(i),'(',num2str(p),',',':',')'])+E_mic(i,:)*(1-hz)+E_mic(i+1,:)*hz; %����΢�������ĵ糡
        %             eval(['E_Li',i_str,'(',p_str,',:',')','=','E_Li',i_str,'(',p_str,',:',')','+E_mic(',i_str,',:',')*(',num2str(1-hz),'+E_mic(',num2str(i+1),',:)','*',num2str(hz)]);
        E_Li1(p,3)=E_Li1(p,3)+E_sta(i)*(1-hz)+E_sta(i+1)*hz;     %��������֮��ľ��糡
        E_Li1(p,:)=E_Li1(p,:)+E_mic(i,:)*(1-hz)+E_mic(i+1,:)*hz; %����΢�������ĵ糡
        %B_e(p,3)=B_e(p,3)+B_ex_z(i)*(1-hz)+B_ex_z(i+1)*hz;   %������ϵľ��ų������ÿ������
        %             eval(['B_Li',i_str,'(',p_str,',3)=1e-4*(2595.7+32.4092*(1000*pos_Li',i_str,'(',p_str,')',')-3.1038*(1000*pos_Li',i_str,'(',p_str,'))^2-0.0130847*(1000*pos_Li'...
        %                 ,i_str,'(',p_str,'))^3+0.00198809*(1000*pos_Li',i_str,'(',p_str,'))^4-2.48996e-5*(1000*pos_Li',i_str,'(',p_str,'))^5+8.9815e-8*(1000*pos_Li',...
        %                 i_str,'(',p_str,'))^6)']);
        B_Li1(p,3)=1e-4*(2595.7+32.4092*(1000*pos_Li1(p))-3.1038*(1000*pos_Li1(p))^2-...   %��ֵ�õ�����λ�ô��ĴŸ�Ӧǿ��
            0.0130847*(1000*pos_Li1(p))^3+0.00198809*(1000*pos_Li1(p))^4-2.48996e-5*(1000*pos_Li1(p))^5+8.9815e-8*(1000*pos_Li1(p))^6);
        if fi<=1.5                                           %�ĸ�if����ȷ������λ�ã��Ӷ�Ϊ֮����΢�������ĴŸ�Ӧǿ��
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
        
        %����BORIS���������һʱ�̵��ٶȺ�λ��
        vel_Li1(p,:)=UpdateVelocity(E_Li1(p,:),B_Li1(p,:),vel_Li1(p,:),dt_i); %�����ٶ�
        pos_Li1(p,1)=pos_Li1(p,1)+vel_Li1(p,3)*dt_i;                          %����λ��
        
        %�������ձ߽�������Li+����໥���ú�������
        if pos_Li1(p,1)<0             %��߽�
            pos_Li1(p,1)=pos_Li1(N_Li1,1); %�����һ��Li+��λ���滻������λ�ô�
            vel_Li1(p,:)=vel_Li1(N_Li1,:); %�����һ��Li+���ٶ��滻������λ�ô�
            N_Li1=N_Li1-1;               %Li+������һ��
        end
        if pos_Li1(p,1)>z_max          %�ұ߽�
            pos_Li1(p,1)=pos_Li1(N_Li1,1); %�����һ��Li+��λ���滻������λ�ô�
            vel_Li1(p,:)=vel_Li1(N_Li1,:); %�����һ��Li+���ٶ��滻������λ�ô�
            N_Li1=N_Li1-1;               %Li+������һ��
        end
        
        p=p+1;
    end
    %     end
    p=1;                             %���ڱ�������Li2+
    while p<=N_Li2
        fi=1+pos_e(p)/dz;  %ʵ�ʸ��λ�ã�Ϊ������
        i=floor(fi);     %��Ӧ�������λ��
        hz=fi-i;         %�������i�����֮��ľ���ռ���벽���ı���
        
        E_Li2(p,3)=E_Li2(p,3)+E_sta(i)*(1-hz)+E_sta(i+1)*hz;     %��������֮��ľ��糡
        E_Li2(p,:)=E_Li2(p,:)+E_mic(i,:)*(1-hz)+E_mic(i+1,:)*hz; %����΢�������ĵ糡
        
        B_Li2(p,3)=1e-4*(2595.7+32.4092*(1000*pos_Li1(p))-3.1038*(1000*pos_Li1(p))^2-...   %��ֵ�õ�����λ�ô��ĴŸ�Ӧǿ��
            0.0130847*(1000*pos_Li1(p))^3+0.00198809*(1000*pos_Li1(p))^4-2.48996e-5*(1000*pos_Li1(p))^5+8.9815e-8*(1000*pos_Li1(p))^6);
        if fi<=1.5                                           %�ĸ�if����ȷ������λ�ã��Ӷ�Ϊ֮����΢�������ĴŸ�Ӧǿ��
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
        
        %����BORIS���������һʱ�̵��ٶȺ�λ��
        vel_Li2(p,:)=UpdateVelocity(E_Li2(p,:),B_Li2(p,:),vel_Li2(p,:),dt_i); %�����ٶ�
        pos_Li2(p,1)=pos_Li2(p,1)+vel_Li2(p,3)*dt_i;                          %����λ��
        
        %�������ձ߽�������Li2+����໥���ú�������
        if pos_Li2(p,1)<0             %��߽�
            pos_Li2(p,1)=pos_Li2(N_Li2,1); %�����һ��Li2+��λ���滻������λ�ô�
            vel_Li2(p,:)=vel_Li2(N_Li2,:); %�����һ��Li2+���ٶ��滻������λ�ô�
            N_Li2=N_Li2-1;               %Li2+������һ��
        end
        if pos_Li2(p,1)>z_max          %�ұ߽�
            pos_Li2(p,1)=pos_Li2(N_Li2,1); %�����һ��Li2+��λ���滻������λ�ô�
            vel_Li2(p,:)=vel_Li2(N_Li2,:); %�����һ��Li2+���ٶ��滻������λ�ô�
            N_Li2=N_Li2-1;               %Li2+������һ��
        end
        
        p=p+1;
    end
    
    
    p=1;                             %���ڱ�������Li3+
    while p<=N_Li3
        fi=1+pos_e(p)/dz;  %ʵ�ʸ��λ�ã�Ϊ������
        i=floor(fi);     %��Ӧ�������λ��
        hz=fi-i;         %�������i�����֮��ľ���ռ���벽���ı���
        
        E_Li3(p,3)=E_Li3(p,3)+E_sta(i)*(1-hz)+E_sta(i+1)*hz;     %��������֮��ľ��糡
        E_Li3(p,:)=E_Li3(p,:)+E_mic(i,:)*(1-hz)+E_mic(i+1,:)*hz; %����΢�������ĵ糡
        
        B_Li3(p,3)=1e-4*(2595.7+32.4092*(1000*pos_Li1(p))-3.1038*(1000*pos_Li1(p))^2-...   %��ֵ�õ�����λ�ô��ĴŸ�Ӧǿ��
            0.0130847*(1000*pos_Li1(p))^3+0.00198809*(1000*pos_Li1(p))^4-2.48996e-5*(1000*pos_Li1(p))^5+8.9815e-8*(1000*pos_Li1(p))^6);
        if fi<=1.5                                           %�ĸ�if����ȷ������λ�ã��Ӷ�Ϊ֮����΢�������ĴŸ�Ӧǿ��
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
        
        %����BORIS���������һʱ�̵��ٶȺ�λ��
        vel_Li3(p,:)=UpdateVelocity(E_Li3(p,:),B_Li3(p,:),vel_Li3(p,:),dt_i); %�����ٶ�
        pos_Li3(p,1)=pos_Li3(p,1)+vel_Li3(p,3)*dt_i;                          %����λ��
        
        %�������ձ߽�������Li3+����໥���ú�������
        if pos_Li3(p,1)<0             %��߽�
            pos_Li3(p,1)=pos_Li3(N_Li3,1); %�����һ��Li3+��λ���滻������λ�ô�
            vel_Li3(p,:)=vel_Li3(N_Li3,:); %�����һ��Li3+���ٶ��滻������λ�ô�
            N_Li3=N_Li3-1;               %Li3+������һ��
        end
        if pos_Li3(p,1)>z_max          %�ұ߽�
            pos_Li3(p,1)=pos_Li2(N_Li3,1); %�����һ��Li3+��λ���滻������λ�ô�
            vel_Li2(p,:)=vel_Li2(N_Li3,:); %�����һ��Li3+���ٶ��滻������λ�ô�
            N_Li3=N_Li3-1;               %Li3+������һ��
        end
        
        p=p+1;
    end
end

filename=strcat(datestr(clock,'yy-mm-dd-HH-MM-SS'),'E_mic.txt');
fid=fopen(filename,'wt');       %��΢���糡����д��E_mic.txt�У����зֱ���ʱEx,Ey,Ez�������������nz�����
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
fid=fopen(filename,'wt');       %��΢���ų�����д��B_mic.txt�У����зֱ���ʱBx,By,Bz�������������nz�����
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
fid=fopen(filename,'wt');       %���������������Ĵų�д��B_ex_z.txt�ļ��У���һ��Ϊÿ������z���꣬�ڶ���Ϊ��Ӧ�����Ÿ�Ӧǿ��
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

%���΢����ų��ֲ���ͼ��
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

%���z-vͼ��
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


