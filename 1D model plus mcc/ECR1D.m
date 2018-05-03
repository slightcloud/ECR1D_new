%1 2/2D ECR plasma simulation using PIC-MCC
clear variables
clc

global q_e m_e

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
lamda=sqrt(EPS0*Te_average/(n0*q_e));    %�°ݳ��ȣ���λ��m
z_min=0;                                 %z�������Сֵ����λ��m
L_real=0.1035;                             %����Դz��������ֵ����λ��m
dz=lamda;                                %���ÿռ���
dt=dz/(2*c);                             %����ʱ������Ϊ��֤ģ�⾫ȷ�ȣ���֤dt<=dz/c��һ����Ϊdz/(2*c)��λ��s
nz=round(L_real/dz)+1;                   %�����
z_max=dz*(nz-1);                         %ģ������z��������ֵ
spwt=1e8;                                %ÿ�������Ӱ�������ʵ������
vth_e=sqrt(2*q_e*0.5/m_e);               %���ӵ����ٶȣ�������ӵ�����Ϊ0.5 eV����λ��m/s
N_e=3000;                                %����N_e���������ڳ������
tol=0.001;                               %�������ľ���ֵ
max_part=3000;                           %�������������ֵ
%B_ex_z=0.0875;                             %�ⲿ�����������ĴŸ�Ӧǿ�ȣ���λ��T
E0=1.5e5;                                %΢���糡ǿ�ȵķ�ֵ����λ��V/m����ȷ����
B0=0.0001;                               %΢���Ÿ�Ӧǿ�ȵķ�ֵ����λ:T����ȷ����
w=2.45;                                  %΢����Ƶ��,��λ:GHz
step_num=25000;                          %�ܹ����е�ʱ�䲽��
N_phi=5e4;                               %�������ĵ�������
z=zeros(nz,1);                           %����ʢ��ÿ����������
z(2:nz,1)=(1:nz-1)'*dz*1000;                  %ÿ����������ֵ����λ��mm
B_ex_z=0.0001*(2595.7+32.4092*z-3.1038*z.^2-0.0130847*z.^3+0.00198809*z.^4-2.48996e-5*z.^5+8.9815e-8*z.^6);   %����ʵ���ϲ�õĴŸ�Ӧǿ��ֵ������ϣ���λ��T
%plot(z,B_ex_z);



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
pos=zeros(max_part,1);         %���ӵ�λ����Ϣ
vel=zeros(max_part,3);         %���ӵ��ٶ���Ϣ



%������ʼ���ӵ�λ�ú��ٶ�
pos(1:N_e)=rand(N_e,1)*z_max;                                                   %��ʼ���ӵ�λ����0��z_max֮��

%vel(1:N_e/2,:)=vth_e*rand(N_e/2,3);
%vel(N_e/2+1:N_e,:)=-vth_e*rand(N_e/2,3);

%vel(1:N_e,:)=1e8;
vel(1:N_e,:)=vth_e*2*(rand(N_e,3)+rand(N_e,3)+rand(N_e,3)-1.5);                 %���ó�ʼ���ӵ��ٶ������ٶȵ�-3����3��֮��
vel(1:N_e,:)=UpdateVelocity(E_e(1:N_e,:),B_e(1:N_e,:),vel(1:N_e,:),-0.5*dt);    %�����ӵ��ٶ���ǰ�ƶ����ʱ�䲽��

%��ʼ��ѭ��������ʱ���ѭ��
for ts=1:step_num                                                               %���е�ʱ�䲽��
    
    chg=zeros(nz,1);                                                            %���ڴ��ÿ�������䵽�ĵ����
    den=zeros(nz,1);                                                            %ÿ�����ĵ���ܶȣ���λ��C/m
    den_vel=zeros(nz,3);                                                        %���ڴ��ÿ������ϵĵ����ܶȷ���
    B_mic_old=B_mic;                                                            %B_mic_old���ں�B_mic��ƽ�����Ӷ���������ٶ�
    pos_half=zeros(N_e,1);                                                      %���ڴ�Ű�ʱ��ʱ������λ��
    
    for p_half=1:N_e                                                            %���������ӵ�λ����ǰ�ƶ�����ʱ�̴�
        pos_half(p_half,1)=pos(p_half,1)-vel(p_half,3)*0.5*dt;
        if pos_half(p_half,1)<0                                                 %����if������ʾ�����Ա߽�����
            pos_half(p_half,1)=pos_half(p_half,1)+z_max;
        end
        if pos_half(p_half,1)>z_max
            pos_half(p_half,1)=pos_half(p_half,1)-z_max;
        end
    end
    % pos_half=pos-vel(:,3)*0.5*dt;                                               %�����ӵ�λ����ǰ�ƶ����ʱ�䲽��������ð�ʱ�̴��ĵ����ܶ�
    
    %���Ƚ����е��ӵĵ�ɰ���������ԭ�����
    
    count=zeros(nz-1,1);            %���ڼ���ÿ����Ԫ���ڵ�������
    count_g=zeros(nz,1);            %���ڼ���ÿ������ϵ��������ܶ�
    
    for p=1:N_e
        fi=1+pos(p)/dz;            %ʵ�ʸ��λ�ã�Ϊ������
        i=floor(fi);               %��Ӧ�������λ��
        hz=fi-i;                   %�������i�����֮��ľ���ռ���벽���ı���
        
        %�����ӵ�ɰ��ձ������䵽�ٽ������������
        chg(i)=chg(i)+(1-hz);      %���䵽��i������ϵĵ�ɱ���
        chg(i+1)=chg(i+1)+hz;      %���䵽��i+1������ϵĵ�ɱ���
        
        fi_half=1+pos_half(p)/dz;  %��ʱ�̴��ĸ��λ��
        i_half=floor(fi_half);     %��ʱ�̴����λ�ö�Ӧ���������λ��
        hz_half=fi_half-i_half;    %��ʱ�̴��������i_half�����֮��ľ���ռ���벽���ı���
        
        count(i_half)=count(i_half)+1;  %ͳ��ÿ����Ԫ���ڵ�������Ŀ
        
        
        
        %���ٶȰ��������䵽�ٽ������������
        den_vel(i_half,:)=den_vel(i_half,:)+spwt*q_e*vel(p,:)*(1-hz_half);  %���䵽��i������ϵ��ٶ�
        den_vel(i_half+1,:)=den_vel(i_half+1,:)+spwt*q_e*vel(p,:)*hz_half;  %���䵽��i+1������ϵ��ٶ�
    end
    
    den=spwt*q_e*chg/dz;       %ÿ������ϵĵ���ܶ�
    %J=spwt*q_e*(-1)*den_vel;  %ÿ������ϵĵ����ܶ�
    %J=den_vel.*den;
    %J=1e9*den_vel;
    
    %����ÿ������ϵ������ܶ�
    for count_i=2:nz-1
        count_g(count_i)=0.5*(count(count_i-1)+count(count_i))/(2*dz);
    end
    count_g(1)=count(1)/dz;      %��һ�����������ܶ�
    count_g(nz)=count(nz-1)/dz;  %���һ�����������ܶ�
    
    J=1e5*den_vel.*count_g;          %����ÿ������ϵĵ����ܶ�
    
    %Ӧ�ñ߽��������ڱ߽��Ͻ���һ��ĳ����й��ף���˱߽��ϵ��ܶ�Ӧ����2
    den(1)=2*den(1);   %��߽�
    den(nz)=2*den(nz); %�ұ߽�
    
    %��ʼ��Ⲵ�ɷ��̣��õ�����ϵĵ��ƣ����÷���Ϊ��˹-���¶���������G-S��
    % phi_tem=phi;
    for j=1:N_phi     %��������
        
        for i=2:nz-1  %�м���ĵ���
            phi_tem(i)=0.5*(den(i)/EPS0*dz*dz+phi(i-1)+phi(i+1));   %�������޲�ַ����õ�
        end
        
        phi_tem(1)=phi(2);     %��ڴ�Ϊ��һ��߽�����
        phi_tem(nz)=phi(nz-1); %���ڴ�ҲΪ��һ��߽�����
        
        if mod(j,10)==0
            R=norm(phi_tem-phi); %������ε����Ĳв�
            if (R<=tol)          %��������
                phi=phi_tem;
                %fprintf('Congratulations! The solver of the electric potential is converged at %d step\n',j);
                break;
            end
        end
        phi=phi_tem;
    end
    if j>=N_phi                  %���������ⲻ����������ֹ���������
        fprintf('The result of the electric potential is not convergent, the program is stopped!\n');
        break;
    end
    
    %���ÿ������ϵľ��糡��������һά�����ֻ�ܵõ�z����ĵ糡ǿ��
    E_sta(2:nz-1)=(phi(1:nz-2)-phi(3:nz))/2/dz;      %�м���ľ��糡
    E_sta(1)=(phi(1)-phi(2))/dz;                     %��ڴ����糡
    E_sta(nz)=(phi(nz-1)-phi(nz))/dz;                %���ڴ����糡
    
    %����΢����ų��ĸ�������
    E_mic_xold=E0*cos(w*1e9*2*pi()*(ts-2)*dt);
    E_mic_yold=E0*sin(w*1e9*2*pi()*(ts-2)*dt);
    E_mic(1,1)=E0*cos(w*1e9*2*pi()*(ts-1)*dt);              %��ڴ�x����ĵ糡ǿ��E_x=E0cos(wt)
    % B_mic_xin_full=B0*sin(w*1e9*2*pi()*(ts-0.5)*dt);       %��ڴ�x����Ÿ�Ӧǿ�ȣ�ʱ��Ϊ��
    B_mic_xin=B0*sin(w*1e9*2*pi()*(ts-0.5)*dt);             %��ڴ�x����ĴŸ�Ӧǿ��ΪB_x=B0sin(wt)����ʱ��������ƶ���t/2
    E_mic(1,2)=E0*sin(w*1e9*2*pi()*(ts-1)*dt);              %��ڴ�y����ĵ糡ǿ��E_y=E0sin(wt)
    %B_mic_yin_full=B0*cos(w*1e9*2*pi()*(ts-0.5)*dt);       %��ڴ�y����Ÿ�Ӧǿ�ȣ�ʱ��Ϊ��
    B_mic_yin=B0*cos(w*1e9*2*pi()*(ts-0.5)*dt);             %��ڴ�y����ĴŸ�Ӧǿ��B_y=B0cos(wt)����ʱ��������ƶ���t/2
    %B_mic(1,1)=B_mic_xin+dt/(2*dz)*(E_mic(2,2)-E_mic(1,2));                   %��Bx�ڿռ��ʱ���϶�����ƶ���z/2�ͦ�t/2
    B_mic(1,1)=B_mic_xin+dz/2*(1/(c^2)*(E_mic(1,2)-E_mic_yold)/dt+MU0*J(1,2)); %��ڴ���Bx����Bx�ڿռ�������ƶ���z/2�������
    B_mic(1,2)=B_mic_yin-dz/2*(1/(c^2)*(E_mic(1,1)-E_mic_xold)/dt+MU0*J(1,1)); %��ڴ���By����By�ڿռ�������ƶ���z/2�������
    %B_mic(1,2)=B_mic_yin-dt/(2*dz)*(E_mic(2,1)-E_mic(1,1));                   %��By�ڿռ��ʱ���϶�����ƶ���z/2�ͦ�t/2
    
    % E_mic(nz,1)=0; %���ڴ�x����ĵ糡ǿ��0
    % E_mic(nz,2)=0; %���ڴ�y����ĵ糡ǿ��0
    
    B_mic_xout=B_mic(nz,1);
    B_mic_yout=B_mic(nz,2);
    
    %B_mic_xout=0;  %���ڴ�x����ĴŸ�Ӧǿ��0
    % B_mic_yout=0; %���ڴ�y����ĴŸ�Ӧǿ��0
    % B_mic(nz-1,1)=B_mic_xout-dz/2*MU0*J(nz,2); %���ڴ���Bx����ǰ���
    % B_mic(nz-1,2)=B_mic_yout+dz/2*MU0*J(nz,1); %���ڴ���By����ǰ���
    
    B_mic(2:nz-1,1)=B_mic(2:nz-1,1)+dt/dz*(E_mic(3:nz,2)-E_mic(2:nz-1,2));   %Bx
    B_mic(2:nz-1,2)=B_mic(2:nz-1,2)-dt/dz*(E_mic(3:nz,1)-E_mic(2:nz-1,1));   %By
    E_mic(2:nz-1,1)=E_mic(2:nz-1,1)-c^2*dt*(MU0*J(2:nz-1,1)+(B_mic(2:nz-1,2)-B_mic(1:nz-2,2))/dz);    %Ex
    E_mic(2:nz-1,2)=E_mic(2:nz-1,2)+c^2*dt*(-1*MU0*J(2:nz-1,2)+(B_mic(2:nz-1,1)-B_mic(1:nz-2,1))/dz); %Ey
    %E_mic(1:nz,3)=E_mic(1:nz,3)-c^2*MU0*dt*J(1:nz,3); %Ez
    
    E_boun(ts,:)=E_mic(nz-1,:);
    B_boun(ts,:)=B_mic(nz-1,:);
    if ts>2
        E_mic(nz,:)=E_boun(ts-2,:);
        B_mic(nz,:)=B_boun(ts-2,:);
    end
    
    %������ϵľ��糡�͵�ų�ȫ�����䵽ÿ�������ϣ����⻹Ҫ�����ⲿ����������Ĵų�
    E_e=zeros(max_part,3);           %ÿ����������λ�ô��ĵ糡ǿ��
    B_e=zeros(max_part,3);           %ÿ����������λ�ô��ĴŸ�Ӧǿ��
    %B_e(:,3)=B_e(:,3)+B_ex_z;        %z����ĴŸ�Ӧǿ����Ҫ�����ⲿ�Ÿ�Ӧǿ��
    B_mic_ave=0.5*(B_mic_old+B_mic); %����ƽ������������ʱ���ֵƽ��������ʱ����
    p=1;                             %���ڱ�����������
    while p<=N_e
        fi=1+pos(p)/dz;  %ʵ�ʸ��λ�ã�Ϊ������
        i=floor(fi);     %��Ӧ�������λ��
        hz=fi-i;         %�������i�����֮��ľ���ռ���벽���ı���
        E_e(p,3)=E_e(p,3)+E_sta(i)*(1-hz)+E_sta(i+1)*hz;     %��������֮��ľ��糡
        E_e(p,:)=E_e(p,:)+E_mic(i,:)*(1-hz)+E_mic(i+1,:)*hz; %����΢�������ĵ糡
        %B_e(p,3)=B_e(p,3)+B_ex_z(i)*(1-hz)+B_ex_z(i+1)*hz;   %������ϵľ��ų������ÿ������
        B_e(p,3)=1e-4*(2595.7+32.4092*(1000*pos(p))-3.1038*(1000*pos(p))^2-...   %��ֵ�õ�����λ�ô��ĴŸ�Ӧǿ��
            0.0130847*(1000*pos(p))^3+0.00198809*(1000*pos(p))^4-2.48996e-5*(1000*pos(p))^5+8.9815e-8*(1000*pos(p))^6);
        if fi<=1.5                                           %�ĸ�if����ȷ������λ�ã��Ӷ�Ϊ֮����Ÿ�Ӧǿ��
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
        vel(p,:)=UpdateVelocity(E_e(p,:),B_e(p,:),vel(p,:),dt); %�����ٶ�
        pos(p,1)=pos(p,1)+vel(p,3)*dt;                          %����λ��
        
        %���������Ա߽�����
        if pos(p,1)<0
            pos(p,1)=pos(p,1)+z_max;
        end
        if pos(p,1)>z_max
            pos(p,1)=pos(p,1)-z_max;
        end
        
        p=p+1;
    end
    
    fprintf('ts: %d\n', ts);
    % plot(1:nz,E_mic(:,1));
    % pause(0.1);
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
            fprintf(fid,'%f\t',pos(i));
        end
        if j>1&&j<4
            fprintf(fid,'%g\t',vel(i,j-1));
        end
        if j==4
            fprintf(fid,'%g\n',vel(i,j-1));
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
scatter(pos(1:N_e),vel(1:N_e,1));
xlabel('z(m)');
ylabel('v_x(m/s)');
subplot(3,1,2);
scatter(pos(1:N_e),vel(1:N_e,2));
xlabel('z(m)');
ylabel('v_y(m/s)');
subplot(3,1,3);
scatter(pos(1:N_e),vel(1:N_e,3));
xlabel('z(m)');
ylabel('v_z(m/s)');


