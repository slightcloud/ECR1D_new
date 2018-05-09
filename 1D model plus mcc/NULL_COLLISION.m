%该子函数用来计算每个单元格内的最大碰撞频率
function nu_t = NULL_COLLISION(n0,n1,n2)
%n0代表的是中性原子在每个单元格内的数密度，一般为一个定值
%n1代表的是Li+在每个单元格内的数密度，随位置和时间变化
%n2代表的是Li2+在每个单元格内的数密度，随位置和时间变化

global m_e q_e c

load sigma_Li0_Li1.txt;   %载入Li0->Li+的截面数据
load sigma_Li1_Li2.txt;   %载入Li+->Li2+的截面数据
load sigma_Li2_Li3.txt;   %载入Li2+->Li3+的截面数据
row = size(n1,1);         %计算密度矩阵的行数row
nu_t = zeros(row,1);        %给nu_t开辟空间
Te = (5:0.1:1500)';          %电子温度的取值范围，eV
ve = 100*sqrt(c^2*(1-(1+Te*q_e/(m_e*c^2)).^(-2)));  %根据相对论运动学，将电子的动能转换为速度，单位：cm/s
sigma = zeros(size(Te,1),3); %该矩阵用来盛放三个电离反应的截面，每一列是一个反应
nu = zeros(size(Te,1),3);    %该矩阵用来盛放
for i = 1:3
    eval(['A = sigma_Li',num2str(i-1),'_Li',num2str(i),'(:,1);']);   %A为外部文件中提取出的电子密度
    eval(['B = sigma_Li',num2str(i-1),'_Li',num2str(i),'(:,2);']);   %B为外部文件中提取出的反应截面
    sigma(:,i)=interp1(A,B,Te,'pichip');                            %通过插值的方法得到每个电子温度对应的反应截面，单位：cm^2
    if i == 2                                                       %将小于Li+->Li2+反应阈能的截面设置为0
        for j = 1:size(Te,1)
            if Te(j,1) <= 76
                sigma(j,i) = 0;
            end
        end
    end
    if i == 3                                                      %将小于Li2+->Li3+反应阈能的截面设置为0
        for j = 1:size(Te,1)
            if Te(j,1) <= 118
                sigma(j,i) = 0;
            end
        end
    end
end

%开始计算每个单元格内最大碰撞频率

for i = 1:row                               %遍历每个单元格
    nu(:,1) = n0*sigma(:,1).*ve;            %电子与中性原子在第i个单元格内的碰撞频率
    nu(:,2) = n1(i,1)*sigma(:,2).*ve;       %电子与Li+在第i个单元格内的碰撞频率
    nu(:,3) = n2(i,1)*sigma(:,3).*ve;       %电子与Li2+在第i个单元格内的碰撞频率
    nu_tem = nu(:,1)+nu(:,2)+nu(:,3);       %三个碰撞频率相加，得到第i个单元格内总的碰撞频率分布
    nu_t(i,1) = max(nu_tem);                %求出第i格单元格内的最大碰撞频率
end

