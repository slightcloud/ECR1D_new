%Diagnostic of simulation results

% z-v of electrons
figure(1);
subplot(3,1,1);
scatter(pos_e(1:N_e,1),vel_e(1:N_e,1));   %z-vx
subplot(3,1,2);
scatter(pos_e(1:N_e,1),vel_e(1:N_e,2));   %z-vy
subplot(3,1,3);
scatter(pos_e(1:N_e,1),vel_e(1:N_e,3));   %z-vz
suptitle('z-v of electrons');

% z-v of Li+
figure(2);
subplot(3,1,1);
scatter(pos_Li1(1:N_Li1,1),vel_Li1(1:N_Li1,1));   %z-vx
subplot(3,1,2);
scatter(pos_Li1(1:N_Li1,1),vel_Li1(1:N_Li1,2));   %z-vy
subplot(3,1,3);
scatter(pos_Li1(1:N_Li1,1),vel_Li1(1:N_Li1,3));   %z-vz
suptitle('z-v of Li+');

% z-N_p

figure(3);
scatter((1:nz-1)*dz,count_e(1:nz-1,1));
hold on
scatter((1:nz-1)*dz,count_Li(1:nz-1,1)/2);
scatter((1:nz-1)*dz,count_Li(1:nz-1,2));
legend('e','Li+','Li2+');
hold off
title('z-particle number');


% microwave

figure(4);
subplot(2,2,1);
plot([0,(1:nz-1)*dz],E_mic(:,1));
subplot(2,2,2);
plot([0,(1:nz-1)*dz],E_mic(:,2));
subplot(2,2,3);
plot([0,(1:nz-1)*dz],B_mic(:,1));
subplot(2,2,4);
plot([0,(1:nz-1)*dz],B_mic(:,2));
suptitle('z-microwave');

% Potential

figure(5);
plot([0,(1:nz-1)*dz],phi(:,1));
title('z-potential');

% vz-vx vz-vy vx-vy

figure(6);
subplot(3,1,1);
scatter(vel_e(:,3),vel_e(:,1));
title('vz-vx');
subplot(3,1,2);
scatter(vel_e(:,3),vel_e(:,2));
title('vz-vy');
subplot(3,1,3);
scatter(vel_e(:,1),vel_e(:,2));
title('vx-vy');

% Export


% filename=strcat(datestr(clock,'yy-mm-dd-HH-MM-SS'),'phi.txt');
% fid=fopen(filename,'wt');       %将电势数据写到phi.txt中
% 
% row = size(phi,1);
% for i=1:1:row
%     for j=1:1:2
%         if (j==1)
%             fprintf(fid,'%f\t',(i-1)*dz);
%         else
%             fprintf(fid,'%f\n',phi(i,j-1));
%         end
%     end
% end
% fclose(fid);


% Energy distribution function of electrons


ve_all=sqrt(sum((vel_e(1:N_e,:))'.^2));    %电子速度的模
Ek_e_all=m_e*c^2*(1./sqrt(1-ve_all.^2/c^2)-1)/q_e; 
figure(7);
scatter(pos_e(1:N_e,1),Ek_e_all);         %不同位置处的电子能量
title('z-Ek_e');
xlabel('z(m)');
ylabel('Ek_e(eV)');

Ek_e_min = min(Ek_e_all);
Ek_e_max = max(Ek_e_all);
Ek_e_dis = linspace(Ek_e_min,Ek_e_max,100);
sta_e = hist(Ek_e_all,Ek_e_dis);     %每个区间的个数，直方图
sta_e_fra = sta_e/length(Ek_e_all);  %每个区间电子个数占总个数的比例
figure(8);
bar(Ek_e_dis,sta_e_fra);             %概率密度分布图
hold on
% Ek_e_fit = polyfit(Ek_e_dis,sta_e_fra,6);
% sta_e_fra_fit = polyval(Ek_e_fit,Ek_e_dis);
Ek_e_dis_smooth = linspace(Ek_e_min,Ek_e_max,1000);
sta_e_fra_smooth = interp1(Ek_e_dis,sta_e_fra,Ek_e_dis_smooth,'spline');
plot(Ek_e_dis_smooth,sta_e_fra_smooth);
title('Ek_e-fraction');
xlabel('Ek_e(eV)');
ylabel('Fraction');
hold off;

% Energy distribution function of Li+

vLi1_all=sqrt(sum((vel_Li1(1:N_Li1,:))'.^2));    %Li+速度的模
Ek_Li1_all=m_Li*c^2*(1./sqrt(1-vLi1_all.^2/c^2)-1)/q_e; 
figure(9);
scatter(pos_Li1(1:N_e,1),Ek_Li1_all);         %不同位置处的Li+能量
title('z-Ek_{Li^+}');
xlabel('z(m)');
ylabel('Ek_{Li^+}(eV)');

Ek_Li1_min = min(Ek_Li1_all);
Ek_Li1_max = max(Ek_Li1_all);
Ek_Li1_dis = linspace(Ek_Li1_min,Ek_Li1_max,100);
sta_Li1 = hist(Ek_Li1_all,Ek_Li1_dis);     %每个区间的个数，直方图
sta_Li1_fra = sta_Li1/length(Ek_Li1_all);  %每个区间Li+个数占总个数的比例
figure(10);
bar(Ek_Li1_dis,sta_Li1_fra);             %概率密度分布图
hold on
% Ek_e_fit = polyfit(Ek_e_dis,sta_e_fra,6);
% sta_e_fra_fit = polyval(Ek_e_fit,Ek_e_dis);
Ek_Li1_dis_smooth = linspace(Ek_Li1_min,Ek_Li1_max,1000);
sta_Li1_fra_smooth = interp1(Ek_Li1_dis,sta_Li1_fra,Ek_Li1_dis_smooth,'spline');
plot(Ek_Li1_dis_smooth,sta_Li1_fra_smooth);
title('Ek_{Li^+}-fraction');
xlabel('Ek_{Li^+}(eV)');
ylabel('Fraction');
hold off;


%Energy distribution at different position

%------------------------z=0.045~0.046---------------------
[row_e,~] = find(pos_e>=0.045 & pos_e<=0.04501);
row_e_n = size(row_e);             %个数
ve_all=sqrt(sum((vel_e(row_e(1):row_e(end),:))'.^2));    %电子速度的模
Ek_e_all=m_e*c^2*(1./sqrt(1-ve_all.^2/c^2)-1)/q_e;
Ek_e_min = min(Ek_e_all);
Ek_e_max = max(Ek_e_all);
Ek_e_dis = linspace(Ek_e_min,Ek_e_max,100);
sta_e = hist(Ek_e_all,Ek_e_dis);     %每个区间的个数，直方图
sta_e_fra = sta_e/length(Ek_e_all);  %每个区间电子个数占总个数的比例
figure(11);
bar(Ek_e_dis,sta_e_fra);             %概率密度分布图
hold on
% Ek_e_fit = polyfit(Ek_e_dis,sta_e_fra,6);
% sta_e_fra_fit = polyval(Ek_e_fit,Ek_e_dis);
Ek_e_dis_smooth = linspace(Ek_e_min,Ek_e_max,1000);
sta_e_fra_smooth = interp1(Ek_e_dis,sta_e_fra,Ek_e_dis_smooth,'spline');
plot(Ek_e_dis_smooth,sta_e_fra_smooth);
title('Ek_e-fraction at z = 0.0455m');
xlabel('Ek_e(eV)');
ylabel('Fraction');
hold off;

%------------------------z=0.080~0.082---------------------
[row_e,~] = find(pos_e>=0.04 & pos_e<=0.0401);
row_e_n = size(row_e);             %个数
ve_all=sqrt(sum((vel_e(row_e(1):row_e(end),:))'.^2));    %电子速度的模
Ek_e_all=m_e*c^2*(1./sqrt(1-ve_all.^2/c^2)-1)/q_e;
Ek_e_min = min(Ek_e_all);
Ek_e_max = max(Ek_e_all);
Ek_e_dis = linspace(Ek_e_min,Ek_e_max,100);
sta_e = hist(Ek_e_all,Ek_e_dis);     %每个区间的个数，直方图
sta_e_fra = sta_e/length(Ek_e_all);  %每个区间电子个数占总个数的比例
figure(12);
bar(Ek_e_dis,sta_e_fra);             %概率密度分布图
hold on
% Ek_e_fit = polyfit(Ek_e_dis,sta_e_fra,6);
% sta_e_fra_fit = polyval(Ek_e_fit,Ek_e_dis);
Ek_e_dis_smooth = linspace(Ek_e_min,Ek_e_max,1000);
sta_e_fra_smooth = interp1(Ek_e_dis,sta_e_fra,Ek_e_dis_smooth,'spline');
plot(Ek_e_dis_smooth,sta_e_fra_smooth);
title('Ek_e-fraction at z = 0.08m');
xlabel('Ek_e(eV)');
ylabel('Fraction');
hold off;