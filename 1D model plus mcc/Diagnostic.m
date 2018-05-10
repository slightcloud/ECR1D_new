%Diagnostic of simulation results

% z-v of electrons
figure(1);
subplot(3,1,1);
scatter(pos_e(1:N_e,1),vel_e(1:N_e,1));   %z-vx
subplot(3,1,2);
scatter(pos_e(1:N_e,1),vel_e(1:N_e,2));   %z-vy
subplot(3,1,3);
scatter(pos_e(1:N_e,1),vel_e(1:N_e,3));   %z-vz

% z-v of Li+
figure(2);
subplot(3,1,1);
scatter(pos_Li1(1:N_Li1,1),vel_Li1(1:N_Li1,1));   %z-vx
subplot(3,1,2);
scatter(pos_Li1(1:N_Li1,1),vel_Li1(1:N_Li1,2));   %z-vy
subplot(3,1,3);
scatter(pos_Li1(1:N_Li1,1),vel_Li1(1:N_Li1,3));   %z-vz

% z-N_p

figure(3);
scatter((100:nz-1)*dz,count(100:nz-1,1));
hold on
scatter((100:nz-1)*dz,count_Li(100:nz-1,1));
scatter((100:nz-1)*dz,count_Li(100:nz-1,2));
legend('e','Li+','Li2+');
hold off

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

% Potential

figure(5);
plot([0,(1:nz-1)*dz],phi(:,1));
