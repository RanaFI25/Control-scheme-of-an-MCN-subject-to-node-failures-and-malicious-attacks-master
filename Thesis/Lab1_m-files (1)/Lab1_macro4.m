% Script that simulates the intersample behaviour for a continuous-time
% system controlled by sampling controllers.

%h=0.1;

Tend=8; % End time for simulation.
Lab1_macro3; % Run script to define systems, parameters etc.

% Simulate step responses for cont.-time system with sampling controllers.
[Ycdb,Tcdb,Ydb,Tdb]=intersample_sim(sys_c,Fr_db,Fy_db,Tend);

%[Ycdb,Tcdb,Ydb,Tdb,Ucdb]=intersample_sim(sys_c,Fr_db,Fy_db,Tend);


% Plot the step responses.
figure(4)
plot(Tcdb,Ycdb)
title('Step response for cont.-time system with sampling controller')
xlabel('Time (seconds)')
ylabel('Output')
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(Tdb,Ydb,'o')
legend('c-t output','d-t output','Location','southeast');
hold off

%
% figure(4)
% plot(Tcdb,Ycdb)
% title('Step response for cont.-time system with sampling controller')
% xlabel('Time (seconds)')
% ylabel('Output')
% hold on
% ax = gca;
% ax.ColorOrderIndex = 1;
% plot(Tdb,Ydb,'o')
% plot(Tcdb,Ucdb)
% legend('c-t output','d-t output','control input','Location','southeast');
% hold off