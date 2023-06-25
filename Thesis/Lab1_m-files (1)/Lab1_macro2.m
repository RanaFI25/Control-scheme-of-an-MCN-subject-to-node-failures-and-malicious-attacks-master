% Script that simulates the intersample behaviour for a continuous-time
% system controlled by sampling controllers.

%h=0.1;

Tend=8; % End time for simulation.
Lab1_macro1; % Run script to define systems, parameters etc.
[Ycont,Tcont]=step(Gc_c,Tend); % Step response for cont.-time controller.

% Simulate step responses for cont.-time system with sampling controllers.
[Yct,Tct,Ydt,Tdt]=intersample_sim(sys_c,Fr_tus,Fy_tus,Tend);
[Ycz,Tcz,Ydz,Tdz]=intersample_sim(sys_c,Fr_zoh,Fy_zoh,Tend);

% Plot the step responses.
figure(2)
plot(Tcont,Ycont,Tct,Yct,Tcz,Ycz)
title('Step response for cont.-time system with sampling controller')
xlabel('Time (seconds)')
ylabel('Output')
hold on
legend('c-t contr','Tustin contr','ZOH contr','Location','southeast');
ax = gca;
ax.ColorOrderIndex = 2;
plot(Tdt,Ydt,'o',Tdz,Ydz,'o')
hold off