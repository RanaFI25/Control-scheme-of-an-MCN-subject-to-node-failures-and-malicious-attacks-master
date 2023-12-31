A1 = [0.8 1 0 0; 0 0.9 1 0; 0 0 0.8 1; 0 0 0 0.9];
B1 = [0.1 0.1; 0.1 0.1; 0.1 0.1; 0.1 0.1];
C1 = eye(4);

A2 = [0.7 1 0 0; 0 0.8 1 0; 0 0 0.7 1; 0 0 0 0.8];
B2 = [0.2 0.2; 0.2 0.2; 0.2 0.2; 0.2 0.2];
C2 = eye(4);

A3 = [0.6 1 0 0; 0 0.7 1 0; 0 0 0.6 1; 0 0 0 0.7];
B3 = [0.3 0.3; 0.3 0.3; 0.3 0.3; 0.3 0.3];
C3 = eye(4);

%x_hat[k+1] = A * x_hat[k] + B * u[k] + L * (y[k] - C * x_hat[k])

%y_hat[k] = C * x_hat[k]
x = [1;1;1;1]; % Initial state
y1 = zeros(100,4);
x_hat = x; % Initial estimated state
for k = 1:100
    x = A1*x;
    y1(k,:) = C1*x;
    x_hat = A1*x_hat ;
    y_hat = C1*x_hat;
    x_hat = x_hat + 0.01 * (y1(k,:) - y_hat); % Luenberger observer
end

x = y1(100,:)';
y2 = zeros(100,4);
x_hat = x; % Initial estimated state
for k = 100:200
    x = A2*x;
    y2(k-99,:) = C2*x;
    x_hat = A2*x_hat ;
    y_hat = C2*x_hat;
    x_hat = x_hat + 0.01 * (y2(k-99,:)' - y_hat); % Luenberger observer
end

x = y2(100,:)';
y3 = zeros(100,4);
x_hat = x; % Initial estimated state
for k = 200:300
    x = A3*x ;
    y3(k-199,:) = C3*x;
    x_hat = A3*x_hat;
    y_hat = C3*x_hat;
    x_hat = x_hat + 0.01 * (y3(k-199,:)' - y_hat); % Luenberger observer
end

t = linspace(0,100,100);
figure;
plot(t,y1,t,y_hat);
xlabel('Time');
ylabel('Output');
title('Model 1 Output (Measured vs. Estimated)');
legend('Measured','Estimated');

t = linspace(100,200,100);
figure;
plot(t,y2,t,y_hat);
xlabel('Time');
ylabel('Output');
title('Model 2 Output (Measured vs. Estimated)');
legend('Measured','Estimated');

t = linspace(200,300,100);
figure;
plot(t,y3,t,y_hat);
xlabel('Time');
ylabel('Output');
title('Model 3 Output (Measured vs. Estimated)');
legend('Measured','Estimated');
