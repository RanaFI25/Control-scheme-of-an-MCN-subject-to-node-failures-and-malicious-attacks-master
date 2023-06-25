clear all
clc
var1 = 0.0146; 
var2 = 0.0253;
var3 = 0.0022;

%q = 1% model 1

A = [-var1 0 0 0; var1 -var1 0 0;0 0 -var2 0; 0 0 var2 -var2 ];
B = transpose([var3 0 0 0;0 0 0 0]); 
C = eye(4);
D = 0;

L = place(A',C',[.5 .7 .5 .7])';

%q = 2% model 2 

A2 = [-var1 0 0 0; var1 -var1 0 0;0 0 -var2 0; 0 0 var2 -var2 ];
B2 = transpose([var3 0 0 0;0 0 var3 0]); 
C2 = eye(4);
D2= 0;



L2 = [7.485 0 0 0; 0.014 0.485 0 0; 0 0 6.477 0;0 0 0.022 0.877];

%q = 3% model 3

A3 = [-0.0976 0 0.0801 0; 0.0175 -0.0175 0 0;0.0801 0 -0.0981 0; 0 0 0.0179 -0.0179];
B3 = transpose([var3 0 0 0;0 0 var3 0]); 
C3 = eye(4);
D3 = zeros(4,2);



L3 = [7.402 0 0.08 0;0.017 0.482 0 0;0.08 0 6.401 0;0 0 0.017 0.882];

%Observer % 


%Luenberger observer example%
% A = [1.8097 -0.8187; 1 0]; 
% B = [0.5; 0]; 
% C = [0.1810 -0.1810]; 
% D = 0; 

% Ob = obsv(A,C);
% unobsv = length(A) - rank(Ob)

% L = place(A',C',[.5 .7])';
% L2 = place(A',C',[.75 .8])';
% L3 = place(A',C',[.4 .5])';

% eig(A+L*C)

% x=[1;1;1;1]; % initial state
% xhat=[1.01;1.01;1.01;1.01]; % initial estimate
% XX=x;
% XXhat=xhat;
% T=100;
% UU=.1*ones(1,T); % input signal
% uvector = [1;1];
% Bu= B*uvector;
% 
% for k=0:T-1,
%     u=UU(k+1);
%     y=C*x; 
%     yhat=C*xhat; 
%     x=A*x+B*u; 
%     %xhat=A*xhat+B*u+L*(y-yhat);
%     xhat = (A-L*C)*xhat+B*u+L*y;
%     XX=[XX,x]; % measured output
%     XXhat=[XXhat,xhat]; % estimated output
% end
% 
% 
% 
% plot(XX(1,:), LineStyle='-.')
% hold on
% plot(XXhat(1,:), LineStyle='--'); %x4 pure  measurement and estimate 
% hold on
% plot([XX(1,:)-XXhat(1,:)],LineStyle=':') % residual (x4 - x^4)

%Measured Output of model 1 in 100 time step x0 = 1. 
% A = [0.8 1 0 0; 0 0.9 1 0; 0 0 0.8 1; 0 0 0 0.9];
% B = [0.1 0.1; 0.1 0.1; 0.1 0.1; 0.1 0.1];
% C = eye(4);

t = 0:100;
% u1 = ones(size(t));
% u2 = ones(size(t))*2;
% u = [u1; u2];

x = [1; 1; 1; 1]; % Initial state
y = zeros(length(t),4);

for k = 1:length(t)
    %x = A*x+B*u(:,k);
    x = A*x;
    y(k,:) = C*x;
end

% figure;
% for i = 1:4
%     subplot(4,1,i);
%     stem(t,y(:,i)); 
%     xlabel('Time (samples)');
%     ylabel(['Output ' num2str(i)]);
% end
figure;
plot(t,y(:,1));
xlabel('Time (samples)');
ylabel(['Output ' num2str(1)]);

x_last = x; % Value of x at the last iteration


%Model 2 measurement from 100 to 200

t1 = 100:200;

x2 = x_last; % Initial state
y2 = zeros(length(t),4);

for k = 1:length(t)
    %x = A*x+B*u(:,k);
    x2 = A2*x2;
    y2(k,:) = C2*x2;
end

% figure;
% for i = 1:4
%     subplot(4,1,i);
%     stem(t,y2(:,i)); 
%     xlabel('Time (samples)');
%     ylabel(['Output ' num2str(i)]);
% end

figure;
plot(t,y2(:,1));
xlabel('Time (samples)');
ylabel(['Output ' num2str(2)]);

x_last2 = x2; % Value of x at the last iteration


%Model 3 initial value of x3 = last value of x2. 

t3 = 200:300;

x3 = x_last2; % Initial state
y3 = zeros(length(t),4);

for k = 1:length(t3)
    %x = A*x+B*u(:,k);
   
    x3 = A3*x3;
    y3(k,:) = C3*x3;
end

% figure;
% for i = 1:4
%     subplot(4,1,i);
%     stem(t,y3(:,i)); 
%     xlabel('Time (samples)');
%     ylabel(['Output ' num2str(i)]);
% end

figure;
plot(t,y3(:,1));
xlabel('Time (samples)');
ylabel(['Output ' num2str(3)]);




