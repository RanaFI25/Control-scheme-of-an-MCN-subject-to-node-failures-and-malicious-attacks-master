%Miroslav Approach%
%A^ = A^(W,H,G) => A = [A B*G;H*C W]%  
clear all 
clc
A = [1 1 0; 0 1 1; 0 0 2]; %Original model 
%Create Different model 
A1 = [0 0 0; 0 1 1; 0 0 2]; % model 1  
A2 = [1 1 0; 0 0 0; 0 0 2]; % model 2
A3 = [1 1 0; 0 1 1; 0 0 2]; % model 3
BO = [0; 0.5; 1];
CO = [1, 0.5 1];
W = [0.228 0.965; -2.872, -1.660];
H = [1;0];
G = [0 1.837];
B = BO*G; 
C = H*CO;
L = W;
Big_A = [A, B; C,  W];
C_Big = eye(5);

r ={};
r{1} = Big_A;
L = {};
L{1} = place(r{1}',eye(5)',[.7 .9 .8 .8 .9])';
for k = 2:6
    Big_A = [A, B; C,  W];
    Big_A(:, k-1) = 0;
    Big_A(k-1,:) = 0;
    r{k} = Big_A;
    L{k} = place(r{k}',eye(5)',[.7 .9 .8 .8 .9])';
end




t = 1:25;
x = [1; 1; 1;1;1]; % Initial state for system 1
xhat = [0; 0;0; 0;0];
XX = zeros(length(t),5);
XXhat = zeros(length(t),5);

measuredResult = {};
estimatedResult = {};
%Measured Output
for i = 1:6
    
    for k = 1:length(t)
        y=C_Big*x;
        yhat=C_Big*xhat;

        x = r{i}*x+B*u;
        xhat=r{i}*xhat+L{i}*(y-yhat);

        XX(k,:) = x;
        XXhat(k,:) = xhat;
    end
    measuredResult{i} = XX;
    estimatedResult{i} = XXhat;
    xhat = x;
end

timeStep = 1:25;
for i = 1:6
    figure;
%     plot(timeStep', measuredResult{i}, 'b-', timeStep', estimatedResult{i}, 'r:.'); 
    plot(timeStep', measuredResult{i} - estimatedResult{i}); 
    xlabel('Time');
    ylabel('Value');
    title(['Plot ' num2str(i)]);
    legend('Residuals')
end

n = 25;
y_random = {};
t = 1:n:300;

for t= 1:12
    random_number = randi([1, 6]);
    y_random{t} = measuredResult{random_number};
end
residuals = [];
r_all = {};
Model_Name = "";
output_Model = [];
for k = 1:12
    for j = 1:6
        diff = estimatedResult{j} - y_random{k}; 
        if all(diff <= 1 & diff >= -1)
            output_Model = [output_Model, j]
            %residuals = [residuals, residual];
            Model_Name = Model_Name+"Output "+k+" correspond to model :"+j+newline
        end
    end
   
end


% treshHoldUpper = 1; 
% treshHoldLower = -1; 
% Model_Name = "";
% output_Model = [];
% for k = 1:12
%      for j = 1:6
%         if all(r_all{k}(:,j) <= treshHoldUpper & r_all{k}(:,j) >= treshHoldLower)
%             disp("in")
%             output_Model = [output_Model, j]
%             Model_Name = Model_Name+"Output "+k+" correspond to model :"+j+newline
%         end
%      end
% end

%Check Controllability: 

% %Method 1 same as before%
% 
% Cm = ctrb(r{1},zeros);
% Om = obsv(r{5},zeros(5,5));
% if rank(Cm) == size(r{1},1)
%     'Controllable'
% else
%     'Not Controllable'
% end
% 
% if rank(Om) == size(r{1},1)
%     'Observable'
% else
%     'Not Observable'
% end



%Hautus Test Method [A-lambda*I,B]%

% System matrices
A = [1 2; 3 4];
B = [1; 1];
C = [1 0];
D = 0;

% Hautus test for controllability
% n = size(A,1);
% m = size(B,2);
% s = sym('s');
% M = sym(zeros(n+m,n+m));
% M(1:n,1:n) = s*eye(n) - A;
% M(1:n,n+1:n+m) = B;
% if rank(M) == n+m
%     disp('System is controllable');
% else
%     disp('System is not controllable');
% end


% Original system
Aqq = [0 1 0; 0 0 1; 0 0 0];
Bqq = [0; 0; 1];
C = [1 0 0];
D = 0;

% Check controllability
if rank(ctrb(A,B)) == size(A,1)
    disp('System is controllable')
else
    disp('System is not controllable')
end

% State transformation
T = [B, A*B, A^2*B];
T_inv = inv(T);
A_tilde = T \ (A * T);
B_tilde = T \ B;

% Check controllability of transformed system
if rank(ctrb(A_tilde,B_tilde)) == size(A_tilde,1)
    disp('Transformed system is controllable')
else
    disp('Transformed system is not controllable')
end

% Design controller for transformed system
p = [-1 -2 -3]; % desired poles
K_tilde = place(A_tilde, B_tilde, p); % pole placement
K = K_tilde * T_inv; % controller gain in original coordinates

% Run F


%%

% for k = 1:length(t)
%     x = A*x + B*u(:,k);
%     y(k,:) = C*x;
% end
% 
% for k = 1:length(t)
%     x = A1*x + B*u(:,k);
%     y1(k,:) = C*x;
% end
% 
% for k = 1:length(t)
%     x = A2*x + B*u(:,k);
%     y2(k,:) = C*x;
% end
% 
% for k = 1:length(t)
%     x = A3*x + B*u(:,k);
%     y3(k,:) = C*x;
% end
% 
% y_All = [y(:,1);y1(:,1);y2(:,1);y3(:,1)];
% 
% x_hat = [1.01; 1.01; 1.01]; % Initial estimated state
% y_hat = zeros(length(t),2);
% for k = 1:length(t)
%     x_hat = A*x_hat+ B*u(:,k) ;
%     y_hat(k,:) = C*x_hat;
%     x_hat = x_hat + W * (y(k,:) - y_hat(k,:))'; % Luenberger observer % Luenberger observer
% end
% 
%Row should be zero as well%

