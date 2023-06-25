clear all
clc
A1 = [0.8 1 0 0; 0 0.9 1 0; 0 0 0.8 1; 0 0 0 0.9];
B1 = [0.1 0.1; 0.1 0.1; 0.1 0.1; 0.1 0.1];
C1 = eye(4);

A2 = [0.7 1 0 0; 0 0.8 1 0; 0 0 0.7 1; 0 0 0 0.8];
B2 = [0.2 0.2; 0.2 0.2; 0.2 0.2; 0.2 0.2];
C2 = eye(4);

A3 = [0.6 1 0 0; 0 0.7 1 0; 0 0 0.6 1; 0 0 0 0.7];
B3 = [0.3 0.3; 0.3 0.3; 0.3 0.3; 0.3 0.3];
C3 = eye(4);
L = place(A1',C1',[.5 .7 .5 .7])';
L1 = place(A2',C2',[.5 .7 .5 .7])';
L2 = place(A3',C3',[.5 .7 .5 .7])';
t = 1:100;
u1 = ones(size(t));
u2 = ones(size(t))*2;
u = [u1; u2];

x1 = [1; 1; 1; 1]; % Initial state for system 1
x_hat = [1.01;1.01;1.01;1.01]; % Initial estimated state
y1 = zeros(length(t),4);
y_hat = zeros(length(t),4);
for k = 1:length(t)
    x1 = A1*x1 + B1*u(:,k);
    y1(k,:) = C1*x1;
    x_hat = A1*x_hat+ B1*u(:,k) ;
    y_hat(k,:) = C1*x_hat;
    x_hat = x_hat + L * (y1(k,:) - y_hat(k,:))'; % Luenberger observer % Luenberger observer
end

x2 = x1; % Initial state for system 2 is the last iteration value of system 1
y2 = zeros(length(t),4);
y_hat2 = zeros(length(t),4);
x_hat2 = x1+0.01;
for k = 1:length(t)
    x2 = A2*x2 + B2*u(:,k);
    y2(k,:) = C2*x2;
    x_hat2 = A2*x_hat2 + B2*u(:,k) ;
    y_hat2(k,:) = C2*x_hat2;
    x_hat2 = x_hat2 + L1 * (y2(k,:) - y_hat2(k,:))'; % Luenberger observer
end

x3 = x2; % Initial state for system 3 is the last iteration value of system 2
y3 = zeros(length(t),4);
x_hat3 = x2+0.01;
y_hat3 = zeros(length(t),4);
for k = 1:length(t)
    x3 = A3*x3 + B3*u(:,k);
    y3(k,:) = C3*x3;
    x_hat3 = A3*x_hat3 + B3*u(:,k);
    y_hat3(k,:) = C3*x_hat3;
    x_hat3 = x_hat3 + L2 * (y3(k,:) - y_hat3(k,:))'; % Luenberger observer '?
end
 
% figure;
% for i = 1:4
%     subplot(4,1,i);
%     stem(0:100,y1(1:101,i), 'r'); hold on;
%     stem(0:100,y_hat(1:101,i),'--b');
%     stem(100:200,y2(1:101,i),'g');
%     stem(100:200,y_hat2(1:101,i),'--r');
%     stem(200:300,y3(1:101,i),'b');
%     stem(200:300,y_hat3(1:101,i),'--r');
%     xlabel('Time (samples)');
%     ylabel(['Output ' num2str(i)]);
%     legend({'System 1 actual','System 1 estimated','System 2 actual','System 2 estimated','System 3 actual','System 3 estimated'});
% end

% Have estimation for each model. y_hat, y_hat2, y_hat3. 
%combine three different output as one output 
 
y = [y1; y2; y3];% Combined output
y_split = mat2cell(y, 25*ones(1,12), size(y,2));% split output in 12 part so that you have 25 time step every time of output. 



%Than run estimators on these splitted part of the output. 
% get first 25 timestep and measure the difference between all estimators. 
% e1 = y_split{2} - y_hat(26:50,:) % first 25 time step with first estimators 25 time step. 
% e2 = y_split{2} - y_hat2(26:50,:)
% e3 = y_split{2} - y_hat3(26:50,:)

% y_split has 12 part determine which model each part correspond to. 

% result = cell(12,4)
% for k = 1:12
%     for i = 0:3
%         e1 = y_split{k} - y_hat(i*25+1:i*25+25,:); % first 25 time step with first estimators 25 time step. 
%         e2 = y_split{k} - y_hat2(i*25+1:i*25+25,:);
%         e3 = y_split{k} - y_hat3(i*25+1:i*25+25,:);
%         result{k,i+1} = [e1(:,1),e2(:,1),e3(:,1)]
%     end
% end

result = cell(12,1);
var = 0; 
treshHoldUpper = 1; 
treshHoldLower = -1; 
output_Model = [];
Model_Name = "";
for k = 1:12
    
    e1 = y_split{k} - y_hat(var*25+1:var*25+25,:); % first 25 time step with first estimators 25 time step. 
    e2 = y_split{k} - y_hat2(var*25+1:var*25+25,:);
    e3 = y_split{k} - y_hat3(var*25+1:var*25+25,:);
    result{k} = [e1(:,1),e2(:,1),e3(:,1)];
    if (var+1 >= 4); var=0; else; var=var+1; end
    for j = 1:3
        if((result{k}(:,j)<=treshHoldUpper) & (result{k}(:,j) >= -treshHoldLower))
            output_Model = [output_Model, j];
            Model_Name = Model_Name+"Output "+k+" correspond to model :"+j+newline;  
        end
    end  
end



for i = 1:12
    figure;
    plot(1:25, result{i})
    xlabel('Time (samples)');
    ylabel('Error');
    legend({'Estimated Model 1','Estimated Model 2','Estimated Model 3'});
    %Analyze the error. 
end


