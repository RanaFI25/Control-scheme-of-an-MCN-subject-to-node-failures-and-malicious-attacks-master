% Different models (Different A's, same B, and Identity C)%

clear all
clc
close all

vars = load('modellllllllll.mat');



all_models ={};
C_matrix = {};
B_matrix = {};


all_models = assing_val(vars,'A');



B_matrix = assing_val(vars,'B');
C_matrix = assing_val(vars,'C');


all_models;% Now all models contains original model, model 1 to model 6 where in each model one sensor has a fault.
len_models = length(all_models);
% Check controllability and observability.

[controllable, observeable] = contAndObsv(all_models,B_matrix,C_matrix);

controllable
observeable


observerGain = {};
controllerGain = {};

% compute state feedback gain using pole placement

closer_to_One_poles =  [.9,.8,.7];
closer_to_zero_poles = [.2,.35,.25];
between_zero_one_poles = [.5,.55,.45];

desired_poles = closer_to_One_poles;



for k = 1:len_models

    if(controllable && observeable)
%         sys = ss(all_models{k},B_matrix{k},C_matrix{k},0);
%         figure;
%         step(sys)
        observerGain{k} = place(all_models{k}', C_matrix{k}', desired_poles)';
        controllerGain{k} = place(all_models{k}, B_matrix{k}, desired_poles);
    else
        disp('Not controllable and observable');
        break;
    end
end



%Observer gain

observerGain; %observerGain 
controllerGain; % controllerGain 

stable = stability_check(all_models, B_matrix,controllerGain)



%Pure state feedback%
t = 1:50;

[estimatedResult,estimatedResult_y,measuredResult, measuredResult_y] = create_LObservers(all_models, observerGain,C_matrix,t,B_matrix, controllerGain);


% figure;
% set(gca,'FontSize',20)
% for i = 1:size(measuredResult, 1)
%     for j = 1:size(measuredResult, 2)
%         subplot(size(measuredResult, 1), size(measuredResult, 2), (i-1)*size(measuredResult, 2) + j);
%         hold on;
%         txt = ['Model: ', num2str(i), ' C: ', num2str(j)];
%         semilogy(estimatedResult_y{i, j}  + eps, 'DisplayName', txt, Color='r');
%         xlabel('Time step')
%         ylabel('Estimated output')
%         
%         legend();
%     end
% end





% x_new = ones(size(all_models{1},1),1);
% 
% XX_new = zeros(length(t),size(all_models{1},1));
% YY_new = zeros(length(t)*length(model),1);
% pos = 1;
% for i = 1:length(model)
%     currentK = controllerGain{ctlr(i)};
%     for j = 1:length(t)
%         %u_new = -currentK*x_new; %control law where u = -K*x
%         y_new=C_matrix{model(i)}*x_new;
%       
% 
%         x_new = (all_models{model(i)}+B_matrix{model(i)}*currentK)*x_new; %Ax +BKx
%         XX_new(pos,:) = x_new;
%         YY_new(pos,:) = y_new;
%         pos = pos+1;
%     end
% end
% hold on;
% plot(estimatedResult_y{1,1}-YY_new(1:50,:));
% 
% legend('First')
% hold off;
% 
% figure;
% hold on;
% plot(estimatedResult_y{2,1}-YY_new(1:50,:));
% 
% legend('Second')
% hold off;
% 
% figure;
% hold on;
% plot(estimatedResult_y{3,1}-YY_new(1:50,:));
% 
% legend('Third')
% hold off;

%plot(estimatedResult_y{1,1}-YY_new(1:50,:))

% treshhold



% start = 1;
% stop = 20;
% start_func = 1;
% stop_func =20;
% pos_remove = 2;
% 
% start_test_var= 1;
% 
% test_var= 20;
% val =20;
% e1=measuredResult_y{3,3}(start_test_var:test_var,:)-estimatedResult_y{1,1}(1:val,:);
% e2=measuredResult_y{3,3}(start_test_var:test_var,:)-estimatedResult_y{2,1}(1:val,:);
% e3=measuredResult_y{3,3}(start_test_var:test_var,:)-estimatedResult_y{3,3}(1:val,:);
% 
% mean_e1= abs(mean(e1));
% var_e1 = var(e1);
% mean_e2 = abs(mean(e2));
% mean_e3 = abs(mean(e3));
% var_e2 = var(e2);
% mean_e1<mean_e2
% var_e1<var_e2






steps = 10;
model = [1, 2, 3]%, 3, 2]%, 1];
ctlr = [1, 1, 2]%, 3, 2]%, 1];

% ctlr = %[1,2][3,3],2,1]%,3,3,2,1,2,2,2];
 %model = %[1,1][2,3],2,1]%,3,1,2,3,1,2,3];
% ctlr = [3,1];
% model = [1,3];
%ctlr = [1,2]%,3,3,2,1,2,2,2];
%model = [1,1]%,3,1,2,3,1,2,3];
begin_var = 1;
last_end = steps;
curr_residual = {};
for i = 1: length(model)
   %subplot(1,2,i)
    colors = ['b','r','c','g','y','k'];
    time = linspace(0,steps);
    for j = 0:(length(t)/steps)-1
        [modelNr,controllerNr,currentResidual]  = closest_Model2(measuredResult_y,estimatedResult_y,model(i),ctlr(i), (j*steps)+1,steps);
%         hold on
%         txt = ['Model: ',num2str(model(i)), ' C: ',num2str(ctlr(i)) ];
%         semilogy(begin_var:last_end,measuredResult_y{model(i),ctlr(i)}((j*steps)+1:((j*steps)+1+steps)-1,:)-estimatedResult_y{model(i),ctlr(i)}((j*steps)+1:((j*steps)+1+steps)-1,:),'DisplayName',txt,'Color',colors(i),'LineStyle','-','LineWidth',1.5,'MarkerSize',10)
       % eigenvalues = eig(all_models{model(i)} - B_matrix{model(i)}*controllerGain{i}); % K is the control law
        %disp(eigenvalues);
       
        title(['Model ', num2str(model(i)), ' Controller ', num2str(ctlr(i))]);
        xlabel('Time step');
        ylabel('Residuals')
        %disp(['Model: ', num2str(modelNr), ' Controller: ', num2str(controllerNr)]);
          begin_var = last_end+1;
        last_end = last_end+steps;
        if(modelNr ~= controllerNr)
            disp('Current Model is running with wrong controller');
            
            ctlr(i) = modelNr;
         txt = ['Eigenvalues ', 'model: ', num2str(model(i)), ' c: ', num2str(ctlr(i)), ' at time step: ', num2str(begin_var)];
        
        eigenvalues1 = eig(all_models{model(i)} - B_matrix{model(i)}*controllerGain{ctlr(i)});
        eigenvalues = eig(all_models{model(i)} - observerGain{model(i)}*C_matrix{model(i)});
        
        scatter(real(eigenvalues), imag(eigenvalues), 'o', 'DisplayName', ['Observer Eigenvalues',' at Time step: ', num2str(begin_var)], 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black');
        hold on;
        scatter(real(eigenvalues1), imag(eigenvalues1), 'o', 'DisplayName', ['Controller Eigenvalues', ' at Time step: ', num2str(begin_var)], 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black', 'SizeData', 100);
        hold on;
        
        viscircles([0, 0], 1, 'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'black');
        
        xlabel('Real part');
        ylabel('Imaginary part');
        title('Eigenvalues of the closed-loop system');
        grid on;
        legend('Location', 'best');
        hold off;



        end
       
       
        curr_residual{i,j+1} = currentResidual;
    end
end
0,
ebegin = 1; 
last = steps;
for i = 1: length(model)
   
  for j = 0:(length(t)/steps)-1
      hold on;
      %plot((j*steps)+1:(steps)+(j*steps),curr_residual{i,j+1});
      plot(begin:last,curr_residual{i,j+1})
      begin = last+1;
      last = last+steps;
    %legend('Measured model '  + num2str(i));  
  end
end

residuals_last = {};
bol = false;
for step = 1:100

    %[current_model,residuals_last] = closest_Model(residuals_last,YY_new(start:stop,:),estimatedResult_y,start_func,stop_func);
    %[current_model,residuals_last]  = closest_Model(residuals_last,YY_new(start_func:stop_func,:),estimatedResult_y,start,stop);

  
    modelnr  = closest_Model2(measuredResult_y(start_func:stop_func,:),estimatedResult_y,start,stop,controller)
    %     bol = within_Treshhold(YY_new(start:stop,:), estimatedResult_y,start_func,stop_func,upper_bound1,lower_bound1,curr_model,curr_controller);
    %      if(~bol)
    %         disp(["In",num2str(start)]);
    %         %not correct model with correct controller
    %         curr_model = model_check(YY_new(start:stop,:),1,estimatedResult_y,start_func,stop_func)
    %         if(curr_model==model(1))
    %             ctlr(1) = curr_model;
    %             curr_controller = ctlr(1);
    %         end
    %         %x = XX(start-1,:);
    %         %XX = generate_output(all_models,model,ctlr,B_matrix,C_matrix,controllerGain,x',false);
    %         %YY = generate_output(all_models,model,ctlr,B_matrix,C_matrix,controllerGain,x',true);
    %
    %         start_func = 0;
    %         stop_func = 0;
    %         upper_bound1 =  max(max(estimatedResult_y{model(1),ctlr(1)}-YY_new(51:100,:)))
    %         lower_bound1 = min(min(estimatedResult_y{model(1),ctlr(1)}-YY_new(51:100,:)))
    %     end
    start = start+10;
    stop = stop+10;
    start_func = start_func+10;
    stop_func = stop_func+10;

    if(stop >50 )%|| start_func==51 || start_func ==101)
        start = 1;
        stop= 20;
        model = model(pos_remove:end);
        ctlr = ctlr(pos_remove:end);
        

  
    end
    if  modelnr ~= controller
        disp("Active model is currently running with wrong controller.")
%         ctlr(1) = modelnr;
%         x_temp = XX_new(start_func,:);
%         %x_temp = ones(size(all_models{1},1),1);
%         YY_new = generate_output(all_models,model,ctlr,B_matrix,C_matrix,controllerGain,x_temp',true);
%         start = 1;
%         stop = 20;
%         start_func = 1;
%         stop_func =20;
    end

end

% % Define the intervals
% intervals = [-10.6519 9.95545; -7.87725 12.099; -17.6081 10.7586;
% -17.9106 10.6162; -15.5885 20.0188];
%
% % Calculate the absolute value of the lower and upper bounds of each
% interval
% abs_lower_bounds = abs(intervals(:,1)); 0.3371 + 0.4344i   0.3371 - 0.4344i  -0.6142 - 0.5919i  -0.6142 + 0.5919i  -0.3497 - 0.0292i  -0.3497 + 0.0292i
% abs_upper_bounds = abs(intervals(:,2));
%
% % Find the interval with the smallest absolute value
% [~, closest_interval_index] = min(min(abs_lower_bounds, abs_upper_bounds));
%
% % Print the closest interval
% closest_interval = intervals(closest_interval_index,:);


function bool = within_Treshhold(chunk,estimated,start,stop,upperbound,lowerbound,model,k)

%for k = 1:length(estimated)
residuals = estimated{model,k}(start:stop,:)-chunk;
bolupper = all(residuals <= upperbound);
bolLower = all(residuals >= lowerbound);

if(all(bolupper) && all(bolLower))
    bool = true;

else
    bool = false;

end
%end
end

function correct_model = model_check(chunk,controller,estimated,start,stop)
upper_bound = [];
lower_bound = [];

for i = 1:3
    upper_bound = [upper_bound,max(max(estimated{i,controller}(start:stop,:)-chunk))];
    lower_bound = [lower_bound, min(min(estimated{i,controller}(start:stop,:)-chunk))];
end


[~,correct_model]=min(upper_bound);
correct_model = correct_model+1;
end

function desired_poles = poles_p(num_poles)


desired_poles = zeros(1, num_poles);


for i = 1:2:num_poles
    pole = rand()*exp(1i*rand()*2*pi);
    while abs(pole) >= 1 % Check if pole is outside the unit circle
        pole = rand()*exp(1i*rand()*2*pi);
    end
    desired_poles(i:i+1) = [pole, conj(pole)];
end
end

function app_model =  model_checker(randomKResult,index,estimatedResult)
list_models = [];

for i = 1:length(estimatedResult)
    for k = 1:5
        %val = abs(randomKResult{index}(:,1)/max(abs(randomKResult{index}(:,1))) - abs(estimatedResult{i,k}(:,1)/max(abs(estimatedResult{i,k}(:,1)))));
        val = randomKResult{index}(:,1) - estimatedResult{i,k}(:,1);
        list_models =[list_models, val];

    end



end

app_model = list_models;
end

function scaled_array = scaled_version(list)
scaling_factor = 1;
scaled_array = scaling_factor *list;
while ~(max(abs(scaled_array)) > 10)

    % Multiply the scaling factor by 10
    scaling_factor = scaling_factor * 10;

    % Multiply the input array by the new scaling factor
    scaled_array = scaling_factor * list;

end
end




function average = Model_closer_to_zero(list, estimatedResult)



mean_abs_A = mean(abs(list)); % Calculate the mean absolute value of elements in A
%mean_abs_B = mean(abs(B)); % Calculate the mean absolute value of elements in B
disp(mean_abs_A);
mean_abs = [];

for i = 1:5
    scaled_estimated = scaled_version(estimatedResult{i,1}(:,1));
    %mean_abs = [mean_abs, mean(abs(list-scaled_estimated))];
    mean_abs = [mean_abs, mean(abs(scaled_estimated))];
end
%[~,index] = min(abs(mean_abs_A-mean_abs));
[~, index] = min(abs(mean_abs - mean_abs_A));
disp(abs(mean_abs - mean_abs_A));
average = index;
end
%me lene
function output = generate_output(all_models,model,ctlr,B_matrix,C_matrix,controllerGain,x,y_data)
t = 1:50;
XX_new = zeros(length(t),size(all_models{1},1));
YY_new = zeros(length(t)*length(model),1);
pos = 1;
for i = 1:length(model)
    currentK = controllerGain{ctlr(i)};
    % bring all controllers controller i.  %k1,k2,k1,k3,k1
    for j = 1:length(t)
        y_new=C_matrix{model(i)}*x;
        x = (all_models{model(i)}+B_matrix{model(i)}*currentK)*x; %Ax+B*u  A0
        XX_new(pos,:) = x;
        YY_new(pos,:)=y_new;
        pos = pos+1;
    end
end


if(y_data)
    output = YY_new;
else
    output = XX_new;
end


end

function assign_values = assing_val(vars,stringToMatch)
fieldNames = fieldnames(vars);
matchedFields = fieldNames(startsWith(fieldNames, stringToMatch));


for i = 1:length(matchedFields)
    varName = matchedFields{i};
    assign_values{i} = vars.(varName);

end

end

function [model,res_cell] = closest_Model(residuals_lastime,currentdata, estimateddata,start,stop)
%create residuals
residuals = {};
total = [0,0,0];
for i = 1:3
    if(length(residuals_lastime)>=1)
        residuals{i} = currentdata-estimateddata{i,1}(start:stop,:);
        residuals{i} = abs(residuals{i})-abs(residuals_lastime{i});

    else
        residuals{i} = abs(currentdata)-abs(estimateddata{i,1}(start:stop,:));

    end

end
residuals{:}

for k = 1:5
    temp = [residuals{1}(k),residuals{2}(k),residuals{3}(k)];
    [~,min_idx] = min(temp);
    total(min_idx) = total(min_idx)+1;

end

[~,model] = max(total);
res_cell = residuals;

end

function [curr_Model,curr_ctrl,curr_residual] = closest_Model2(currentdata, estimateddata,model,controller,crr_step,steps)

smallest_mean = 99999999;
curr_Model = -1;
curr_ctrl = -1;
%begin_var = 1;
%last_end = 10;
for i = 1: size(estimateddata,1)
     %Color = ['b','r','c'];
    for j = 1:size(estimateddata,2)
        residual=currentdata{model,controller}(crr_step:(crr_step+steps)-1,:)-estimateddata{i,j}(crr_step:(crr_step+steps)-1,:);
        mean_e1= abs(mean(residual));
%         hold on
%         txt = ['Model: ',num2str(i), ' C: ',num2str(j) ];
%         semilogy(crr_step:(crr_step+steps)-1,residual,'DisplayName',txt,LineWidth=2)
%         legend();
        %var_e1 = var(residual);
        if(mean_e1 < smallest_mean)
            curr_residual = residual;
            smallest_mean = mean_e1;
            curr_Model = i;
            curr_ctrl = j;
           
        end
    end
end

% e1=currentdata{model,controller}(start:stop,:)-estimateddata{1,controller}(start:stop,:);
% e2=currentdata-estimateddata{2,controller}(start:stop,:);
% e3=currentdata-estimateddata{3,controller}(start:stop,:);
% mean_e1= abs(mean(e1));
% var_e1 = var(e1);
% mean_e2 = abs(mean(e2));
% var_e2 = var(e2);
% mean_e3 = abs(mean(e3));
% var_e3 = var(e3);
% 
% if mean_e1<mean_e2 && mean_e1<mean_e3 && var_e1<var_e2 && var_e1<var_e3 && mean_e1~=mean_e2
%     model = 1;
% elseif ~(mean_e1<mean_e2) && mean_e2<mean_e3 || ~(var_e1<var_e2) && var_e2<var_e3
%     model = 2;
% elseif ~(mean_e2<mean_e3) && ~(mean_e1<mean_e3) || ~(var_e1<var_e3) && ~(var_e2<var_e3)
%     model = 3;
% else
%     model = -1;
% end


end
function stable = stability_check(all_models, B_matrix,controllerGain)
len_models = length(all_models);
% Compute the closed-loop system matrix A
    for i = 1:len_models
         Am = all_models{i} - B_matrix{i}*controllerGain{i};
%         syscl = ss(Am,B_matrix{i},C_matrix{i},0);
%         figure;
%         step(syscl);
        %Pcl = pole(syscl);
        % Check if the closed-loop system is stable for discrete time
        if all(abs(eig(Am)) < 1)
            %disp('The system is stable');
            stable = 1;
        else
            %disp('The system is unstable');
            stable = 0;
            break;
        end
    end

end

function [estimatedResult,estimatedResult_y,measuredResult, measuredResult_y] = create_LObservers(all_models, observerGain,C_matrix,t,B_matrix,controllerGain)
    len_models = length(all_models);
    x = ones(size(all_models{1},1),1); %initial x
    XX = zeros(length(t),size(all_models{1},1));
    yy = zeros(length(t),1);
    measuredResult = {};
    measuredResult_y = {}; 
    
    %observer based state feedback
    
    xhat = zeros(size(all_models{1},1),1);
    XXhat = zeros(length(t),size(all_models{1},1));
    yyhat = zeros(length(t),1);
    estimatedResult = {};
    estimatedResult_y = {};
    for i = 1:len_models % Current model
        for k = 1:len_models % Current controller
            for j = 1:length(t)
    
                u = -controllerGain{i}*xhat; %control law where u = -K*x^
                %u = -controllerGain{i}*xhat;%u(t) = −Lˆx(t)
    
                y=C_matrix{i}*x;
                yhat=C_matrix{i}*xhat;
    
                x = all_models{i}*x+ B_matrix{i}*u;
    
                %xhat = all_models{i}*xhat + observerGain{i}*(y-C_matrix{i}*xhat); %Axhat + L(y-Cxhat)
                
                %xhat = (all_models{i}-observerGain{k}*C_matrix{i})*xhat +observerGain{k}*y; %(A-LC)xhat + Ly.
                xhat=all_models{i}*xhat+B_matrix{i}*u+observerGain{i}*(y-yhat); % x(t) = Aˆx(t) + Bu(t) + L(y(t) −C ˆx(t))
                %xhat= (all_models{i}-observerGain{i}*C_matrix{i}-B_matrix{i}*controllerGain{k})*xhat + observerGain{i}*C_matrix{i}*x;    %(A −KC)ˆx(t) + Bu(t) + Ky(t)%
                
                XX(j,:) = x;
                XXhat(j,:) = xhat; %observer based
                yy(j,:) = y;
                yyhat(j,:) = yhat;
            end
            measuredResult{i,k} = XX; %Model, controller e.g Model1controller1, model1controller2 ..
            estimatedResult{i,k} = XXhat;
            measuredResult_y{i,k} = yy;
            estimatedResult_y{i,k} = yyhat;
            
           %xhat = x;
        end
    end
end
%mode detection
%multimode diagnosis on switched affine system with noisy measurement
