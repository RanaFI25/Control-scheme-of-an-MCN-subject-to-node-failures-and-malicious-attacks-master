%% for j = 1:sched_length(2)
%            sz1 = size(sched_funcb{j}); % size of current list in cell like [1 2];[2 3];[3 4] gives sz1 3 3 and we need to take sz(1)
%            %route_check{j} = sched_funcb{j}(1:min(d,sz1(1)),:);
%            %sz = size(route_check{j});
%            for k = 1:sz(1) 
%                route_final{current} = route_check{j}(k,:);
%                current = current +1 ; 
%            end
%         end
%         current = 1; 
%         if length(route_final) == maxDR
%             % means we have now enough elements to go from source to
%             % destination now question is they are in order. 
%            
%            delay = d; 
%               curr = 1;
%             while ~isempty(route_final)
%                 if all(pathsCopy{1} == route_final{curr})
%                     pathsCopy(1) = [];
%                     route_final(curr) = [];
%                      curr = 1;
%                      iteration = true; 
%                 else
%                     if iteration || delay == d
%                         delay = delay +1;
%                         iteration = false;
%                     end
%                     curr = curr+1;
%                         
%                 end
%             end
%                            
%             
%                 
%            
%             break;
%         end



% n = [1 2 1];
% d = [2 2 1];
% g = tf(n,d);
% bode(g)
% nyquist(g);
% 
% %% Define the transfer function in the z-domain
% num = [0 3/5 2/5];
% den = [1 0 0];
% sys_z = tf(num, den, 1, 'Variable', 'z');
% 
% 
% % Compute the frequency response using freqz
% [h, w] = freqz(num, den);
% 
% % Plot the magnitude and phase response
% subplot(2, 1, 1);
% plot(w, abs(h));
% title('Magnitude Response');
% xlabel('Frequency (radians/sample)');
% ylabel('Magnitude (dB)');
% grid on;
% 
% subplot(2, 1, 2);
% plot(w, angle(h));
% title('Phase Response');
% xlabel('Frequency (radians/sample)');
% ylabel('Phase (radians)');
% grid on;
% 
% 
% % Compute the step response using stepz
% n = 0:50; % Number of samples
% h = stepz(num, den, n);
% 
% % Plot the step response
% figure;
% plot(n,h);
% title('Step Response');
% xlabel('Samples');
% ylabel('Amplitude');
% grid on;

%Magnitude response: how the amplitude of the system's output varies with input frequnecy. 
%Phase Response: a phase shift of 180 degrees or more at a certain
%frequency indicates that the system is unstable at that frequency. 
%Cutoff frequency: The cutoff frequency of a system is the frequency at which the magnitude response drops by 3 dB (or 0.707 times) from its maximum value. In the z-domain, the cutoff frequency is related to the pole locations of the transfer function. If the poles are close to the unit circle, the cutoff frequency is high, 
% indicating that the system responds quickly to changes in the input. If the poles are far from the unit circle, the cutoff frequency is low, indicating that the system responds slowly to changes in the input.


% upper_bound1 =  max(max(estimatedResult_y{1,1}-YY_new(1:50,:)));
% lower_bound1 = min(min(estimatedResult_y{1,1}-YY_new(1:50,:)));
% upper_bound2 =  max(max(estimatedResult_y{2,1}-YY_new(51:100,:)));
% lower_bound2 = min(min(estimatedResult_y{2,1}-YY_new(51:100,:)));
% upper_bound3 =  max(max(estimatedResult_y{3,1}-YY_new(101:150,:)));
% lower_bound3 = min(min(estimatedResult_y{3,1}-YY_new(101:150,:)));
%curr_model = 1;
%residuals_last = {};
%last_start = 1;
%last_stop = 20;
%curr_controller = 1;



delaychecker;

% x = 1:25;
% var2 = 5; 
% var3 = 5; 
% 
% y1 = measuredResult{1}(:,var2);
% y2 = measuredResult{2}(:,var2);
% y3 = measuredResult{3}(:,var2);
% y4 = measuredResult{4}(:,var2);
% y5 = measuredResult{5}(:,var2);
% 
% x = 1:25;
% e1 = estimatedResult{1}(:,var3);
% e2 = estimatedResult{2}(:,var3);
% e3 = estimatedResult{3}(:,var3);
% e4 = estimatedResult{4}(:,var3);
% e5 = estimatedResult{5}(:,var3);
% 
% % create the figure
% figure
% 
% % plot the first subplot in the first row
% subplot(2, 3, 1)
% plot(x, y1) 
% hold on
% plot(x, e1)
% title('Plot 1')
% 
% % plot the second subplot in the first row
% subplot(2, 3, 2)
% plot(x, y2) 
% hold on
% plot(x, e2)
% title('Plot 2')
% 
% % plot the third subplot in the first row
% subplot(2, 3, 3)
% plot(x, y3) 
% hold on
% plot(x, e3)
% title('Plot 3')
% 
% % plot the fourth subplot in the second row
% subplot(2, 3, 4)
% plot(x, y4) 
% hold on
% plot(x, e4)
% title('Plot 4')
% 
% % plot the fifth subplot in the second row
% subplot(2, 3, 5)
% plot(x, y5) 
% hold on
% plot(x, e5)
% title('Plot 5')
% 
% figure
% 
% % plot the first subplot in the first row
% subplot(2, 3, 1)
% plot(x, y1-e1) 
% title('Plot 1')
% 
% % plot the second subplot in the first row
% subplot(2, 3, 2)
% plot(x, y2-e2) 
% 
% title('Plot 2')
% 
% % plot the third subplot in the first row
% subplot(2, 3, 3)
% plot(x, y3-e3) 
% 
% title('Plot 3')
% 
% % plot the fourth subplot in the second row
% subplot(2, 3, 4)
% plot(x, y4-e4) 
% 
% title('Plot 4')
% 
% % plot the fifth subplot in the second row
% subplot(2, 3, 5)
% plot(x, y5-e5) 
% 
% title('Plot 5')
% 
% %Create a random output with different K's




%If model already exist as one of the 6 estimated%

% for i = 1:6
%     for k = 1:6
%         diff = randomKResult{1} - estimatedResult{i,k}; %A1K1 - E(A1K1 k2 k3 k4 k5..)
%         if diff <= 1 & diff >= -1
%              disp([num2str(i), ' ' ,num2str(k), ' Exisiting in estimated'])
%              break
%         end
%     end
% end
% 
% 
% t = 1:25;
% x = [1; 1; 1;1;1]; %initial x
% 
% XX = zeros(length(t),5);
% randomKResult2 = {}; 
% controllernr = 1; 
% modelnr = 1; 
%  
%   while true
%     currentK = controllerGain{controllernr}; % bring controller i. 
%        for j = 1:length(t)
%             u = -currentK*x; %control law where u = -K*x
%             y=C*x;           
%             x = all_models{modelnr}*x + B*u; %Ax+B*u  
%             XX(j,:) = x;
%        end
%        x = [1; 1; 1;1;1]; %initial x
%        randomKResult2{1} = XX;
%        diff = randomKResult{1} - randomKResult2{1};
%        if diff <= 1 & diff >= -1
%             disp([num2str(modelnr), ' ' ,num2str(controllernr), ' OK'])
%             break
%        end
%        controllernr = controllernr + 1;
%        if controllernr > 6
%            if modelnr+1 ~= 7 
%             modelnr = modelnr+1;
%            end
%            controllernr = 1;
%        end
%   end
%Made state estimation for each model with different controllers i.e state estimation for A1 with k1 to k6%
%During state estimation x and xhat reset at each iteration because having last x value from previous model as initial value leads to high residual%
%This was check by generating random model with random controller several time and check if difference between them is within treshold i.e -1 to 1%
%Don't reset x and xhat%  klar
%Generate output start with A0K0 then A1K0%




% Controllabiltiy  

% if rank(ctrb(A,B)) == size(A,1)
%     disp(['Model ' , num2str(1) ,' is controllable']);
% else
%     disp(['Model ' , num2str(1) ,' not controllable']);
% end




% Compute the frequency response using freqz

% [h, w] = freqz(num_m, den_m);
% phi = 180*unwrap(angle(h))/pi;
% % Plot the magnitude and phase response
% subplot(2, 1, 1);
% 
% plot(w, abs(h));
% title('Magnitude Response');
% xlabel('Frequency (radians/sample)');
% ylabel('Magnitude (dB)');
% grid on;
% 
% subplot(2, 1, 2);
% plot(w, phi);
% title('Phase Response');
% xlabel('Frequency (radians/sample)');
% ylabel('Phase (radians)');
% grid on;
% 
% 
% % Compute the step response using stepz
% n = 0:50; % Number of samples
% h = stepz(num_m, den_m, n);
% 
% % Plot the step response
% figure;
% plot(n,h);
% title('Step Response');
% xlabel('Samples');
% ylabel('Amplitude');
% grid on;




% maxDR = max_num_hops;
% pathsCopy = paths{1};
% sched_length = length(sched_func);
% delay = 1;
% iteration = false;
% route_check = {};
% route_final ={} ;
% current = 1;
% sep_delay = {};
% delay_all_paths = {};
% 
% if sched_length == maxDR % if the order of scheduling function is same as flow of path then in this case each hop can be done in different periods-
%     curr = 1;
%     while ~isempty(sched_func)
%         if all(pathsCopy{1} == sched_func{curr})
%             pathsCopy(1) = [];
%             sched_func(curr) = [];
%             curr = 1;
%             iteration = true;
%         else
%             if iteration || delay == 1
%                 delay = delay +1;
%                 iteration = false;
%             end
%             curr = curr+1;
% 
%         end
%     end
% 
%     % if one period contains more than one route e.g sched_funcb = {[1 2 ; 2 3; 3 4], [4 5; 5 6]} sched_funcc = {[2 3; 1 2], [3 4 ;4 5], [5 6]};
%     %[sep_delay{:}] = deal(3) assigning same value to all
% else
%     for i = 1: length(paths)
%         pathsCopy = route_path(paths{i});
%         delay_all_paths{i} = [delaychecker(pathsCopy,sched_func), i]  ;
%     end
% end
% 
% 
% % Extract the first elements of each cell
% first_elems = cellfun(@(x) x(1), delay_all_paths);
% 
% % Sort the indices based on the first elements
% [~, idx] = sort(first_elems);
% 
% % Sort the cell array using the sorted indices
% sorted_A = delay_all_paths(idx);
% 
% least_delay = sorted_A{1}(1);
% max_delay = sorted_A{length(sorted_A)}(1);
% dic_delays = dictionary;
% dic_weight = dictionary;
% for curr_delay = least_delay:max_delay
%     first_elements = cellfun(@(x) x(1), sorted_A);
%     output_array = sorted_A(first_elements == curr_delay); %first = [3,2],[3,4][3,10]
% 
%     %now output array correspond to the path we could take in
%     %curr_delay/curr_hops. we want to check what paths.
%     second_element = cellfun(@(x) x(2), output_array);
%     nr_paths = length(second_element);
%     curr_paths = {};
%     for k = 1: nr_paths
%         curr_paths{k} = route_path(paths{second_element(k)});
%        
%     end
% 
% 
% 
%     seen_before = false(length(curr_paths),length(curr_paths{1}));
%     sz3 = size(seen_before);
%     for i = 1: sz3(1)
%         len = sz3(2);
%         for j = 1: len
%             if i <= length(curr_paths)
%                 curr = curr_paths{i}(1,j);
%                 curr_value = 1;
%                 ny_length = length(curr_paths);
%                 for kk = 1: ny_length
% 
%                     if (any(cellfun(@(c) isequal(c, curr{1}), curr_paths{curr_value})) && curr_value ~= i)
%                         seen_before(i,j) = true;
%                         curr_paths(curr_value) = [];
%                         curr_value= curr_value-1;
%                     else
%                         %not exist
%                     end
% 
%                     curr_value= curr_value+1;
%                 end
% 
%             end
% 
% 
% 
%         end
% 
%     end
% 
% 
%     dic_delays(curr_delay) = length(curr_paths);
%     %product of weight for p. 
%     total_weight = 0; 
%     for len_paths = 1:length(curr_paths)
%         curr_weight = 1; 
%         for len_p = 1:length(curr_paths{len_paths})
%              curr_route = curr_paths{len_paths}(len_p);
%              edge_idx = findedge(G, curr_route{1}(1), curr_route{1}(2));
%              edge_weight = G.Edges.Weight(edge_idx);
%              curr_weight = curr_weight * edge_weight;
%                 
% 
%         end
%         total_weight = total_weight+curr_weight;
%     end
% dic_weight(curr_delay) = total_weight;
%     
% end
% z = tf('z');
% transfer_function = 0;
% % iterate over keys
% possible_path = length(paths); % number of edges connected to destination node. 
% keys = dic_delays.keys;
% for i = 1:numel(keys)
%     key = keys(i);
%     value = dic_delays(key);
%     weight_value = dic_weight(key);
%     transfer_function = transfer_function + (weight_value/z^key)*(value/possible_path);
% end
% 
% 
% 
% transfer_function = minreal(transfer_function);

%sched_func = {[1 2; 1 3], [2 4 ; 2 5], [3 4 ;3 6], [4 5; 4 6], [5 7], [6 7]};
%different scheduling functions
%scheduling function 1: {(v1, v2)}, {(v2, v3)}, {(v3, v4)},{(v4, v5)} ,{(v5, v6)}
%sched_func = {[2 3], [3 4], [4 5], [1 2], [5 6]};
%sched_func = {[1 2],[2 3], [3 4], [4 5], [5 6]};
%scheduling function 2: {(v1, v2),(v2, v3),(v3, v4)},{(v4, v5),(v5, v6)}
%sched_func = {[1 2 ;1 3], [2 4;2 5;3 2],[3 5;4 5;4 6;5 6]}; latest

%sched_func = {[1 2;1 3;2 4; 3 4]};

% highlight(p, {'Vuc', 'V1', 'V3', 'V4'}, {'V1', 'V3', 'V4', 'Va'}, "EdgeColor","green");
% highlight(p, {'Vuc', 'V1', 'V3', 'V5'}, {'V1', 'V3', 'V5', 'Va'}, "EdgeColor","cyan");
% highlight(p, {'Vuc', 'V1', 'V4'}, {'V1', 'V4', 'Va'}, "EdgeColor","red");
% highlight(p, {'Vuc', 'V2', 'V3', 'V4'}, {'V2', 'V3', 'V4', 'Va'}, "EdgeColor","magenta");
% highlight(p, {'Vuc', 'V2', 'V3', 'V5'}, {'V2', 'V3', 'V5', 'Va'}, "EdgeColor","blue");
% highlight(p, {'Vuc', 'V2', 'V5'}, {'V2', 'V5','Va'}, "EdgeColor","yellow");
% highlight(p, {'Vc', 'V2', 'V4', 'V5'}, {'V2', 'V4', 'V5', 'Vu'}, "EdgeColor","green");
% highlight(p, {'Vc', 'V2', 'V4'}, {'V2', 'V4', 'Vu'}, "EdgeColor","blue");

%weight = 2*rand(1,10) - 1;
%weightR = [-0.020197222975552,-0.664145708635487,0.957361299282318,0.425388943357828,9.432483096860622e-04,-0.057823250916121,-0.880762264840722,0.363943808298126,-0.915137724998517,-0.857109070798715];
%weight = [100,-3, 1, 20];
%weight = [0.1660,   -0.4964,   -0.4191 ,   0.2342, 0.5670];
%weight = [0.043299684928567,-0.806539948438266,0.636297107719249,0.635094184158573,0.444879184733685,-0.700269115044066,0.319210505816614,0.037189885021076,0.945949109527725,0.297982985424712];

%names = {'V1','V2', 'V3', 'V4', 'V5','V6', 'V7','V8'}; 
%names = {'Vc', 'V2', 'V3','Vu'};

%s = [1,1,1,1,1,1,2,3,4,5,6,7]
%d = [2,3,4,5,6,7,8,8,8,8,8,8]; 

% 
% s = [1,1,2,2,3,3,4,4,5,6];
% d = [2,3,4,5,4,6,5,6,7,7];
%s = [1,1,2,3];
%d = [2,3,4,4];


% minimum_val = [1000    1000    1000    1000    1000    1000];  -0.0158 + 0.1812i  -0.0158 - 0.1812i   0.0955 + 0.1098i   0.0955 - 0.1098i  -0.7625 - 0.4174i  -0.7625 + 0.4174i
% old_max = -9999;
%     for i = 1:5
%         for k = 1:5
%             curr_max = sum(residual_max{i,k} < minimum_val);
%
%            if(curr_max >old_max)
%                 old_max = curr_max;
%                 minimum_val = residual_max{i,k};
%                 mod = i;
%                 contr = k;
%            end
%         end
%
%     end
%correct_model = [mod, contr];

%end



% for i = 1:10
%     figure;
%     hold on;
%     plot(abs(randomKResult{i}(:,1))/max(abs(randomKResult{i}(:,1))),'o', 'MarkerSize', 8, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black'); %A1K1-A1K
%     plot(abs(estimatedResult{model(i),ctlr(i)}(:,1))/max(abs(estimatedResult{model(i),ctlr(i)}(:,1))),'*', 'MarkerSize', 10, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black');
%     legend('Random', 'Estimated')
%     hold off;
% end
%
%
% for i = 1:10
%      figure;
% %     hold on;
%      random_Model_Controller = scaled_version(randomKResult{i}(:,1));
%     corr_estimated =  scaled_version(estimatedResult{5,2}(:,1));
%     disp(Model_closer_to_zero(random_Model_Controller,estimatedResult));
%     %disp([random_Model_Controller,corr_estimated]);
% %     plot(random_Model_Controller,'-', 'Color','blue'); %A1K1-A1K
% %     plot(corr_estimated,'-', 'Color','black');
%      plot(random_Model_Controller-corr_estimated,'-', 'Color','black');
% %     legend('Random', 'Estimated')
% %     hold off;
% end




% for i = 1:5
%     for k = 1:5
%         figure;
%
%         scaled_estimated = scaled_version(estimatedResult{i,k}(:,1));
%         scaled_estimated1 = scaled_version(measuredResult{i,k}(:,1));
%         plot(abs(scaled_estimated1)-abs(scaled_estimated),'-', 'Color','blue');
%
%
%     end
% end
%
% for i = 1:5
%     figure;
%     scaled= scaled_version(randomKResult{i}(:,1));
%     scaled_estimated = scaled_version(estimatedResult{i,2}(:,1));
%     plot(scaled-scaled_estimated);
% end

%model = [1,3,2,5,4,5,5,4,3,2];


% current_model = randomKResult{10};
% estimated_model_1 = scaled_version(estimatedResult{model(1),ctlr(2)});
% estimated_model_2 = scaled_version(estimatedResult{model(2),ctlr(2)});
% estimated_model_3 = scaled_version(estimatedResult{model(3),ctlr(2)});
% estimated_model_4 = scaled_version(estimatedResult{model(4),ctlr(2)});
% estimated_model_5 = scaled_version(estimatedResult{model(5),ctlr(2)});
%
%
% corrcoef_1 = corrcoef(current_model, estimated_model_1);
% corrcoef_2 = corrcoef(current_model, estimated_model_2);
% corrcoef_3 = corrcoef(current_model, estimated_model_3);
% corrcoef_4 = corrcoef(current_model, estimated_model_4);
% corrcoef_5 = corrcoef(current_model, estimated_model_5);
%
%
% fprintf('Correlation coefficients with current model:\n')
% fprintf('Estimated Model 1: %.3f\n', corrcoef_1(1,2))
% fprintf('Estimated Model 2: %.3f\n', corrcoef_2(1,2))
% fprintf('Estimated Model 3: %.3f\n', corrcoef_3(1,2))
% fprintf('Estimated Model 4: %.3f\n', corrcoef_4(1,2))
% fprintf('Estimated Model 5: %.3f\n', corrcoef_5(1,2))


% TF 

% function [tf_empty, transfer_function] = transferFunctionCalculate(paths,sched_func,G)
% 
% if isempty(paths) % if there is no path then no need of computing the transfer function
%     tf_empty = 1;
%     return;
% else
%     max_num_hops = max(cellfun(@(path) numel(path)-1, paths)); 
%     maxDR = max_num_hops;
%     pathsCopy = paths{1};
%     sched_length = length(sched_func);
%     delay = 1;
%     iteration = false;
%     % route_check = {};
%     % route_final ={} ;
%     % current = 1;
%     % sep_delay = {};
%     delay_all_paths = {};
%     
% 
%     if sched_length == maxDR % if the order of scheduling function is same as flow of path then in this case each hop can be done in different periods. 
%         curr = 1;
%         while ~isempty(sched_func)
%             if all(pathsCopy{1} == sched_func{curr})
%                 pathsCopy(1) = [];
%                 sched_func(curr) = [];
%                 curr = 1;
%                 iteration = true;
%             else
%                 if iteration || delay == 1
%                     delay = delay +1;
%                     iteration = false;
%                 end
%                 curr = curr+1;
% 
%             end
%         end
% 
%         % if one period contains more than one route e.g sched_funcb = {[1 2 ; 2 3; 3 4], [4 5; 5 6]} sched_funcc = {[2 3; 1 2], [3 4 ;4 5], [5 6]};
%         %[sep_delay{:}] = deal(3) assigning same value to all
%     else
%         % Loop through each path of the graph. 
%         for i = 1: length(paths)
%             pathsCopy = route_path(paths{i}); % Formating paths from 1 4 5 to cells {[1 4]} {[4 5]} corresponding to edges instead. 
%             delay_all_paths{i} = [delaychecker(pathsCopy,sched_func), i]; % Computing delays/Hops for current paths by using defined scheduling function
%         end
%     end
% 
%     %Sorting paths from least delay to maximum delay. 
%     first_elems = cellfun(@(x) x(1), delay_all_paths);
%     [~, idx] = sort(first_elems);
%     sorted_A = delay_all_paths(idx);
%     least_delay = sorted_A{1}(1);
%     max_delay = sorted_A{length(sorted_A)}(1);
%     
%     %Creating dictionary which includes delays and weights. 
%     dic_delays = dictionary;
%     dic_weight = dictionary;
%     for curr_delay = least_delay:max_delay
%         first_elements = cellfun(@(x) x(1), sorted_A);
%         output_array = sorted_A(first_elements == curr_delay); 
% 
%         %Now output array correspond to the path we could take in
%         %curr_delay/curr_hops. we want to check what paths has been taken to compute product of weights.
%         second_element = cellfun(@(x) x(2), output_array);
%         nr_paths = length(second_element);
%         curr_paths = {};
%         for k = 1: nr_paths
%             curr_paths{k} = route_path(paths{second_element(k)});
% 
%         end
% 
%         % As mentioned in report that each path can be taken only once
%         % thats why this section is created to check we don't have taken
%         % any duplicate paths beacuse it will otherwise impact the weights.
%         
%         seen_before = false(length(curr_paths),length(curr_paths{1}));
%         sz3 = size(seen_before);
%         for i = 1: sz3(1)
%             len = sz3(2);
%             for j = 1: len
%                 if i <= length(curr_paths)
%                     curr = curr_paths{i}(1,j);
%                     curr_value = 1;
%                     ny_length = length(curr_paths);
%                     for kk = 1: ny_length
% 
%                         if (any(cellfun(@(c) isequal(c, curr{1}), curr_paths{curr_value})) && curr_value ~= i)
%                             seen_before(i,j) = true;
%                             curr_paths(curr_value) = [];
%                             curr_value= curr_value-1;
%                         else
%                             %not exist
%                         end
% 
%                         curr_value= curr_value+1;
%                     end
% 
%                 end
% 
%             end
%         end
%         dic_delays(curr_delay) = length(curr_paths);
%         % Computes the product of weight for current path.
%         total_weight = 0;
%         for len_paths = 1:length(curr_paths)
%             curr_weight = 1;
%             for len_p = 1:length(curr_paths{len_paths})
%                 curr_route = curr_paths{len_paths}(len_p);
%                 edge_idx = findedge(G, curr_route{1}(1), curr_route{1}(2));
%                 edge_weight = G.Edges.Weight(edge_idx);
%                 curr_weight = curr_weight * edge_weight;
%             end
%             total_weight = total_weight+curr_weight; %Sum the weight of all paths.  
%         end
%         dic_weight(curr_delay) = total_weight;
% 
%     end
% 
%     z = tf('z');
%     transfer_function = 0;
%     possible_path = length(paths);  
%     keys = dic_delays.keys;
%     for i = 1:numel(keys)
%         key = keys(i);
%         value = dic_delays(key);
%         weight_value = dic_weight(key);
%         transfer_function = transfer_function + (weight_value/z^key)*(value/possible_path); % Defining the transfer function where numerator is the total weight and denominator is the z^delay. 
%     end
%     transfer_function = minreal(transfer_function); % Minimal realization of tranfer function.  
%     tf_empty = 0;
% end
% 
% end

for i = 1: length(model)
   subplot(2,3,i)
    time = linspace(0,steps);
    for j = 0:(length(t)/steps)-1
        [modelNr,controllerNr,currentResidual]  = closest_Model2(measuredResult_y,estimatedResult_y,model(i),ctlr(i), (j*steps)+1,steps);
        hold on
        txt = ['Model: ',num2str(i), ' C: ',num2str(j+1) ];
        %plot(begin_var:last_end,measuredResult_y{model(i),ctlr(i)}((j*steps)+1:((j*steps)+1+steps)-1,:)-estimatedResult_y{model(i),ctlr(i)}((j*steps)+1:((j*steps)+1+steps)-1,:),'DisplayName',txt)
        begin_var = last_end+1;
        last_end = last_end+steps;
        %title(['Model ', num2str(model(i)), ' Controller ', num2str(ctlr(i))]);
        xlabel('Time steps');
        ylabel('Residuals')
        disp(['Model: ', num2str(modelNr), ' Controller: ', num2str(controllerNr)]);
        if(modelNr ~= controllerNr)
            disp('Current Model is running with wrong controller');
            ctlr(i) = modelNr;
        end
        curr_residual{i,j+1} = currentResidual;
    end
end

ctlr = [1,2,3,1,3,1]%,3,3,2,1,2,2,2];
model = [1,1,3,2,1,3]%,3,1,2,3,1,2,3];
curr_residual = {};
for i = 1: length(model)
    subplot(2,3,i)
    time = linspace(0,steps);
    for j = 0:(length(t)/steps)-1
        [modelNr,controllerNr,currentResidual]  = closest_Model2(measuredResult_y,estimatedResult_y,model(i),ctlr(i), (j*steps)+1,steps);
        title(['Model ', num2str(model(i)), ' Controller ', num2str(ctlr(i))]);
        xlabel('Time steps');
        ylabel('Residuals')
        disp(['Model: ', num2str(modelNr), ' Controller: ', num2str(controllerNr)]);
        if(modelNr ~= controllerNr)
            disp('Current Model is running with wrong controller');
            
        end
        curr_residual{i,j+1} = currentResidual;
    end
end
begin = 1; 
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
for i = 1: size(estimateddata,1)
     %Color = ['b','r','c'];
    for j = 1:size(estimateddata,2)
        residual=currentdata{model,controller}(crr_step:(crr_step+steps)-1,:)-estimateddata{i,j}(crr_step:(crr_step+steps)-1,:);
        mean_e1= abs(mean(residual));
        hold on
        txt = ['Model: ',num2str(i), ' C: ',num2str(j) ];
        plot(crr_step:(crr_step+steps)-1,residual,'DisplayName',txt)
        
        legend();
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
            
            xhat = x;
        end
    end
end
%mode detection
%multimode diagnosis on switched affine system with noisy measurement


