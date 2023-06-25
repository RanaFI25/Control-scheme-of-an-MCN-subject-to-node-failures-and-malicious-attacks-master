
clear all
clc
close all

 s = [1,1,1,2,3,4]; %  Edges from for example from 1 to 2 then one edge is created from node 1 to 2.
 d = [2,3,4,5,5,5]; % Edges to


%s = [1,1,2,2,2,3,3,4,4,5,5,5,6,6,7,8,8,8,9]; %  Edges from for example from 1 to 2 then one edge is created from node 1 to 2.
%d = [2,3,3,4,5,5,6,5,7,7,8,9,5,9,10,7,10,9,10]; % Edges to

 weight_R = [-0.0264   -0.1283   -0.1064   -0.3873    0.0170    0.0215];
 weight_O = [0.6353    0.5897    0.2886   -0.2428    0.6232    0.0657];
% 
 [paths_R,G_R] = createGraph(s,d,[4,3],weight_R);
 [paths_O,G_O] = createGraph(s,d,50,weight_O);

%[paths_R,G_R] = createGraph(s,d,[2,3]); %Creating reachability graph where node 2 is attacked.
%[paths_O,G_O] = createGraph(s,d,[3,4]); %Creating observability graph where node 3 is attacked.

% Different scheduling function for assigns transmission to nodes.
sched_func = {[1,2],[1,3],[1,4],[2,5],[3,5],[4,5]};
%sched_func = {[1,2],[1,3],[2,3],[2,4],[2,5],[3,5],[3,6],[4,5],[4,7],[5,7],[5,8],[5,9],[6,5],[6,9],[7,10],[8,7],[8,9],[8,10],[9,10]};
%sched_func = {[1,2 ;1,3;2,3; 2,4; 2,5; 3,5; 3,6;4,5;4,7;5,7;5,8;5,9;6,5;6,9;7,10;8,7;8,9;8,10;9,10]};
%sched_func_2 = {[1 2 ;1 3;1 4;2 5; 3 5;4 5]};
%sched_func_3 = {[1 2 ;1 3;1 4], [2 5], [3 5],[4 5]};

%start with delay equal to one and check if in thi  s delay message can be delivered to dest: delay 1 means that in every
%period there should be maximum one hop in above scheduling functions only first satisfy that%


%Computing transfer function for both graphs
[tf_empty_R , transfer_function_R] = transferFunctionCalculate(paths_R,sched_func,G_R);
[tf_empty_O , transfer_function_O] = transferFunctionCalculate(paths_O,sched_func,G_O);

%Computing transfer function of M(z).
if ~tf_empty_R && ~tf_empty_O % if transfer function is not zero
    z = tf('z');
    %Accessing transfer function data in order to check the zeros and poles
    %of it.
    [num_R, den_R] = tfdata(transfer_function_R, 'v');
    [num_O, den_O] = tfdata(transfer_function_O, 'v');
    % R_tf = tf(num_R, den_R, 1, 'Variable', 'z');
    num_p = 0.5;
    den_p = [1, -1.1];

    %Defining numerator and denominator of the plant which can than be
    %coverted to the discrete time transfer function.
    %num_p =  [0     0     1];
    %den_p = [1    1    5];
    plant_tf = tf(num_p, den_p, 1, 'variable', 'z');
    % [A_0,B_0,C_0,D_0] = tf2ss(num_p,den_p);

    %R_tf = -0.00513/z;


    R_tf = transfer_function_R;
    O_tf = transfer_function_O;
    %O_tf = transfer_function;

    % num_O = [1];
    % den_O = [1 0];
    % O_tf = tf(num_O, den_O, 1, 'Variable', 'z');
    % zeros and poles
    [z_O,p_O,k_O] = tf2zp(num_O,den_O);
    [z_R,p_R,k_R] = tf2zp(num_R,den_R);
    [z_P,p_P,k_P] = tf2zp(num_p,den_p);

    %Creating transfer function of Multihop control network
     M_z_tf = R_tf*plant_tf*O_tf;
    [num_m, den_m] = tfdata(M_z_tf, 'v');
    [z_m,p_m,k_m] = tf2zp(num_m,den_m);

    % Converting transfer function to state space representation.
    [A3,B3,C3,D3] = tf2ss(num_m,den_m);

end




