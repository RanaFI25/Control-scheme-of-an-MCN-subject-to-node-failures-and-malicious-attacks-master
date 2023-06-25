function [Yc,Tc,Yd,Td,Uc]=intersample_sim(sys_c,Fr_d,Fy_d,Tend)
%
% [Yc,Tc,Yd,Td,Uc]=intersample_sim(sys_c,Fr_d,Fy_d,Tend)
%
% Simulates the step response of a continuous-time system controlled
% by a sampling controller, including the intersample behaviour.
% Based on the control law u = Fr_d(q)r - Fy_d(q)y, and ZOH sampling.
% 
% Outputs:
% Yc - vector with the cont.-time outputs, at time instants in Tc
% Tc - vector with cont.-time instants
% Yd - vector with sampled outputs
% Td - vector with the samplings instants
% Uc - vector with control input, at time instants in Tc
% 
% Inputs:
% sys_c - LTI object representing cont.-time system
% Fr_d - LTI object representing the discr.-time controller filter Fr_d(q)
% Fy_d - LTI object representing the discr.-time controller filter Fy_d(q)
% Tend - Simulation time, default value = 8 [s]

% HN 2015-06-18
if nargin<4
    Tend=8;
end


F_d=minreal([Fr_d -Fy_d]); % Merge discr.-time controller into one LTI object 
[F,G,H,D,h]=ssdata(F_d);   % and extract its parameters.

N=floor(Tend/h); % Number of samples to be simulated.
Tend=N*h;        % Set end time to integer number of sample periods.

dt=0.01;
P=floor(h/dt);
dt=h/P;

sys_cd=c2d(sys_c,dt);
[AA,BB,CC]=ssdata(sys_cd);

% Initialize
xsys=zeros(size(BB));
Yc=[];
Yd=0;
Tc=[];
Uc=[];
xz=zeros(2,1);
r=1;
y=0;

% Loop for every sampling interval
for k=0:N-1
    u=H*xz+D*[r;y]; % Compute control input,
    xz=F*xz+G*[r;y]; % and update the controller.
    for l=1:P % Simulate the cont.-time system during the sampling interval.
        xsys=AA*xsys+BB*u;
        y=CC*xsys;
        Yc=[Yc; y];
        Tc=[Tc;k*h+l*dt];
        Uc=[Uc; u];
    end
    Yd=[Yd y]; % Store the sampled output.
end
Td=0:h:Tend; % Sampling time vector.