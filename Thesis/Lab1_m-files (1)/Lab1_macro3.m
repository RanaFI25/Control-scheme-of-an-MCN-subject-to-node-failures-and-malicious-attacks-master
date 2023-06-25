% Set system parameters for cont.-time state space model:
A=[-1 -4; 1 0];
B=[1; 0];
C=[0 1];

sys_c=ss(A,B,C,0); % LTI-object for cont.-time system.


%h=0.1; % Sampling period.

% Discr.-time model by zero-order hold (ZOH):
sys_d=c2d(sys_c,h); 


% Method 2: Sample the system, by ZOH-sampling, and (re-) design the
% controller based on the discr.-time model:
pd_L=[]; % Desired closed loop poles in discr.-time,
pd_K=[]; % and corresponding observer poles.

[Ad,Bd,Cd]=ssdata(sys_d); % Extract discr.-time model parameters.
Ld=place(Ad,Bd,pd_L); % State feedback gain, by pole placement,
Kd=place(Ad',Cd',pd_K)'; % and observer gain.

Gc_d=ss(Ad-Bd*Ld,Bd,Cd,0,h); % Theoretical discr.-time closed loop system
md=1/dcgain(Gc_d); % Gain to achieve static gain = 1 from r to y

% Create LTI objects for feedback and pre-filters:
Fy_db=ss(Ad-Bd*Ld-Kd*Cd,Kd,Ld,0,h);
Fr_db=md*ss(Ad-Bd*Ld-Kd*Cd,Bd,-Ld,1,h);

Gc_db=minreal(Fr_db*feedback(sys_d,Fy_db)); % Discr.-time closed loop system.

% Plot step responses of closed loop systems
figure(3)
step(Gc_db,8);
