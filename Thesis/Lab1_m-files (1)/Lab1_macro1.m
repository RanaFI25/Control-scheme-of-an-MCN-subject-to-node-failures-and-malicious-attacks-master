% Set system parameters for cont.-time state space model:
A=[-1 -4; 1 0];
B=[1; 0];
C=[0 1];

sys_c=ss(A,B,C,0); % LTI-object for cont.-time system.


%h=0.1; % Sampling period.

% Discr.-time model by zero-order hold (ZOH):
sys_d=c2d(sys_c,h); 

% Cont.-time controller deisgn: Observer based state feedback:
% u = -Lc*xhat + mc*r, xhat estimate of x, from observer.
pc=2*[-1+i -1-i]; % Desired poles for cont.-time closed loop system.
Lc=place(A,B,pc); % State feedback gain, obtained by pole placement.
Gc=ss(A-B*Lc,B,C,0); % Theoretical cont.-time closed loop system.
mc=1/dcgain(Gc); % Gain to achieve static gain =1 from r to y.
Gc_c=mc*Gc; % Modified closed loop system, with static gain = 1.

Kc=place(A',C',2*pc)'; % Observer gain, by pole placement.
                       % Observer poles chosen twice as fast 
                       % as closed loop poles.

% The control law can be represented as
% U(s) = Fr(s)*R(s) - Fy(s)*Y(s).
Fy_c=ss(A-B*Lc-Kc*C,Kc,Lc,0); % Feedback filter.
Fr_c=mc*ss(A-B*Lc-Kc*C,B,-Lc,1); % Pre-filtering of reference.

% Construct sampling, discr.-time controllers.
% Method 1: Discretize the cont.-time controller:
Fy_tus=c2d(Fy_c,h,'tustin'); % Tustin's method used here.
Fr_tus=c2d(Fr_c,h,'tustin'); % You can try other methods as well ...
Gc_tus=Fr_tus*feedback(sys_d,Fy_tus); % Discr.-time closed loop system

% Method 2: Sample the system, by ZOH-sampling, and (re-) design the
% controller based on the discr.-time model:
pd_L=exp(h*pc); % Desired closed loop poles in discr.-time,
pd_K=exp(2*h*pc); % and corresponding observer poles.

[Ad,Bd,Cd]=ssdata(sys_d); % Extract discr.-time model parameters.
Ld=place(Ad,Bd,pd_L); % State feedback gain, by pole placement,
Kd=place(Ad',Cd',pd_K)'; % and observer gain.

Gc_d=ss(Ad-Bd*Ld,Bd,Cd,0,h); % Theoretical discr.-time closed loop system
md=1/dcgain(Gc_d); % Gain to achieve static gain = 1 from r to y

% Create LTI objects for feedback and pre-filters:
Fy_zoh=ss(Ad-Bd*Ld-Kd*Cd,Kd,Ld,0,h);
Fr_zoh=md*ss(Ad-Bd*Ld-Kd*Cd,Bd,-Ld,1,h);

Gc_zoh=Fr_zoh*feedback(sys_d,Fy_zoh); % Discr.-time closed loop system.

% Plot step responses of closed loop systems
figure(1)
step(Gc_c,Gc_tus,Gc_zoh,8);
legend('show','Location','southeast');