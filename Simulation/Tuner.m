clear all;
clc;

%===================================================
% Tuner for the loops using linearised state space system
%===================================================

% ----------------------------------------------------------------------------------------
%INITIAL DATA HELICOPTER
g=9.81;	
cla=2*pi; %Airfoil theory [1/rad]
volh=.0722;	%blade solidity	
lok=9.334;
cds=2.415;
mass=4536;
W = mass*g;
rho=1.225;
diam=13.41;
omega = 324*2*pi/60;
vtip= omega*(diam/2);
iy=(1/12)*mass*(4.1^2 + 16.14^2); % Rectangular box approximation with a = 4.1 and b = 16.41 [kg/m^4]
mast=2.285; %(m) Vertical distance between rotor CG and rotor hub
area=pi/4*diam^2;
tau=.1;		%time constant in dynamiCs inflow!!!
collect0=6*pi/180; % (rad) Collective pitch
longit0=0*pi/180; % (rad) Longitudinal cyclic pitch angle
[Cdf, S] = fuselage_dragS();

% ---------------------------------Parameters------------------------------------
%Initial values;
t0=0; %(sec) Setting initial time 
u0=90*0.51444; %(m/sec) Setting initial helicopter airspeed component along body x-axis 
w0=0; %(m/sec) Setting initial helicopter airspeed component along body z-axis
q0=0; %(rad/sec) Setting initial helicopter pitch rate 
pitch0=-3*pi/180; %(rad) Setting initial helicopter pitch angle
x0=0; %(m) Setting initial helicopter longitudinal position 
V0 = sqrt(u0^2 + w0^2);
labi0=lambda_i(V0, omega, diam/2, rho, S, Cdf, W); %Initialization non-dimensional inflow in hover!!!

%===================================================
% Linearised state space
%===================================================

% Symbolic variables
syms u w pitch q collect longit labi real

% !! PARAMETERS !!
% mu and lambda_c
qdiml=q/omega;
vdiml=sqrt(u^2+w^2)/vtip;
phi=atan(w/u);
vv_sym = sqrt(u^2 + w^2);
alfc=longit-phi;
mu=vdiml*cos(alfc); %Vcos(alphac)/omega*R
labc=vdiml*sin(alfc); %Vsin(alphac)/omega*R

ctelem=cla*volh/4*(2/3*collect*(1+1.5*mu^2)-(labc+labi));
labidot=ctelem; 
thrust=labidot*rho*vtip^2*area;
stap = 0.01;

%a1 Flapping calculi
teller=-16/lok*qdiml+8/3*mu*collect-2*mu*(labc+labi);
a1=teller/(1-.5*mu^2);

%Thrust coefficient from Glauert
alfd=alfc-a1; % alpha_d
ctglau=2*labi*sqrt((vdiml*cos(alfd))^2+(vdiml*sin(alfd)+labi)^2);

helling=longit-a1;
vv=vdiml*vtip; %it is 1/sqrt(u^2+w^2)

% State vector including labi
x_sym = [u; w; pitch; q; labi];  % Including labi as a state

% Symbolic non-linear dynamics
udot_sym     = -g*sin(pitch) - cds/mass * 0.5*rho*u*vv + thrust/mass*sin(helling) - q*w;
wdot_sym     =  g*cos(pitch) - cds/mass * 0.5*rho*w*vv - thrust/mass*cos(helling) + q*u;
pitchdot_sym = q;
qdot_sym     = -thrust*mast/iy * sin(helling);
labi_dot_sym = (ctelem - ctglau) / tau;  % Dynamics for labi

% Dynamics vector including labi_dot
f = [udot_sym; wdot_sym; pitchdot_sym; qdot_sym; labi_dot_sym];  % Updated dynamics with labi

% Jacobians
A_sym = jacobian(f, x_sym);  % Jacobian with respect to x including labi
B_sym = jacobian(f, [collect; longit]);  % Jacobian with respect to inputs (collective, longitudinal cyclic)

% Trim values (substitution list)
subs_list = [u; w; pitch; q; collect; longit; labi];  % Include labi in the substitution list
vals_list = [u0; w0; pitch0; q0; collect0; longit0; labi0];  % Trim values including labi

% Substitute trim values into the Jacobians
A_lin = double(subs(A_sym, subs_list, vals_list));  % Linearized A matrix with labi
B_lin = double(subs(B_sym, subs_list, vals_list));  % Linearized B matrix with labi

% Output matrices
C = eye(5);
D = zeros(5,2);

sys_lin = ss(A_lin, B_lin, C, D);

% % Optional: display matrices
disp('A matrix:'); disp(A_lin);
disp('B matrix:'); disp(B_lin);

output_idx = 4;
input_idx = 2;

sys_q = tf(sys_lin(output_idx, input_idx));  % System transfer function

% Define the PI controller
K_d = 1;  % Derivative gain
K_i = 0;    % Integral gain
C = K_d;  % Integral controller: K_d + K_i/s

figure;
rlocus(-sys_q); 

% % Plot the step response
% step(CL_q_PI);
% title('Step Response with PI Controller');
% grid on;

