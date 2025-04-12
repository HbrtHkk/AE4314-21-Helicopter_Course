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
R= diam/2;
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
cdes0 = u0*sin(pitch0)-w0*cos(pitch0);
c0 = cdes0;
%=================
% ==================================
% Linearised state space
%===================================================

% % Symbolic variables
% syms u w pitch q collect longit labi c real
% 
% % !! PARAMETERS !!
% % mu and lambda_c
% qdiml=q/omega;
% vdiml=sqrt(u^2+w^2)/vtip;
% phi=atan(w/u);
% vv_sym = sqrt(u^2 + w^2);
% alfc=longit-phi;
% mu=vdiml*cos(alfc); %Vcos(alphac)/omega*R
% labc=vdiml*sin(alfc); %Vsin(alphac)/omega*R
% 
% ctelem=cla*volh/4*(2/3*collect*(1+1.5*mu^2)-(labc+labi));
% labidot=ctelem; 
% thrust=labidot*rho*vtip^2*area;
% stap = 0.01;
% 
% %a1 Flapping calculi
% teller=-16/lok*qdiml+8/3*mu*collect-2*mu*(labc+labi);
% a1=teller/(1-.5*mu^2);
% 
% %Thrust coefficient from Glauert
% alfd=alfc-a1; % alpha_d
% ctglau=2*labi*sqrt((vdiml*cos(alfd))^2+(vdiml*sin(alfd)+labi)^2);
% 
% helling=longit-a1;
% vv=vdiml*vtip; %it is 1/sqrt(u^2+w^2)
% 
% % State vector including labi
% x_sym = [u; w; pitch; q; labi; c];  % Including labi as a state
% inputs = [collect; longit];
% % Symbolic non-linear dynamics
% udot_sym     = -g*sin(pitch) - cds/mass * 0.5*rho*u*vv + thrust/mass*sin(helling) - q*w;
% wdot_sym     =  g*cos(pitch) - cds/mass * 0.5*rho*w*vv - thrust/mass*cos(helling) + q*u;
% pitchdot_sym = q;
% qdot_sym     = -thrust*mast/iy * sin(helling);
% labi_dot_sym = (ctelem - ctglau) / tau;  % Dynamics for labi
% corrdot=cdes0-c;


%% Linearize non linear system equations
syms u w q theta_f lambda_i theta_0 theta_c

V = sqrt(u^2 + w^2);
alpha_c = theta_c - atan(w/u); % not sure about this
mu = V/(omega*R)*cos(alpha_c);
lambda_c = V*sin(alpha_c)/(omega*R);
a1 = (8/3*mu*theta_0 - 2*mu*(lambda_c + lambda_i)-16/lok*q/omega)/(1-1/2*mu^2);
C_Tglau = 2*lambda_i*sqrt((V/(omega*R)*cos(alpha_c-a1))^2+(V/(omega*omega)*sin(alpha_c-a1)+lambda_i)^2);
C_TBEM = 1/4*cla*volh*(2/3*theta_0*(1+3/2*mu^2)-(lambda_c+lambda_i));
T = C_TBEM*rho*(omega*omega)^2*pi*omega^2;
D = Cdf*1/2*rho*V^2*S;

X = T/mass*sin(theta_c-a1)-D/mass*u/V-g*sin(pitch0);
dX = jacobian(X, [u, w, q, theta_f, lambda_i, theta_0, theta_c]);
dX_new = double(subs(dX, [u, w, q, theta_f, lambda_i, theta_0, theta_c], ...
    [u0, w0, q0, pitch0, labi0, collect0, longit0]));

Z = -T/mass*cos(theta_c-a1)-D/mass*w/V+g*cos(pitch0);
dZ = jacobian(Z, [u, w, q, theta_f, lambda_i, theta_0, theta_c]);
dZ_new = double(subs(dZ, [u, w, q, theta_f, lambda_i, theta_0, theta_c], ...
    [u0, w0, q0, pitch0, labi0, collect0, longit0]));

M = -T/iy*mast*sin(theta_c-a1);
dM = jacobian(M, [u, w, q, theta_f, lambda_i, theta_0, theta_c]);
dM_new = double(subs(dM, [u, w, q, theta_f, lambda_i, theta_0, theta_c], ...
    [u0, w0, q0, pitch0, labi0, collect0, longit0]));

tau = 0.1;
Lambda = 1/tau*(C_TBEM-C_Tglau)/(omega*omega);
dLambda = jacobian(Lambda, [u, w, q, theta_f, lambda_i, theta_0, theta_c]);
dLambda_new = double(subs(dLambda, [u, w, q, theta_f, lambda_i, theta_0, theta_c], ...
    [u0, w0, q0, pitch0, labi0, collect0, longit0]));

%% Linear system
A = [dX_new(1) dX_new(2) (dX_new(3)-w0) (dX_new(4)) dX_new(5);
    dZ_new(1) dZ_new(2) (dZ_new(3)+u0) (dZ_new(4)) dZ_new(5);
    dM_new(1) dM_new(2) dM_new(3) dM_new(4) dM_new(5);
    0 0 1 0 0;
    dLambda_new(1) dLambda_new(2) dLambda_new(3) dLambda_new(4) dLambda_new(5)];
B = [dX_new(6) dX_new(7);
    dZ_new(6) dZ_new(7);
    dM_new(6) dM_new(7);
    0 0; 
    dLambda_new(6) dLambda_new(7)];
C = eye(5,5);
mdl = ss(A,B,C,0);

% Output matrix to select pitch rate (q = 4th state)
C_q = [0 0 0 1 0];

% D matrix for SISO system (just one output and one input)
D_q = 0;

% Create SISO system: cyclic input → pitch rate
sys_q = ss(A, B(:,2), C_q, D_q);

figure;
K =64.5;
rlocus(-K*sys_q);

% %% Verifying linear model
% t_end = 50;
% t = linspace(0,t_end, t_end/0.1+1);
% u = [ones(1,10)*1*pi/180, zeros(1,length(t)-10); zeros(1,length(t))];
% y = lsim(mdl, u, t);
% 
% t = linspace(0,t_end, t_end/0.1+1);
% u = [ones(1,10)*1*pi/180, zeros(1,length(t)-10); zeros(1,length(t))];
% y = lsim(mdl, u, t);
% 
% figure;
% subplot(2, 3, 1);
% plot(t, y(:, 1));
% legend('u')
% 
% subplot(2, 3, 2);
% plot(t, y(:, 2));
% legend('w')
% 
% subplot(2, 3, 3);
% plot(t, y(:, 3)*180/pi);
% legend('q [deg/s]')
% 
% subplot(2, 3, 4);
% plot(t, y(:, 4)*180/pi);
% legend('Pitch [deg]')
% 
% subplot(2, 3, 5);
% plot(t, y(:, 5));
% legend('labi')


%----------------------------------------------------------------------OLD
% % Dynamics vector including labi_dot
% f = [udot_sym; wdot_sym; pitchdot_sym; qdot_sym; labi_dot_sym; corrdot];  % Updated dynamics with labi
% 
% % Jacobians
% A_sym = jacobian(x_sym, x_sym);  % Jacobian with respect to x including labi
% B_sym = jacobian(inputs, [collect; longit]);  % Jacobian with respect to inputs (collective, longitudinal cyclic)
% 
% % Trim values (substitution list)
% subs_list = [u; w; pitch; q; collect; longit; labi; c];  % Include labi in the substitution list
% vals_list = [u0; w0; pitch0; q0; collect0; longit0; labi0; c0];  % Trim values including labi
% 
% % Substitute trim values into the Jacobians
% A_lin = double(subs(A_sym, subs_list, vals_list));  % Linearized A matrix with labi
% B_lin = double(subs(B_sym, subs_list, vals_list));  % Linearized B matrix with labi
% 
% % Output matrices
% C = eye(6);
% D = zeros(6,2);



