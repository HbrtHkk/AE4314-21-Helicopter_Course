clear; 

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
collect(1)=6*pi/180; % (rad) Collective pitch
longit(1)=0*pi/180; % (rad) Longitudinal cyclic pitch angle
[Cdf, S] = fuselage_dragS();

% ---------------------------------Parameters------------------------------------
%Trim settings;
t0=0; %(sec) Setting initial time 
u0=90*0.51444; %(m/sec) Setting initial helicopter airspeed component along body x-axis 
w0=0; %(m/sec) Setting initial helicopter airspeed component along body z-axis
q0=0; %(rad/sec) Setting initial helicopter pitch rate 
pitch0=-2*pi/180; %(rad) Setting initial helicopter pitch angle
x0=0; %(m) Setting initial helicopter longitudinal position 
V0 = sqrt(u0^2 + w0^2);
labi0= lambda_i(V0, omega, diam/2, rho, S, Cdf, W); %Initialization non-dimensional inflow in hover!!!

t(1)=t0;
u(1)=u0;
w(1)=w0;
q(1)=q0;
pitch(1)=pitch0;
x(1)=x0;
labi(1)=labi0;
z(1)=0;
V_tot(1) = V0;

% INTEGRATION 
aantal=1000;
teind=20;
stap=(teind-t0)/aantal;
int_error = 0;
int_errorc =0;
cdes = u0*sin(pitch0)-w0*cos(pitch0);

%  -------------------------------Start of Simulation------------------------------------
for i=1:aantal 
   if t(i)>=0.5 & t(i)<=1 
       longit(i)=1*pi/180;
       collect(i)= 0*pi/180 + collect(1);
       
   % elseif t(i)>120
   %     % Updating trim position
   %     longit(1) = longit(i-1);
   %     collect(1) = collect(i-1);
   % 
   %     % New forward velocity input
   %     u(1) = 70*0.5144;
   % elseif t(i)>240
   %     % Updating trim position
   %     longit(1) = longit(i-1);
   %     collect(1) = collect(i-1);
   % 
   %     % New forward velocity input
   %     u_des(1) = 90*0.5144;  
   % 
   % elseif t(i)>360
   %     % Updating trim position
   %     longit(1) = longit(i-1);
   %     collect(1) = collect(i-1);
   % 
   %     % New forward velocity input
   %     u0 = 110*0.5144;  
   else 
       longit(i) = 0*pi/180;
       collect(i) = collect(1);
   end
   %===================================================
   % Controller always on
   %===================================================
   if t(i)>=1

       % Compute the error term
       error = pitch0-pitch(i);
       errorc = cdes-c(i);
       int_error = int_error + error * stap;  % Integral of the error

       K1 = 0; % Proportional term for angle
       K2 = -8.75;  % Derivative term (rate) - reduced from -0.85
       K3 = 0; % Integral term - start small
       longit(i) = (K1 + K3*int_error)*(0-(pitch0 - pitch(i))) + (K2)*(0-(q0 - q(i))); %PD 	%In rad

       % K4 = 0.0075;
       % int_errorc = int_errorc + errorc * stap;
       % K5 = 0;  %-0.0001;
       % %cdes = K5*(hdes - h);
       % collect(i)= collect(1) + (K4 + K5*int_error)*(cdes - c(i));
   end   
 
    c(i)=u(i)*sin(pitch(i))-w(i)*cos(pitch(i));
    h(i)=-z(i);

   %===================================================
   % Parameters
   %===================================================
    qdiml(i)=q(i)/omega;
    vdiml(i)=sqrt(u(i)^2+w(i)^2)/vtip;
    
    % If u velocity component is 0 
    % phi is the inflow angle!!!
    if u(i)==0 	
        if w(i)>0 	
            phi(i)=pi/2;
        else 
            phi(i)=-pi/2;
        end
    else
        phi(i)=atan(w(i)/u(i));
    end
    
    % Accounting for backwards velocity?
    if u(i)<0
        phi(i)=phi(i)+pi;
    end 
    
    % Angle of attack 
    alfc(i)=longit(i)-phi(i);
    
    % mu and lambda_c
    mu(i)=vdiml(i)*cos(alfc(i)); %Vcos(alphac)/omega*R
    labc(i)=vdiml(i)*sin(alfc(i)); %Vsin(alphac)/omega*R
    
    %a1 Flapping calculi
    teller(i)=-16/lok*qdiml(i)+8/3*mu(i)*collect(i)-2*mu(i)*(labc(i)+labi(i));
    a1(i)=teller(i)/(1-.5*mu(i)^2);
    
    %Thrust coefficient 
    ctelem(i)=cla*volh/4*(2/3*collect(i)*(1+1.5*mu(i)^2)-(labc(i)+labi(i)));

    %Thrust coefficient from Glauert
    alfd(i)=alfc(i)-a1(i); % alpha_d
    ctglau(i)=2*labi(i)*sqrt((vdiml(i)*cos(alfd(i)))^2+(vdiml(i)*sin(alfd(i))+labi(i))^2);
    
   %===================================================
   % Equations of Motion
   %===================================================
    labidot(i)=ctelem(i); 
    thrust(i)=labidot(i)*rho*vtip^2*area;
    helling(i)=longit(i)-a1(i);
    vv(i)=vdiml(i)*vtip; %it is 1/sqrt(u^2+w^2)
    D = 0.5*rho*(V_tot(i)^2)*cds*S;

    % udot
    udot(i)=-g*sin(pitch(i))-cds/mass*.5*rho*u(i)*vv(i)+...
    thrust(i)/mass*sin(helling(i))-q(i)*w(i);
    
    % wdot
    wdot(i)=g*cos(pitch(i))-cds/mass*.5*rho*w(i)*vv(i)-...
    thrust(i)/mass*cos(helling(i))+q(i)*u(i);
    
    %qdot 
    qdot(i)=-thrust(i)*mast/iy*sin(helling(i));
    
    % theta_f
    pitchdot(i)=q(i);
    
    % Change in longitudinal position
    xdot(i)=u(i)*cos(pitch(i))+w(i)*sin(pitch(i));
    
    % Change in altitude
    zdot(i)=-c(i);

    % Change in lambda_i
    labidot(i)=(ctelem(i)-ctglau(i))/tau;
    labi(i+1)=labi(i)+stap*labidot(i);
    u(i+1)=u(i) + stap*udot(i);
    w(i+1)=w(i) + stap*wdot(i);

    corrdot(i)=cdes-c(i);

    c(i+1) = c(i) + stap*corrdot(i);
    q(i+1)= q(i) + stap*qdot(i);
    pitch(i+1)=pitch(i) + stap*pitchdot(i);
    x(i+1)=x(i) +stap*xdot(i);
    z(i+1)=z(i) + stap*zdot(i);
    t(i+1)=t(i) +stap;
    V_tot(i+1) = sqrt(u(i+1)^2 + w(i+1)^2);
end;

figure;
subplot(2, 2, 1);
plot(t,c),xlabel('t (s)'),ylabel('c'),grid;
subplot(2, 2, 2);
plot(t(1:aantal),collect*180/pi),xlabel('t (s)'),ylabel('collective (deg)'),grid;
subplot(2, 2, 3);
plot(t, V_tot),xlabel('t (s)'),ylabel('Total airspeed V (m/s)'),grid;
subplot(2, 2, 4);
plot(t, w),xlabel('t (s)'),ylabel('w (m/s)'),grid, pause;

figure;
subplot(1, 3, 1);
plot(t,pitch*180/pi),xlabel('t (s)'),ylabel('Pitch (deg)'),grid;
subplot(1, 3, 2);
plot(t,q*180/pi),xlabel('t (s)'),ylabel('Pitch rate q (deg/s)'),grid;
subplot(1, 3, 3);
plot(t, u),xlabel('t (s)'),ylabel('u'),grid, pause;

plot(t(1:aantal),longit*180/pi),xlabel('t (s)'),ylabel('Cyclic input (deg)'),grid; 

% figure;
% subplot(1, 2, 1);
% plot(t,u),xlabel('t (s)'),ylabel('Horizontal airspeed u (m/s)'), grid;
% subplot(1, 2, 2);
% plot(t,q*180/pi),xlabel('t (s)'),ylabel('Pitch rate (deg/s)'),grid, pause;

% figure;
% subplot(1, 2, 1);
% %Longitudinal cyclic 

%subplot(1, 2, 2);
%plot(t(1:aantal),collect*180/pi),xlabel('t (s)'),ylabel('Collective input (deg)'),grid, pause;
% 
% figure;
% subplot(1, 3, 1);
% plot(t,u),xlabel('t (s)'),ylabel('Horizontal airspeed u (m/s)'), grid;
% 
% subplot(1, 3, 2);
% plot(t,w),xlabel('t (s)'),ylabel('Vertical airspeed w (m/s)'),grid;
% 
% subplot(1, 3, 3);
% plot(t, V_tot),xlabel('t (s)'),ylabel('Total airspeed V (m/s)'),grid;
% 
% figure;
% subplot(1, 2, 1);
% plot(t,-z),xlabel('t (s)'),ylabel('Altitude h (m)'),grid;
% subplot(1, 2, 2);
% plot(t,w),xlabel('t (s)'),ylabel('Vertical airspeed w (m/s)'),grid, pause; 
% 
% figure;
% plot(t,x),xlabel('t (s)'),ylabel('Horizontal distance covered x(m)'),grid,pause;

