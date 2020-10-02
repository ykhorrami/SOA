%%% In The Name of Allah
%fourth order Runge-kutta finite difference solution for SOA  
clc;
clear all;
close all;
t0 =0.05e-12;      % initial time (starting time)
H0 =0;             % initial condition (H(t0))
h =0.05e-12;       %simulation steps
tp=5e-12;          %FWHM
tau0=tp/1.7627;
tauc=200e-12;       %carrier recovery time
alfa=5;            %linewidth enhancement factor
esat=10e-12;   
ein=15e-15;         %input energy
L=6.38336677e-3;   %length of SOA
g0=1096.6/2;         %linear gain(1/m)
h0=g0*L;              %unsaturated gain
tend=40;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hr(1) = h0;
tr(1) = t0;
iter=3200;         %iteration
for n = 1:iter;
t = tr(n);
H = Hr(n);         % %total integrated gain factor of SOA
t=t-tend*1e-12/2;                 
tparam=t./tau0;
tparam1=(t+(h/2))./tau0;
tparam2=(t+h)./tau0;
 
Pin=sqrt(ein/(2*tau0))*(sech(-tparam+3.52))+sqrt(ein/(2*tau0))*(sech(-tparam+10.57))+sqrt(ein/(2*tau0))*(sech(-tparam+17.65))+sqrt(ein/(2*tau0))*(sech(-tparam+45.83));
Pin2=sqrt(ein/(2*tau0))*(sech(-tparam1+3.52))+sqrt(ein/(2*tau0))*(sech(-tparam1+10.57))+sqrt(ein/(2*tau0))*(sech(-tparam1+17.65))+sqrt(ein/(2*tau0))*(sech(-tparam1+45.83));
Pin3=sqrt(ein/(2*tau0))*(sech(-tparam2+3.52))+sqrt(ein/(2*tau0))*(sech(-tparam2+10.57))+sqrt(ein/(2*tau0))*(sech(-tparam2+17.65))+sqrt(ein/(2*tau0))*(sech(-tparam2+45.83));
 
k1 =(h0-H)/tauc-(abs(Pin).^2).*(exp(H)-1)/esat;
k2 = (h0-(H+0.5*h*k1))/tauc-(abs(Pin2).^2).*(exp(H+0.5*h*k1)-1)/esat;
k3 = (h0-(H+0.5*h*k2))/tauc-(abs(Pin2).^2).*(exp(H+0.5*h*k2)-1)/esat;
k4 = (h0-(H+h*k3))/tauc-(abs(Pin3).^2).*(exp(H+h*k3)-1)/esat;
Hr(n+1) = Hr(n) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
tr(n+1) = tr(n) + h;
 
pin(n)=Pin;
tf(n)=t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain=exp(Hr);           %gain of SOA
phi=-0.5*alfa*Hr;       %phase difference os SOA
 
Pout=(abs(pin).^2).*gain(1:n);      %output power of SOA
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
subplot(3,1,1)
plot(tf/1e-12,(abs(pin).^2)*1000),grid
ylabel('Input power(mW)')
xlabel('(a) time(ps)')
ylim([0 3])
subplot(3,1,2)
plot(tf/1e-12,Hr(1:n)),grid
ylabel('h(t)')
xlabel('(b) time(ps)')
%ylim([4 7])
subplot(3,1,3)
plot(tf/1e-12,(Pout*1000),'b'),grid
ylabel('output power(mW)')
xlabel('(c) time(ps)')
ylim([0 89])
figure(2)
subplot(3,1,1)
plot(tf/1e-12,(abs(pin).^2)*1000),grid
ylabel('Input power(W)')
xlabel('time(ps)')
subplot(3,1,2)
plot(tf/1e-12,gain(1:n)),grid
ylabel('gain')
xlabel('time(ps)')
subplot(3,1,3)
plot(tf/1e-12,phi(1:n)),grid
ylabel('phi')
xlabel('time(ps)')

figure(3)
plot(tf/tau0,phi(1:n)/pi),grid
ylabel('phase/pi')
xlabel('time')
% Create an eye diagram scope object
h = commscope.eyediagram('PlotType', '2D Line', ...
                         'NumberOfStoredTraces', 100);
                    
update(h,Pout);
plot(h)
eyescope
 
 
 
