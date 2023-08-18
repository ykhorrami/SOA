%%%%%%%%%%%%%%%%%%%%%%%%%%% In The Name of Allah %%%%%%%%%%%%%%%%%%%%%%%%%%
%{
                                 SOA Simulator
            Fourth-order Runge-kutta finite difference simulator for 
          Semiconductor Optical Amplifier (SOA) for 1110001 bit stream
                        %% Written by Yaser Khorrami %%
                            % Y.khorrami@gmail.com %
                                 %% version 1.1 Jan 2022 %%
%}

clc;
clear all;
close all;

%% Initialize The Parameters
t0   = 0.05e-12;          % initial time (starting time)
H0   = 0;                 % initial condition (H(t0))
h    = 0.05e-12;          %simulation steps
tp   = 5e-12;             %FWHM
tau0 = tp/1.7627;
tauc = 200e-12;           %carrier recovery time
alfa = 5;                 %linewidth enhancement factor
esat = 10e-12;   
ein  = 15e-15;            %input energy
L    = 6.38336677e-3;     %length of SOA
g0   = 1096.6/2;          %linear gain(1/m)
h0   = g0*L;              %unsaturated gain
tend = 40;
bitnum = 7;
%% Start the iteration
Hr(1) = h0;
tr(1) = t0;
iter = 3200;            %iteration
for n = 1:iter
t = tr(n);
H = Hr(n);              % %total integrated gain factor of SOA
t = t-tend*1e-12/2;                 
tparam = t./tau0;
tparam1 = (t+(h/2))./tau0;
tparam2 = (t+h)./tau0;
 
Pin   = input(tparam,tau0,ein);
Pin1  = input(tparam1,tau0,ein);
Pin2  = input(tparam2,tau0,ein);

k1 = (h0-H)/tauc-(abs(Pin).^2).*(exp(H)-1)/esat;
k2 = (h0-(H+0.5*h*k1))/tauc-(abs(Pin1).^2).*(exp(H+0.5*h*k1)-1)/esat;
k3 = (h0-(H+0.5*h*k2))/tauc-(abs(Pin1).^2).*(exp(H+0.5*h*k2)-1)/esat;
k4 = (h0-(H+h*k3))/tauc-(abs(Pin2).^2).*(exp(H+h*k3)-1)/esat;

Hr(n+1) = Hr(n) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
tr(n+1) = tr(n) + h;
 
pin(n) = Pin;
tf(n) = t;
end
%% Gain Calculation
gain = exp(Hr);                   %Gain of SOA
phi = -0.5*alfa*Hr;               %Phase difference os SOA
 
Pout = (abs(pin).^2).*gain(1:n);  %Output power of SOA
%% Plot Outputs
figure(1)
subplot(3,1,1)
plot(tf/1e-12,(abs(pin).^2)*1000,'linewidth',1.5),grid
ylabel('P_{in}(mW)'); xlabel('Time(ps)'); title('(a)')
ylim([0 3])
set(gca,'linewidth', 1.5,'fontsize',12,'fontname','Times New Roman')
subplot(3,1,2)
plot(tf/1e-12,Hr(1:n),'linewidth',1.5),grid
ylabel('h(t)'); xlabel('Time(ps)'); title('(b)')
set(gca,'linewidth', 1.5,'fontsize',12,'fontname','Times New Roman')
subplot(3,1,3)
plot(tf/1e-12,(Pout*1000),'b','linewidth',1.5),grid
ylabel('P_{out}(mW)'); xlabel('Time(ps)'); title('(c)')
ylim([0 89])
set(gca,'linewidth', 1.5,'fontsize',12,'fontname','Times New Roman')

figure(2)
subplot(3,1,1)
plot(tf/1e-12,(abs(pin).^2)*1000,'linewidth',1.5),grid
ylabel('P_{in}(W)'); xlabel('Time(ps)')
set(gca,'linewidth', 1.5,'fontsize',12,'fontname','Times New Roman')
subplot(3,1,2)
plot(tf/1e-12,gain(1:n),'linewidth',1.5),grid
ylabel('Gain'); xlabel('Time(ps)')
set(gca,'linewidth', 1.5,'fontsize',12,'fontname','Times New Roman')
subplot(3,1,3)
plot(tf/1e-12,phi(1:n),'linewidth',1.5),grid
ylabel('\phi'); xlabel('Time(ps)')
set(gca,'linewidth', 1.5,'fontsize',12,'fontname','Times New Roman')

figure(3)
plot(tf/tau0,phi(1:n)/pi,'linewidth',1.5),grid
ylabel('\phi/\pi'); xlabel('t_{f}/\tau_{0}')
set(gca,'linewidth', 1.5,'fontsize',12,'fontname','Times New Roman')
 
 
 
