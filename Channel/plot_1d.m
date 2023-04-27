clc
clear 
addpath(genpath('function'));
%% plot received power 1D plot
figure
subplot(2,3,1);
transmit_pos=[35,95];
minimumd=10;
k=(0-95)/(80-35); % eastern corner （80，0）
b=95-k*35;
syms x
equ=(x-35)^2+(k*x+b-95)^2;
x1=double(solve(equ==minimumd)); %%  Minimal distance to the base station: 10 m 
x=x1(2):0.5:80;
y=k*x+b;

receive_pos=[x',y'];

for i=1:length(receive_pos)
[POWER_ONED_d0(i),Rice_factor(i),delay_spread(i)]= powercalculate(transmit_pos,receive_pos(i,:));
end
POWER_ONED_d0=10*log10(POWER_ONED_d0);
distance=log10((sqrt(45^2+95^2)/length(POWER_ONED_d0))*(1:length(POWER_ONED_d0)));
plot(distance,POWER_ONED_d0);
title('Received Power (dbw)');
xlabel('distance  log d')
ylabel('Received Power db')
hold on
%% linear regression  average power
p = polyfit(distance,POWER_ONED_d0,1);
average_power=polyval(p,distance);
plot(distance,average_power,'r');
hold off

%% SNR PATHLOSS 1d plot
subplot(2,3,6);
K=1.379*10^(-23);
T=290;
B=200*10^6;
Noise_P=12+10*log10(K*T*B);

plot(distance,POWER_ONED_d0-Noise_P);
title('SNR 1d (db)');
xlabel('distance log d')
ylabel('SNR db')
%% rice factor 1d plot
subplot(2,3,4);
plot(distance,Rice_factor);
title('rice factor 1d (db)');
xlabel('distance log d')
ylabel('rice factor in db')
%% delay spread 1d plot
subplot(2,3,5);
plot(distance,delay_spread);
title('delay spread 1d (S)');
xlabel('distance log d')
ylabel('delay (S)')


%% fading variability 
subplot(2,3,2);
hold on
p_p_av= POWER_ONED_d0-average_power;
plot(distance,p_p_av);
plot(distance,average_power-average_power,'r');

%% add fade margin 
sigma=sqrt(var(p_p_av));
n=-1*p(1)/10;%path loss exponent
connection_p=[0.5,0.7,0.9,0.95,0.99];  %% probability connection
for i=1:length(connection_p)
fade_margin(i)=erfcinv((1-connection_p(i))*2)*sigma*sqrt(2);
end
plot(distance,average_power-average_power-fade_margin(4),'k');
title('p - p_av (db) connection probability = 0.95');
xlabel('distance log d')
ylabel('Received Power db')
legend('Received Power','average power','sensitivity')
hold off

%% Cell range as a function of connection probability at cell edge
subplot(2,3,3)
plot(distance,POWER_ONED_d0);
hold on
distance1=log10((sqrt(45^2+95^2)/length(POWER_ONED_d0))*(1:length(POWER_ONED_d0)+30));
SNR_TARGET=5;  %5db
p_sensitivity= SNR_TARGET+Noise_P;
p_sensitivity=p_sensitivity*ones(1,length(distance1));
average_power1=polyval(p,distance1);
plot(distance1,p_sensitivity);
plot(distance1,average_power1,'k');
plot(distance1,p_sensitivity+fade_margin(4),'b');
title(' connection probability = 0.95');
xlabel('distance log d')
ylabel('Received Power db')
legend('Received Power','sensitivity','average power','fade margin')
hold off

%% path loss model 
transmit_gain = (16/(3*pi))*(sin(pi/2))^3;             
EIRP=0.25;
P_TX=10*log10(EIRP/transmit_gain);

figure
hold on
L=P_TX-POWER_ONED_d0;
plot(distance,L)
d0=10;
L0_d0=P_TX-POWER_ONED_d0(1);
L0_d= L0_d0 + 10*n*(distance);
plot(distance,L0_d,'r')
title(' path loss model ');
xlabel('distance log d')
ylabel('Received Power db')
hold off
%% Cell range as a function of connection probability at cell edge
L_MAX=P_TX-(SNR_TARGET+Noise_P);
connection_p=0.5:0.01:0.99;  %% probability connection
p_sensitivity= SNR_TARGET+Noise_P;
for i=1:length(connection_p)
fade_margin(i)=erfcinv((1-connection_p(i))*2)*sigma*sqrt(2);
end
%% prx_d=prx_d0-10*n*log(d)  logd= (prx_d-prx_d0)/10*n
for i=1:length(connection_p)
Cell_range(i)=((p_sensitivity+fade_margin(i) - p(2))/p(1));
end
figure
plot(connection_p,Cell_range,'r')
title(' Cell range as a function of p ');
xlabel('connection p')
ylabel('Cell_range log d')
%% connection probability therough the whole cell
fade_margin_v=1:20;
for i=1:length(fade_margin_v)
a=fade_margin_v(i)/((sigma*sqrt(2)));
b=(1/(sigma*sqrt(2)))*10*n*log(exp(1));
Fu(i)=1-0.5*erfc(a)+0.5*exp(2*a/b + 1/(b^2))*erfc(a+1/b);
end
figure
plot(fade_margin_v,Fu)
title('connection probability therough the whole cell: sigma=3.5491');
xlabel('fade margin')
ylabel('Fu')
