
%% Update path information matrix:
path_matrix(:,2) = path_matrix(:,1)/light_speed; %% T1= d1/c  T2=d2/c... delay
LOS_received_E=sqrt(transmit_gain*receive_gain)*(exp((-1j*2*pi*frequency_carrier/light_speed)*LOS_distance)/LOS_distance); %% LOS electric field
if isempty(LOS_distance)
    LOS_received_E=0;
    Ground_reflection=0;
end
path_matrix(:,6) = sqrt(transmit_gain*receive_gain).*(path_matrix(:,5).*(exp((-1j*2*pi*frequency_carrier/light_speed).*path_matrix(:,1))./path_matrix(:,1))); %% electric field
%% total power
tot_Received_P=(light_speed/(frequency_carrier*4*pi))^2*P_TX*abs(sum(path_matrix(:,6))+Ground_reflection+diffrection_E+LOS_received_E)^2;

%% real electric field
path_matrix(:,6)=(path_matrix(:,6)./sqrt(transmit_gain*receive_gain)).*sqrt(60*EIRP);
if LOS_received_E~=0
    LOS_received_E=(LOS_received_E./sqrt(transmit_gain*receive_gain)).*sqrt(60*EIRP);
end
if Ground_reflection~=0
    Ground_reflection=(Ground_reflection./sqrt(transmit_gain_gnd*receive_gain_gnd)).*sqrt(60*EIRP);
end
if diffrection_E~=0
    diffrection_E=(diffrection_E./sqrt(transmit_gain*receive_gain)).*sqrt(60*EIRP);
end
%% MPC
LOS_RICE=abs(LOS_received_E)^2;
MPC_RICE=sum(abs(path_matrix(:,6)).^2)+abs(diffrection_E)^2+abs(Ground_reflection)^2;
%% Rice
Rice_factor=db(LOS_RICE/MPC_RICE);
%% delay spread
max_T=max([LOS_distance/light_speed,Ground_reflection_dis/light_speed,(path_matrix(:,1)/light_speed)',diffrection_distance/light_speed]);
min_T=min([LOS_distance/light_speed,Ground_reflection_dis/light_speed,(path_matrix(:,1)/light_speed)',diffrection_distance/light_speed]);
delay_spread=(max_T-min_T);
if isempty(path_matrix)
   delay_spread=0;
end

%% impulse response
if isempty(LOS_distance)
    LOS_received_E=[];
    Ground_reflection=[];
end
if isempty(diffrection_distance)
    diffrection_E=[];
end
E_vector=[LOS_received_E,LOS_distance/light_speed;Ground_reflection,Ground_reflection_dis/light_speed;diffrection_E,diffrection_distance/light_speed;path_matrix(:,6),path_matrix(:,2)];
E_vector(~any(E_vector(:,1),2),:) = [];

N= 1000;
Ts = max_T/N;
fs = 1/Ts;
T=Ts*(1:N);%% time vector
impulse_response=zeros(1,length(T));
delay_vector = E_vector(:,2);
for i=1:length(T)
    impulse_response(i)=impulse_calculate(T(i),E_vector(:,1),delay_vector);
end
T1=T(find(impulse_response~=0));
impulse_response(impulse_response==0)=[];
figure 
hold on
title('physical impulse response')
xlabel('T')
ylabel('|h_t|')
stem(T1,impulse_response)
hold off
figure
f1=fs/length(impulse_response)*(1:length(impulse_response));
plot(f1,db(fft(impulse_response)))
xlabel('f Hz')
ylabel('|h_f| db')

%% TDL Impulse response
bandwidth=[1*10^6,50*10^6,100*10^6,200*10^6];
TDL_response=zeros(length(bandwidth),length(T));
for j=1:length(bandwidth)
for i=1:length(T)
    TDL_response(j,i)=TDL_calc(T(i),E_vector(:,1),delay_vector,bandwidth(j));
end
end

TDL_response_2mhz=TDL_response(1,:);
TDL_response_2mhz(TDL_response_2mhz==0)=nan;
TDL_response_50mhz=TDL_response(2,:);
TDL_response_50mhz(TDL_response_50mhz==0)=nan;
TDL_response_100mhz=TDL_response(3,:);
TDL_response_100mhz(TDL_response_100mhz==0)=nan;
TDL_response_200mhz=TDL_response(4,:);
TDL_response_200mhz(TDL_response_200mhz==0)=nan;
figure 
subplot(2,2,1)
hold on
title('TDL impulse response bandwidth=2MHZ Δt=250ns')
xlabel('T')
ylabel('|TDL|')
stem(T,TDL_response_2mhz)
hold off
subplot(2,2,2)
hold on
title('TDL impulse response bandwidth=50MHZ Δt=10ns')
xlabel('T')
ylabel('|TDL|')
stem(T,TDL_response_50mhz)
hold off
subplot(2,2,3)
hold on
title('TDL impulse response bandwidth=100MHZ Δt=5ns')
xlabel('T')
ylabel('|TDL|')
stem(T,TDL_response_100mhz)
hold off
subplot(2,2,4)
hold on
title('TDL impulse response bandwidth=200MHZ Δt=2.5ns')
xlabel('T')
ylabel('|TDL|')
stem(T,TDL_response_200mhz)
hold off
%% Uncorrelated scattering TDL impulse response 
% bandwidth=[2*10^6,50*10^6,100*10^6,200*10^6];
% TDL_response=zeros(length(bandwidth),length(T));
% for j=1:length(bandwidth)
% for i=1:length(T)
%     TDL_response(j,i)=US_TDL_calc(T(i),E_vector(:,1),delay_fs,bandwidth(j));
% end
% end
% 
% TDL_response_2mhz=TDL_response(1,:);
% TDL_response_2mhz(TDL_response_2mhz==0)=nan;
% TDL_response_50mhz=TDL_response(2,:);
% TDL_response_50mhz(TDL_response_50mhz==0)=nan;
% TDL_response_100mhz=TDL_response(3,:);
% TDL_response_100mhz(TDL_response_100mhz==0)=nan;
% TDL_response_200mhz=TDL_response(4,:);
% TDL_response_200mhz(TDL_response_200mhz==0)=nan;
% 
% figure 
% subplot(2,2,1)
% hold on
% title('US TDL impulse response bandwidth=2MHZ Δt=250ns')
% xlabel('T')
% ylabel('|TDL|')
% stem(T,TDL_response_2mhz)
% hold off
% subplot(2,2,2)
% hold on
% title('US TDL impulse response bandwidth=50MHZ Δt=10ns')
% xlabel('T')
% ylabel('|TDL|')
% stem(T,TDL_response_50mhz)
% hold off
% subplot(2,2,3)
% hold on
% title('US TDL impulse response bandwidth=100MHZ Δt=5ns')
% xlabel('T')
% ylabel('|TDL|')
% stem(T,TDL_response_100mhz)
% hold off
% subplot(2,2,4)
% hold on
% title('US TDL impulse response bandwidth=200MHZ Δt=2.5ns')
% xlabel('T')
% ylabel('|TDL|')
% stem(T,TDL_response_200mhz)
% hold off
