clc
clear 
addpath(genpath('data'));
load('matlab10.mat')
load('matlab11.mat')
load('matlab12.mat')
%% plot received power heat map
figure
colormap('hot');
h=imagesc(log_heatmap);
set(h,'alphadata',~isnan(log_heatmap))
background_color = [0.5 0.2 0.1]; % 如设置为灰色
set(gca,'color',background_color)
colorbar;
title('power heatmap (db)');
set(gca,'YDir','normal');
xlabel('X position')
ylabel('Y position')

%% SNR heatmap
K=1.379*10^(-23);
T=290;
B=200*10^6;
Noise_P=12+10*log10(K*T*B);
[row,col]=find(log_heatmap~=-inf);
receive_pos=zeros(length(row),2);
SNR_heatmap=log_heatmap;
for i=1:length(row)
    for j=1:length(col)
      SNR_heatmap(row(i),col(j))=log_heatmap(row(i),col(j))-Noise_P;
    end
end
SNR_heatmap(95,35)=inf;
figure
colormap('hot');
h=imagesc(SNR_heatmap);
set(h,'alphadata',~isnan(SNR_heatmap))
background_color = [0.5 0.2 0.1]; % 如设置为灰色
set(gca,'color',background_color)
colorbar;
title('SNR heatmap (db)');
set(gca,'YDir','normal');
xlabel('X position')
ylabel('Y position')
%% rice factor heatmap
figure
colormap('hot');
h=imagesc(Rice_factor_hp);
set(h,'alphadata',~isnan(Rice_factor_hp))
background_color = [0.5 0.2 0.1]; % 如设置为灰色
set(gca,'color',background_color)
colorbar;
title('rice factor heatmap (db)');
set(gca,'YDir','normal');
xlabel('X position')
ylabel('Y position')
%% delay spread heatmap
figure
colormap('hot');
delay_spread_hp(95,35)=nan;
h=imagesc(delay_spread_hp);
set(h,'alphadata',~isnan(delay_spread_hp))
background_color = [0.3 0.3 0.3]; % 如设置为灰色
set(gca,'color',background_color)
colorbar;
title('delay spread heatmap (S)');
set(gca,'YDir','normal');
xlabel('X position')
ylabel('Y position')
