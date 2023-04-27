%% Geometric properties of building:
close all
clear
clc
addpath(genpath('function'));

building_geometry(1,1:4) = [[0,130],[70,100]];
building_geometry(2,1:4) = [[0,90],[30,35]];
building_geometry(3,1:4) = [[0,25],[30,10]];
building_geometry(4,1:4) = [[80,130],[110,75]];
building_geometry(5,1:4) = [[80,65],[110,35]];
building_geometry(6,1:4) = [[80,25],[110,10]];

[building_number,~] = size(building_geometry);

transmit_pos = [35 , 95];

x = 1:110;
y = 1:130;
[X,Y] = meshgrid(x,y);
heat_map=zeros(size(X))*nan;

for i=1:building_number
    x1=building_geometry(i,1)+1;
    y1=building_geometry(i,2);
    x2=building_geometry(i,3);
    y2=building_geometry(i,4);
    for ii=y2:y1
        for jj=x1:x2
            X(ii,jj)=0;
        end
    end 
end

% r=10;
% indicator = (X-35).^2 + (Y-95).^2 < r^2 ;
% X(indicator)=0;
[row,col]=find(X~=0);
receive_pos=zeros(length(row),2);
for i=1:length(row)
   receive_pos(i,1)= X(row(i),col(i));
   receive_pos(i,2)= Y(row(i),col(i));
end
power=zeros(1,length(row));
Rice_factor=zeros(1,length(row));
delay_spread=zeros(1,length(row));
for i=1:length(row)
    [power(1,i),Rice_factor(1,i),delay_spread(1,i)]=powercalculate(transmit_pos,receive_pos(i,:));
end
Rice_factor_hp=heat_map;
delay_spread_hp=heat_map;
for i=1:length(row)
    heat_map(row(i),col(i))=power(1,i);
    Rice_factor_hp(row(i),col(i))=Rice_factor(1,i);
    delay_spread_hp(row(i),col(i))=delay_spread(1,i);
end

heat_map(95,35)=inf;
Rice_factor_hp(95,35)=inf;
delay_spread_hp(95,35)=-inf;
log_heatmap=10*log10(heat_map);
