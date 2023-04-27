function [tot_Received_P,Rice_factor,delay_spread] = powercalculate(transmit_pos,receive_pos)
parameters
environment
% Number of possible path combinations (excluding LOS):
path_combination_number = 2^object_number-1;
% Binary vector representing path combinations:
path_combinations = de2bi(0:path_combination_number,object_number);

% Determine maximum number of rays to cast:
% ray_collisions = maximum times of reflection for single ray 
ray_number = 0;
for i = 1:ray_collisions
    if i <= object_number
        ray_number = ray_number + factorial(object_number)/factorial(object_number-i);
    end
end
% Account for line-of-sight ray:
ray_number = ray_number + 1;

%% Initialise path information matrix:
path_matrix = zeros(ray_number,6);
path_matrix(:,5) = 1;
Ground_reflection=0;
diffrection_E=0;
LOS_distance=[];
Ground_reflection_dis=[];
diffrection_distance=[];
max_old_path_length=0;
min_old_path_length=1000;
% Initialise valid path counter:
m = 1;
% Switch back to ray tracing plot:
% figure(figure_counter)

% For each path combination:
for i = 1:path_combination_number+1
    % Get a path combination:
    path_combination_current = find(path_combinations(i,:)~=0);
    % i=1   path_combinations(i,:)=0 0 0
    %path_combination_current =  emepty
    % As long as path is below maximum collisions numbers:
    if length(path_combination_current) <= ray_collisions
        %maximum length(path_combination_current) = object number
        % All path permutations of particular path combination:
        path_permutations = unique(perms(path_combination_current),'rows');
        % Get number of permutations of particular path combination:
        [path_permutation_number,~] = size(path_permutations);
        % For each path permutation of particular path combination:
        for j = 1:path_permutation_number
            % Get a path permutation:
            path = path_permutations(j,:);
            % Number of collisions along path:
            path_collisions = length(path);
            % if path is empty , path_collisions=0  direclty jump to
            % Line-of-sight   (156)
            % Initialise mirrored points matrix:
            mirror_matrix = zeros(path_collisions+1,2);
            % Set first point to receiver position:
            mirror_matrix(1,:) = receive_pos;
            %% METHOD OF IMAGES %%%
            % Initialise validity flag:
            flag_ray = 0;
            % For each object in path:
            for k = 1:path_collisions
                % Get point to mirror:
                x = mirror_matrix(k,1);
                y = mirror_matrix(k,2);
                % Get object line :
                x1 = object_geometry(path(end+1-k),1);
                y1 = object_geometry(path(end+1-k),2);
                x2 = object_geometry(path(end+1-k),3);
                y2 = object_geometry(path(end+1-k),4);
                % Perpendicular distance from point to line:
                point_difference  = point_to_line_distance(x,y,x1,y1,x2,y2);
                % Mirrored point:

                    if y1-y2==0
                         mirror_matrix(k+1,1) = x;
                         mirror_matrix(k+1,2) = y+2*point_difference;
                    else                
                         mirror_matrix(k+1,1) = x+2*point_difference;
                         mirror_matrix(k+1,2) = y;
                    end

            end
            %% PATH POINTS %%%
            % Initialise path points matrix:
            points_matrix = zeros(path_collisions+2,2);
            % Set origin to transmitter position:
            points_matrix(1,:) = transmit_pos;
            % Set termination to receiver position:
            points_matrix(end,:) = receive_pos;

            %%% For each intersection point along the path:
            for k = 1:path_collisions
                % Get ray segment origin:
                A = points_matrix(k,:);
                % Get ray segment termination:
                B = mirror_matrix(end+1-k,:);
                [doesInt,inter]=checkSegmentIntersection([A ; B],[object_geometry(path(k),1),object_geometry(path(k),2);object_geometry(path(k),3),object_geometry(path(k),4)]);
                % If ray is parallel or collinear to wall:
                if doesInt==false 
                    % Set invalidity flag:
                    flag_ray = 1;
                    break
                end
                % Intersection point with object (not considering bounds):
                points_matrix(k+1,:) = [inter(1),inter(2)];
            end
                %%% Cross building  CHECK %%%
              %% For each building:
             for k = 1:path_collisions+1
                    A = points_matrix(k,:);
                    % Get termination of ray segment:
                    B = points_matrix(k+1,:);
                    sharkpoint=[0,0];
                    cross_building=0;
                    [LOS,~]=size(points_matrix);
                    for obj = 1:building_number
                            line(1:2,1:2)=[building_geometry(obj,1),building_geometry(obj,2);building_geometry(obj,3),building_geometry(obj,2)];%up wall
                            line(3:4,1:2)=[building_geometry(obj,1),building_geometry(obj,4);building_geometry(obj,3),building_geometry(obj,4)];%down wall
                            line(5:6,1:2)=[building_geometry(obj,3),building_geometry(obj,2);building_geometry(obj,3),building_geometry(obj,4)];%right wall
                            line(7:8,1:2)=[building_geometry(obj,1),building_geometry(obj,2);building_geometry(obj,1),building_geometry(obj,4)];%left wall
                            cross_wall_count=1;
                            cross_wall=zeros(4,2);                            
                            for objj=0:3
                                [doesInt,~]=checkSegmentIntersection([A;B],line(2*objj+1:(objj+1)*2,:));
                                if doesInt
                                    cross_building=cross_building+1;
                                    if LOS==2
                                        cross_wall(2*cross_wall_count+1:(cross_wall_count+1)*2,:)=line(2*objj+1:(objj+1)*2,:);
                                        cross_wall_count=cross_wall_count-1;
                                    end
                                end
                                
                            end

                         if cross_building >1
                               flag_ray = 1;
                               if LOS==2
                                   [~,sharkpoint]=checkSegmentIntersection(cross_wall(1:2,:),cross_wall(3:4,:));
                                   if isnan(sharkpoint)
                                       if points_matrix(end,1)<60
                                         sharkpoint=[max(cross_wall(:,1)),min(cross_wall(:,2))];
                                       else
                                         sharkpoint=[min(cross_wall(:,1)),min(cross_wall(:,2))];
                                       end
                                   end
                               end                               
                         end 
                         
                         cross_building=0;
                         
                    end
                   if LOS==2 && sharkpoint(1)~=0 && sharkpoint(2)~=0
                       %% draw_diffraction_ray
                       S1=rayDistance(A,sharkpoint);
                       S2=rayDistance(sharkpoint,B);
                       distance= rayDistance(A,B);
                       diffrection_distance=distance;
                       delta_r=(S1+S2)-distance;
                       V=sqrt((4*frequency_carrier/light_speed)*delta_r);
                       E=sqrt(transmit_gain*receive_gain)*(exp((-1j*2*pi*frequency_carrier/light_speed)*distance)/distance);
                       diffrection_E=E*((1+1j)/2)*((0.5 - fresnelc(V)) - (0.5 - fresnels(V))*1j);       
                        if isnan(diffrection_E) 
                            diffrection_E=0;
                        end                       
                   end
            end
            % If ray path may be valid:
            if flag_ray == 0
                
                % For each ray segment:
                for k = 1:path_collisions+1 
                    % Get origin of ray segment:
                    A=points_matrix(k,:);
                    B=points_matrix(k+1,:);
                    %%% OBJECT BOUNDS CHECK %%%
                    % Line-of-sight path:
                   if isempty(path) 
                        distance= rayDistance(A,B);
                        LOS_distance = distance; %% LOS distance
                      %% ground reflection
                        [Ground_reflection,Ground_reflection_dis,transmit_gain_gnd,receive_gain_gnd]=GroundReflection(distance);
                    % Multi-path segments:
                   elseif k ~= path_collisions+1         
                       distance= rayDistance(A,B);
                       path_matrix(m,1) = path_matrix(m,1)+distance;
                       % Update attenuation due to collisions:
                      [~,reflection_fector] = reflection_factor(A,B,object_geometry(path(k),1:2));
                      path_matrix(m,5) = path_matrix(m,5)*reflection_fector;
                       
                    % Terminating at receiver segment:
                   elseif k == path_collisions+1 && isempty(path)==0
                            % Draw ray segment: 
    %                                 Update path distance:
                           distance= rayDistance(A,B);
                           path_matrix(m,1) = path_matrix(m,1)+distance;
                        if mod(length(path),2) == 1
                            path_matrix(m,6) = -1;
                        end

                    end
                end
             %% find out the maximum and minimum length rays               
                if max_old_path_length<path_matrix(m,1)
                   max_old_path_length=path_matrix(m,1);
                end
                
               if path_matrix(m,1) ~=0 
                if min_old_path_length > path_matrix(m,1) 
                   min_old_path_length= path_matrix(m,1);
                end 
               end                
                                     
            end
        end
        % Update valid path counter
        m = m + 1;
    end
end

%% Clean path information matrix:
path_matrix(~any(path_matrix(:,1),2),:) = [];
%% Update path information matrix: 
path_matrix(:,2) = path_matrix(:,1)/light_speed; %% T1= d1/c  T2=d2/c... delay
LOS_received_E=sqrt(transmit_gain*receive_gain)*(exp((-1j*2*pi*frequency_carrier/light_speed)*LOS_distance)/LOS_distance); %% LOS electric field
if isempty(LOS_distance)
    LOS_received_E=0;
    Ground_reflection=0;
end
%% calculate received power  (change E component format ： since the gain of ground reflection is different with building reflection )
path_matrix(:,6) = sqrt(transmit_gain*receive_gain).*(path_matrix(:,5).*(exp((-1j*2*pi*frequency_carrier/light_speed).*path_matrix(:,1))./path_matrix(:,1))); %% electric field component
tot_Received_P=(light_speed/(frequency_carrier*4*pi))^2*P_TX*abs(sum(path_matrix(:,6))+Ground_reflection+diffrection_E+LOS_received_E)^2;


%% change E back to real electric field  format  : E= sqrt(60*EIRP)*(E^(-jβd)/d)
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
Rice_factor=10*log10(LOS_RICE/MPC_RICE);

%% delay spread
    
max_T=max([LOS_distance/light_speed,Ground_reflection_dis/light_speed,(path_matrix(:,1)/light_speed)',diffrection_distance/light_speed]);
min_T=min([LOS_distance/light_speed,Ground_reflection_dis/light_speed,(path_matrix(:,1)/light_speed)',diffrection_distance/light_speed]);
delay_spread=(max_T-min_T);
if isempty(path_matrix)
   delay_spread=0;
end
end

