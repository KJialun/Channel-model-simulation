figure_counter = 1;
figure(figure_counter)
hold on
xlabel('x')
ylabel('y')
plot(transmit_pos(1),transmit_pos(2),'x','MarkerSize',8,'LineWidth',3,'Color',blue)
plot(receive_pos(1),receive_pos(2),'x','MarkerSize',8,'LineWidth',3,'Color',red)
%% draw wall
% for i = 1:object_number
%     
%   draw_wall(object_geometry(i,1),object_geometry(i,2),object_geometry(i,3),object_geometry(i,4))
% 
% end

%% draw building:
for i = 1:building_number
    x1=building_geometry(i,1);
    y1=building_geometry(i,2);
    x2=building_geometry(i,3);
    y2=building_geometry(i,4);
    x=[x1,x2,x2,x1];
    y=[y1,y1,y2,y2];
    fill(x,y,'k')
end


%% ray_tracing

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
max_old_path_point=[];

double_reflection=[];
single_reflection=[];

% Initialise valid path counter:
m = 1;
num_single=0;
num_double=0;
% kkk=0;
% iii=0;
% Switch back to ray tracing plot:
figure_counter=1;
figure(figure_counter)

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
                            line(1:2,1:2)=[building_geometry(obj,1),building_geometry(obj,2);building_geometry(obj,3),building_geometry(obj,2)];
                            line(3:4,1:2)=[building_geometry(obj,1),building_geometry(obj,4);building_geometry(obj,3),building_geometry(obj,4)];
                            line(5:6,1:2)=[building_geometry(obj,3),building_geometry(obj,2);building_geometry(obj,3),building_geometry(obj,4)];
                            line(7:8,1:2)=[building_geometry(obj,1),building_geometry(obj,2);building_geometry(obj,1),building_geometry(obj,4)];
                            
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
                       draw_diffraction_ray(A,sharkpoint,ray_resolution)
                       draw_diffraction_ray(sharkpoint,B,ray_resolution)
                       S1=rayDistance(A,sharkpoint);
                       S2=rayDistance(sharkpoint,B);
                       distance= rayDistance(A,B);
                       diffrection_distance=distance;
                       delta_r=(S1+S2)-distance;
                       V=sqrt((4*frequency_carrier/light_speed)*delta_r);
                       E=sqrt(transmit_gain*receive_gain)*(exp((-1j*2*pi*frequency_carrier/light_speed)*distance)/distance);
                       diffrection_E=E*((1+1j)/2)*((0.5 - fresnelc(V)) - (0.5 - fresnels(V))*1j);                  
                   end
            end
            % If ray path may be valid:
            if flag_ray == 0
                
                if path_collisions==1
                    num_single=num_single+1;
                    single_reflection(:,:,num_single)=points_matrix;
                end
                if path_collisions==2
                     num_double=num_double+1;
                    double_reflection(:,:,num_double)=points_matrix;                  
                end                                
                % For each ray segment:
                for k = 1:path_collisions+1 
                    % Get origin of ray segment:
                    A=points_matrix(k,:);
                    B=points_matrix(k+1,:);
                    reflection_fector=1;
                    %%% OBJECT BOUNDS CHECK %%%
                    % Line-of-sight path:
                   if isempty(path) 
                        distance= rayDistance(A,B);
                        LOS_distance = distance; %% LOS distance
                      %% ground reflection
                        [Ground_reflection,Ground_reflection_dis,transmit_gain_gnd,receive_gain_gnd]=GroundReflection(distance);
                         Gound_reflection_ray=[A;B];
                    % Multi-path segments:
                   elseif k ~= path_collisions+1                             
                       distance= rayDistance(A,B);
                       path_matrix(m,1) = path_matrix(m,1)+distance;
                       % Update attenuation due to collisions:
                      [Theta,reflection_fector] = reflection_factor(A,B,object_geometry(path(k),1:2));
                      path_matrix(m,5) = path_matrix(m,5)*reflection_fector;
                       
                    % Terminating at receiver segment:
                   elseif k == path_collisions+1 && isempty(path)==0
   
                            % Update path distance:
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
                   max_old_path_point=points_matrix;
                end
                
               if path_matrix(m,1) ~=0 
                if min_old_path_length > path_matrix(m,1) 
                   min_old_path_length= path_matrix(m,1);
                   min_old_path_point=points_matrix;
                end 
               end                
                       
            end
        end
        % Update valid path counter
        m = m + 1;
    end
end


 %% draw LOS and ground reflection ray
if ~isempty(Gound_reflection_ray)
 draw_ground_reflection(Gound_reflection_ray(1,:),Gound_reflection_ray(2,:),ray_resolution)   
end

%% draw single reflection ray

if ~isempty(single_reflection)
    for j=1:num_single
        single_reflection1=single_reflection(:,:,j);
        for i=1:(length(single_reflection1)/2+1)
            drawRay(single_reflection1(i,:),single_reflection1(i+1,:),ray_resolution)   
        end
    end
end

%% draw double reflection ray
if ~isempty(double_reflection)
  for j=1:num_double
     double_reflection1=double_reflection(:,:,j);
    for i=1:(length(double_reflection1)/2+1)
        draw_double_reflection_ray(double_reflection1(i,:),double_reflection1(i+1,:),ray_resolution)   
    end
  end
end
%% draw maximum length ray
if ~isempty(max_old_path_point)
for i=1:(length(max_old_path_point)/2+1)
    draw_maxlength_ray(max_old_path_point(i,:),max_old_path_point(i+1,:),ray_resolution)   
end
end
legend('Tx',' Rx')
%% Clean path information matrix:
path_matrix(~any(path_matrix(:,1),2),:) = [];