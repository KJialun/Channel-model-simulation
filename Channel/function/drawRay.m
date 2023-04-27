function drawRay(origin,termination,step_size)

% Name:         drawRay
% Description:  Draws a ray using start and end point and determines ray
%               distance
% Inputs:       origin: start point of ray
%               termination: end point of ray
%               step_size: number of points to draw line
% Outputs:      distance: distance from start and end point of ray

direction = termination-origin;
distance = sqrt(direction(1)^2+direction(2)^2);
direction = direction/distance;

for i = 0:step_size:distance
    ray_point = origin+i*direction;
    plot(ray_point(1),ray_point(2),'.g','MarkerSize',2)
end

end