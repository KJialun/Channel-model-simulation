function draw_ground_reflection(origin,termination,step_size)
direction = termination-origin;
distance = sqrt(direction(1)^2+direction(2)^2);
direction = direction/distance;

for i = 0:step_size:distance
    ray_point = origin+i*direction;
    plot(ray_point(1),ray_point(2),'.b','MarkerSize',2)
end

end

