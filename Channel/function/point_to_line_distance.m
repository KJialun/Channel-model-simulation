function d = point_to_line_distance(point_x,point_y,x1,y1,x2,y2)
if (point_x==x1) && (point_x==x2)
d=0;
return
end
if (point_y==y1) && (point_y==y2)
d=0;
return
end
if x1-x2==0 
d=x1-point_x;
return
end
if y1-y2==0
d=y1-point_y;
return
end

% line= polyfit(x,y,1);
% a=(y2-y1)/(x2-x1);
% c=y2-(a*x2);
% b=-1;
% d=(a*point_x+b*point_y+c)/sqrt(a*a+b*b);
end
