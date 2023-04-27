function [doesIntersect, intersection] = checkSegmentIntersection(A, B)

% Check if two line segments (A and B) intersect in 2D space.
% 
% Usage: [doesIntersect, intersection] = checkSegmentIntersection(A, B, plotResults)
% 
% Inputs:
%  A and B are setup with two rows and 2 columns for the 2 points which
%  define the line segments and 2 columns for 2 dimensions
% 
% Outputs:
%  doesIntersect - provides a boolean about whether the segments
%    intersect.
%  intersection - The second output argument provides the point at which the intersection
%    occurs.

% initialize
doesIntersect = false;
intersection = NaN;

% Error check
if any(size(A) ~=2) || any(size(B) ~=2)
    error('Line segments are not properly defined.');
end

% Given two line segments:
% segment1 = p + t*r
% segment2 = q + u*s
% Can solve intersection of two segments and values of u and t where
% intersection occurs (where x is cross product)
% t = (q-p) x s/(r x s)
% u = (q-p) x r/(r x s)

% Solve for all of these variables given the two line segments
p = A(1,:);
r = A(2,:)-A(1,:);
q = B(1,:);
s = B(2,:)-B(1,:);
% Solve for cross products
r_cross_s = cross(r, s);
q_p_cross_s = cross(q-p, s);
q_p_cross_r = cross(q-p, r);
% solve for t and u
t = q_p_cross_s / r_cross_s;
u = q_p_cross_r / r_cross_s;

% Five possibilities
% 1) Collinear and overlapping
% 2) Collinear and disjoint
% 3) Parallel and non-intersecting
% 4) Not parallel and intersect
% 5) Not parallel and non-intersecting

%% First Possibility
if r_cross_s == 0
    if q_p_cross_r == 0
        t0 = dot(q-p,r)/dot(r,r);
        if t0 >= 0 && t0 <= 1
            doesIntersect = false;
            % return a line segment where intersection occurs
%             intersection = [p; p+t0*r];
           intersection = NaN;
        end
    end
%% Fourth possibility
else
    if t >= 0 && t <= 1 && u >=0 && u <= 1
        doesIntersect = true;
        intersection = p + t * r;
    end    
end
%% Second, third, and fifth possibilities return initialized values

% %% plot and check results
% if plotResults
%     figure; hold on;
%     plot(A(:,1), A(:,2), 'ro-')
%     plot(B(:,1), B(:,2), 'bo-')
%     if doesIntersect
%         if numel(intersection(:,1)) == 2
%             plot(intersection(:,1), intersection(:,2), 'kx-', 'lineWidth', 2)
%         else
%             plot(intersection(1), intersection(2), 'ks')
%         end
%         legend('A', 'B', 'Intersection');
%     else
%         legend('A', 'B');
%     end
% end

%% Cross product returning z value
function z = cross(a, b)
    z = a(1)*b(2) - a(2)*b(1);