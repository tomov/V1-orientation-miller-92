function area = circle_intersect( d, r, R )
% area of intersection of two circles of radii r and R, separated by d
% from http://mathworld.wolfram.com/Circle-CircleIntersection.html
%

circle_helper = @(d, R) R^2 * acos(d / R) - d * sqrt(R^2 - d^2);
if d <= min(r, R)
    area = pi * min(r, R)^2;
else
    % if d > min(r, R)
    area = circle_helper((d^2 - r^2 + R^2) / (2*d), R) - circle_helper((d^2 - R^2 + r^2) / (2*d), r);
end



