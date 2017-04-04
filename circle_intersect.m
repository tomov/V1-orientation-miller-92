function area = circle_intersect( d, r, R )
% area of intersection of two circles of radii r and R, separated by d
% from http://mathworld.wolfram.com/Circle-CircleIntersection.html
%
circle_helper = @(d, R) R^2 * acos(d / R) - d .* sqrt(R^2 - d.^2);

area = zeros(size(d));
too_close = d <= min(r, R);
too_far = d >= r + R;
just_right = ~too_close & ~too_far;

area(too_close) = pi * min(r, R)^2;
area(too_far) = 0;
area(just_right) = circle_helper((d(just_right).^2 - r^2 + R^2) ./ (2 * d(just_right)), R) - ...
                   circle_helper((d(just_right).^2 - R^2 + r^2) ./ (2 * d(just_right)), r);


