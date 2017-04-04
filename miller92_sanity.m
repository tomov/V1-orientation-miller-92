R = 5.5;
G = @(x, r) exp(- norm(x) / (r * R)^2);
kronecker = @(x) real(norm(x) == 0);

%{
[X,Y] = meshgrid(-20:.5:20);
Z = zeros(size(X));
for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        x = [X(i, j), Y(i, j)];
        Z(i, j) = G(x, r1);
    end
end
figure;
surf(Z,'EdgeColor','None');
%}


a = 1/2;
r1 = 0.4;
k = 0; % 'excit': 0, 'excit/inhibt': 1/9
x1 = 2.5; % 'excit': 2.5, 'excit/inhibt': 7.5
I = @(x) (a + (1 - a) * kronecker(x)) .* (G(x, r1) - k * G(x, 3 * r1)) .* (norm(x) <= x1);

%{
[X,Y] = meshgrid(-x1:.1:x1);
Z = zeros(size(X));
for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        x = [X(i, j), Y(i, j)];
        Z(i, j) = I(x);
    end
end
figure;
surf(Z,'EdgeColor','None');
%}


rc = 0.28; % 'excit': 0.28, 'excit/inhib': 0.2
C_ON_ON = @(x) G(x, rc) - (1/9) * G(x, 3*rc);
C_ON_OFF = @(x) -0.5; % ???

%{
[X,Y] = meshgrid(-x1 * 2 : .1 : x1 * 2);
Z = zeros(size(X));
for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        x = [X(i, j), Y(i, j)];
        Z(i, j) = C_ON_ON(x);
    end
end
figure;
surf(Z,'EdgeColor','None');
%}

% from http://mathworld.wolfram.com/Circle-CircleIntersection.html
A = @(x) circle_intersect(norm(x), 5, 2.5) * (norm(x) <= 5.5);

%{
[X,Y] = meshgrid(-6:.1:6);
Z = zeros(size(X));
for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        x = [X(i, j), Y(i, j)];
        Z(i, j) = A(x);
        assert(isreal(Z(i, j)));
    end
end
figure;
surf(Z,'EdgeColor','None');
%}

lambda_dt = 0.0012; % 'excit': 0.0012; 'excit/inhib': 0.0019
