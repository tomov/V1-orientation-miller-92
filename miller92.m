miller92_sanity; % include lambdas

maxX = 10; % dimension of grid

[X,Y] = meshgrid(1:1:maxX);
id_to_coords = [X(:) Y(:)]; % mapping of coordinates to index
coords_to_id = zeros(size(X)); % mapping of index to coordinates
for id = 1:size(id_to_coords, 1)
    coords_to_id(X(id), Y(id)) = id;
    assert(id_to_coords(coords_to_id(X(id), Y(id)), 1) == X(id));
    assert(id_to_coords(coords_to_id(X(id), Y(id)), 2) == Y(id));
end
N = numel(X); % number of neurons in each layer


S_ON = zeros(N, N);
S_OFF = zeros(N, N);

rng('default');

for x_id = 1:N
    x = id_to_coords(x_id, :);
    for alpha_id = 1:N
        alpha = id_to_coords(alpha_id, :);
        S_ON(x_id, alpha_id) = (rand * 0.4 + 0.8) * A(x - alpha);
        S_OFF(x_id, alpha_id) = (rand * 0.4 + 0.8) * A(x - alpha);
    end
end



[ds_on, ds_off] = LS([3 4], [2 5], S_ON, S_OFF, id_to_coords)

[ds_on, ds_off] = dS([3 4], [3 4], S_ON, S_OFF, id_to_coords)


%{
[X,Y] = meshgrid(1 : maxX);
Z_ON = zeros(size(X));
Z_OFF = zeros(size(X));
for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        x = [X(i, j), Y(i, j)];
        [ls_on, ls_off] = LS([3 4], [i j], S_ON, S_OFF, id_to_coords)
        Z_ON(i, j) = ls_on;
        Z_OFF(i, j) = ls_off;
        assert(isreal(Z_ON(i, j)));
        assert(isreal(Z_OFF(i, j)));
    end
end
figure;
surf(Z_OFF,'EdgeColor','None');
%}


[X,Y] = meshgrid(1 : maxX);
Z_ON = zeros(size(X));
Z_OFF = zeros(size(X));
for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        x = [X(i, j), Y(i, j)];
        [ds_on, ds_off] = dS([3 4], [i j], S_ON, S_OFF, id_to_coords)
        Z_ON(i, j) = ds_on;
        Z_OFF(i, j) = ds_off;
        assert(isreal(Z_ON(i, j)));
        assert(isreal(Z_OFF(i, j)));
    end
end
figure;
surf(Z_ON,'EdgeColor','None');



% technically LS * dt
%
function [LS_ON, LS_OFF] = LS(x, alpha, S_ON, S_OFF, id_to_coords)
    maxX = 10;
    N = maxX * maxX;
    
    R = 5.5;
    G = @(x, r) exp(- norm(x) / (r * R)^2);
    kronecker = @(x) real(norm(x) == 0);

    a = 1/2;
    r1 = 0.4;
    k = 0; % 'excit': 0, 'excit/inhibt': 1/9
    x1 = 2.5; % 'excit': 2.5, 'excit/inhibt': 7.5
    I = @(x) (a + (1 - a) * kronecker(x)) .* (G(x, r1) - k * G(x, 3 * r1)) .* (norm(x) <= x1);

    rc = 0.28; % 'excit': 0.28, 'excit/inhib': 0.2
    C_ON_ON = @(x) G(x, rc) - (1/9) * G(x, 3*rc);
    C_ON_OFF = @(x) -0.5; % ???
    C_OFF_OFF = C_ON_ON; % ???????????
    C_OFF_ON = C_ON_OFF; % ???????????

    A = @(x) circle_intersect(norm(x), 5, 2.5) * (norm(x) <= 5.5);

    lambda_dt = 0.0012; % 'excit': 0.0012; 'excit/inhib': 0.0019

    sum_ON = 0;
    sum_OFF = 0;
    for y_id = 1:N
        y = id_to_coords(y_id, :);
        for beta_id = 1:N
            beta = id_to_coords(beta_id, :);
            sum_ON = sum_ON + I(x - y) * (C_ON_ON(alpha - beta) * S_ON(y_id, beta_id) + ...
                                          C_ON_OFF(alpha - beta) * S_OFF(y_id, beta_id));
            sum_OFF = sum_OFF + I(x - y) * (C_OFF_OFF(alpha - beta) * S_OFF(y_id, beta_id) + ...
                                            C_OFF_ON(alpha - beta) * S_ON(y_id, beta_id));
        end
    end
    LS_ON = lambda_dt * A(x - alpha) * sum_ON;
    LS_OFF = lambda_dt * A(x - alpha) * sum_OFF;
end



function [dS_ON, dS_OFF] = dS(x, alpha, S_ON, S_OFF, id_to_coords)
    maxX = 10;
    N = maxX * maxX;
    
    A = @(x) circle_intersect(norm(x), 5, 2.5) * (norm(x) <= 5.5);
    
    sum_A = 0;
    for beta_id = 1:N
        beta = id_to_coords(beta_id, :);
        sum_A = sum_A + A(x - beta);
    end
    
    sum_LS = 0;
    for beta_id = 1:N
        beta = id_to_coords(beta_id, :);
        [LS_ON, LS_OFF] = LS(x, beta, S_ON, S_OFF, id_to_coords);
        sum_LS = sum_LS + LS_ON + LS_OFF;
    end
    
    [LS_ON, LS_OFF] = LS(x, alpha, S_ON, S_OFF, id_to_coords);
    dS_ON = LS_ON - (A(x - alpha) / (2 * sum_A)) * sum_LS;
    dS_OFF = LS_OFF - (A(x - alpha) / (2 * sum_A)) * sum_LS;
end
