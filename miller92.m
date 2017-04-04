maxX = 31; % dimension of grid

[X,Y] = meshgrid(1:1:31);
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

for x_id = 1:N
    x = id_to_coords(x_id, :);
    for alpha_id = 1:N
        alpha = id_to_coords(alpha_id, :);
        S_ON(x_id, alpha_id) = (rand * 0.4 + 0.8) * A([y1 y2]);
        S_OFF(y1, y2) = (rand * 0.4 + 0.8) * A([y1 y2]);
    end
end



x = LS([3 4], [1 2], S_ON, S_OFF);

function [LS_ON, LS_OFF] = LS(x, alpha, S_ON, S_OFF)
    maxX = 31;
    
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
    for y1 = 1:maxX
        for y2 = 1:maxX
            y = [y1 y2];
            for beta1 = 1:maxX
                for beta2 = 1:maxX
                    beta = [beta1 beta2];
                    sum_ON = sum_ON + C_ON_ON(alpha - beta) * S_ON(y, beta) + ...
                                      C_ON_OFF(alpha - beta) * S_OFF(y, beta);
                    sum_OFF = sum_OFF + C_OFF_OFF(alpha - beta) * S_OFF(y, beta) + ...
                                        C_OFF_ON(alpha - beta) * S_ON(y, beta);
                end
            end
        end
    end
    LS_ON = lambda_dt * A(x - alpha) * sum_ON;
    LS_OFF = lambda_dt * A(x - alpha) * sum_OFF;
end

function [dS_ON, dS_OFF] = dS(x, alpha, LS_ON, LS_OFF)
    maxX = 31;
    
    A = @(x) circle_intersect(norm(x), 5, 2.5) * (norm(x) <= 5.5);
    
    sum_norm_A = 0;
    for beta1 = 1:maxX
        for beta2 = 1:maxX
            beta = [beta1 beta2];
            sum_norm_A = sum_norm_A + A(x - beta);
        end
    end
    
    sum_norm_LS = 0;
    for beta1 = 1:maxX
        for beta2 = 1:maxX
            beta = [beta1 beta2];
            [LS_ON, LS_OFF] = LS(x, beta);
            sum_norm_LS = sum_norm_LS + LS_ON + LS_OFF;
        end
    end
    
    dS_ON = LS_ON(x, alpha) - (A(x - alpha) / 2 * sum_norm_A) * sum_norm_LS;
    dS_OFF = LS_OFF(x, alpha) - (A(x - alpha) / 2 * sum_norm_A) * sum_norm_LS;   
end
