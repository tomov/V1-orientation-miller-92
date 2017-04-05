classdef Simulator
    
    properties (Constant = true)
        maxX = 30;
        R = 5.5;    
        a = 1/2;
        r1 = 0.4;
        k = 0; % 'excit': 0, 'excit/inhibt': 1/9
        x1 = 2.5; % 'excit': 2.5, 'excit/inhibt': 7.5
        rc = 0.28; % 'excit': 0.28, 'excit/inhib': 0.2
        lambda_dt = 0.0012; % 'excit': 0.0012; 'excit/inhib': 0.0019
        
        norm = @(x) sqrt(sum(x.^2, 2)); % norm of each row of matrix
        kronecker = @(x) real(Simulator.norm(x) == 0); % delta function for >= 1D
    end
    
    methods (Static)    
        function res = G(x, r) % G from miller92
            R = Simulator.R;
            res = exp(- Simulator.norm(x) / (r * R)^2);
        end        
            
        function res = I(x) % I from miller92
            a = Simulator.a;
            r1 = Simulator.r1;
            k = Simulator.k;
            x1 = Simulator.x1;
            res = (a + (1 - a) * Simulator.kronecker(x)) .* (Simulator.G(x, r1) - k * Simulator.G(x, 3 * r1)) .* (Simulator.norm(x) <= x1);
        end
    end

    
    properties (Access = public)  
       X; % 0..30
       Y; % 0..30
       id_to_coords; % mapping of coordinates (offset by 1) to index
       coords_to_id; % mapping of index to coordinates
       coords_offset; % offset of coordinates when used as indices
       id_to_diff; % mapping of coordinate differences (offset by maxX + 1) to index
       diff_to_id; % mapping of index to coordinate differences
       diff_offset; % offset of differences when used as indices
       C_ON_ON;
       C_ON_OFF;
       C_OFF_ON;
       C_OFF_OFF;
       A;
       N; % number of neurons in each layer
    end
    
            
    methods
        function self = Simulator()
            maxX = Simulator.maxX;
            coords_offset = 1;
            diff_offset = maxX + 1;
            rc = Simulator.rc;
            
            % precalculate coordinates
            %
            [X,Y] = meshgrid(0:maxX); % MUST start form 0 -- b/c we're using it for differences too
            id_to_coords = [X(:) Y(:)];
            coords_to_id = zeros(size(X));
            for id = 1:size(id_to_coords, 1)
                coords_to_id(X(id) + coords_offset, Y(id) + coords_offset) = id;
                assert(id_to_coords(coords_to_id(X(id) + coords_offset, Y(id) + coords_offset), 1) == X(id));
                assert(id_to_coords(coords_to_id(X(id) + coords_offset, Y(id) + coords_offset), 2) == Y(id));
            end
            self.N = numel(X);
            self.X = X;
            self.Y = Y;
            self.id_to_coords = id_to_coords;
            self.coords_to_id = coords_to_id;
            self.coords_offset = coords_offset;
            
            % precalculate differences of coordinates
            %
            assert(id_to_coords(1, 1) == 0 && id_to_coords(1, 2) == 0);
            id_to_diff = [  id_to_coords(:,1),   id_to_coords(:,2);
                          - id_to_coords(:,1),   id_to_coords(:,2);
                            id_to_coords(:,1), - id_to_coords(:,2);
                          - id_to_coords(:,1), - id_to_coords(:,2)];
            diff_to_id = zeros(size(X) * 2); % offset by maxX + 1
            %{
            for id1 = 1:size(id_to_coords, 1)
                id1
                for id2 = 1:size(id_to_coords, 1)
                    a = id_to_coords(id1, :);
                    b = id_to_coords(id2, :);
                    diff = a - b;
                    id = find(ismember(id_to_diff, diff, 'rows'));
                    assert(length(id) == 1 || diff(1) == 0 || diff(2) == 0);
                    diff_to_id(diff(1) + diff_offset, diff(2) + diff_offset) = id(1);
                    assert(id_to_diff(diff_to_id(diff(1) + diff_offset, diff(2) + diff_offset), 1) == diff(1));
                    assert(id_to_diff(diff_to_id(diff(1) + diff_offset, diff(2) + diff_offset), 2) == diff(2));
                end
            end
            save(['Simulator_', num2str(maxX), '.mat']);
            %}
            % assumes it's been precalc'd -- takes a long time
            %
            load(['Simulator_', num2str(maxX), '.mat']);
            
            self.id_to_diff = id_to_diff;
            self.diff_to_id = diff_to_id;
            self.diff_offset = diff_offset;
            
            % precalculate C_ON_ON, etc
            %
            C_ON_ON = Simulator.G(id_to_diff, rc) - (1/9) * Simulator.G(id_to_diff, 3 * rc);
            C_ON_OFF = -0.5 * ones(size(C_ON_ON)); % ????????????
            self.C_ON_ON = C_ON_ON;
            self.C_OFF_OFF = C_ON_ON; % ??????????????
            self.C_ON_OFF = C_ON_OFF;
            self.C_OFF_ON = C_ON_OFF; % ??????????????
            
            % precalculate A
            %
            A = circle_intersect(Simulator.norm(id_to_diff), 5, 2.5) .* (Simulator.norm(id_to_diff) <= 5.5);
            self.A = A;
        end
        
        % A = C_ON_ON, A, or any other row vector with values corresponding
        %     to coordinate differences (see id_to_diff and diff_to_id)
        % diffs = differences (row vector) whose values we want to find in A
        %
        function res = get_diff(self, A, diffs)
            which = sub2ind(size(self.diff_to_id), diffs(:,1) + self.diff_offset, diffs(:,2) + self.diff_offset);
            assert(immse(self.id_to_diff(self.diff_to_id(which),:), diffs) == 0);
            res = A(self.diff_to_id(which));
        end
        
        % Initialize S_ON and S_OFF
        %
        function [S_ON, S_OFF] = initialize(self)
            N = self.N;
            rng('default');
            
            [x_ids, alpha_ids] = meshgrid(1:N);
            xs = self.id_to_coords(x_ids, :);
            alphas = self.id_to_coords(alpha_ids, :);
            
            S_ON = (rand(size(xs, 1), 1) * 0.4 + 0.8) .* self.get_diff(self.A, xs - alphas);
            S_OFF = (rand(size(xs, 1), 1) * 0.4 + 0.8) .* self.get_diff(self.A, xs - alphas);
            
            S_ON = reshape(S_ON, [N N]);
            S_OFF = reshape(S_OFF, [N N]);
        end
        
        % LS for given x, alpha
        % technically LS * dt
        %
        function [LS_ON, LS_OFF] = LS(self, x, alpha, S_ON, S_OFF)
            ys = self.id_to_coords;
            betas = self.id_to_coords;
            
            I_x_y = self.I(x - ys); % vector of I(x - y) for all y
            
            C_ON_ON_a_b = self.get_diff(self.C_ON_ON, alpha - betas); % vector of C(a - b) for all b
            C_ON_OFF_a_b = self.get_diff(self.C_ON_OFF, alpha - betas); % vector of C(a - b) for all b
            C_OFF_OFF_a_b = self.get_diff(self.C_OFF_OFF, alpha - betas); % vector of C(a - b) for all b
            C_OFF_ON_a_b = self.get_diff(self.C_OFF_ON, alpha - betas); % vector of C(a - b) for all b
            
            summands_ON = (I_x_y * C_ON_ON_a_b') .* S_ON + (I_x_y * C_ON_OFF_a_b') .* S_OFF;
            summands_OFF = (I_x_y * C_OFF_OFF_a_b') .* S_OFF + (I_x_y * C_OFF_ON_a_b') .* S_ON;
            
            sum_ON = sum(summands_ON(:));
            sum_OFF = sum(summands_OFF(:));
            
            LS_ON = self.lambda_dt * self.get_diff(self.A, x - alpha) * sum_ON;
            LS_OFF = self.lambda_dt * self.get_diff(self.A, x - alpha) * sum_OFF;
        end
        
        % LS for all x, alpha
        %
        function [LS_ON, LS_OFF] = LS_all(self, S_ON, S_OFF)
            N = self.N;
            
            ys = self.id_to_coords;
            betas = self.id_to_coords;
            
            for x_id = 1:N
                x = self.id_to_coords(x_id, :);
                I_x_y{x_id} = self.I(repmat(x, [N 1]) - ys); % vector of I(x - y) for all y
            end
            
            for alpha_id = 1:N
                alpha = self.id_to_coords(alpha_id, :);
                C_ON_ON_a_b{alpha_id} = self.get_diff(self.C_ON_ON, repmat(alpha, [N 1]) - betas); % vector of C(a - b) for all b
                C_ON_OFF_a_b{alpha_id} = self.get_diff(self.C_ON_OFF, repmat(alpha, [N 1]) - betas); % vector of C(a - b) for all b
                C_OFF_OFF_a_b{alpha_id} = self.get_diff(self.C_OFF_OFF, repmat(alpha, [N 1]) - betas); % vector of C(a - b) for all b
                C_OFF_ON_a_b{alpha_id} = self.get_diff(self.C_OFF_ON, repmat(alpha, [N 1]) - betas); % vector of C(a - b) for all b
            end
            
            LS_ON = zeros(size(S_ON));
            LS_OFF = zeros(size(S_OFF));
            for x_id = 1:N
                x = self.id_to_coords(x_id, :)
                for alpha_id = 1:N
                    alpha = self.id_to_coords(alpha_id, :);
                    
                    summands_ON = (I_x_y{x_id} * C_ON_ON_a_b{alpha_id}') .* S_ON + (I_x_y{x_id} * C_ON_OFF_a_b{alpha_id}') .* S_OFF;
                    summands_OFF = (I_x_y{x_id} * C_OFF_OFF_a_b{alpha_id}') .* S_OFF + (I_x_y{x_id} * C_OFF_ON_a_b{alpha_id}') .* S_ON;

                    sum_ON = sum(summands_ON(:));
                    sum_OFF = sum(summands_OFF(:));

                    A_x_a = self.get_diff(self.A, x - alpha);
                    LS_ON(x_id, alpha_id) = self.lambda_dt * A_x_a * sum_ON;
                    LS_OFF(x_id, alpha_id) = self.lambda_dt * A_x_a * sum_OFF;
                end
            end
        end
        
        % dS for given x, alpha
        %
        function [dS_ON, dS_OFF] = dS(self, x, alpha, S_ON, S_OFF)
            x_id = self.coords_to_id(x(1) + self.coords_offset, x(2) + self.coords_offset);
            alpha_id = self.coords_to_id(alpha(1) + self.coords_offset, alpha(2) + self.coords_offset);
                   
            [LS_ON, LS_OFF] = self.LS_all(S_ON, S_OFF);

            betas = self.id_to_coords;
            A_a_b = self.get_diff(self.A, x - betas); % vector of A(x - b) for all b
            sum_A = sum(A_a_b);
            
            sum_LS = sum(LS_ON(x_id, :)) + sum(LS_OFF(x_id, :));
            
            A_x_a = self.get_diff(self.A, x - alpha);
            dS_ON = LS_ON(x_id, alpha_id) - (A_x_a / (2 * sum_A)) * sum_LS;
            dS_OFF = LS_OFF(x_id, alpha_id) - (A_x_a / (2 * sum_A)) * sum_LS;
        end
        
        % dS for all x, alpha
        %
        function [dS_ON, dS_OFF] = dS_all(self, S_ON, S_OFF)
            [LS_ON, LS_OFF] = self.LS_all(S_ON, S_OFF);

            xs = self.id_to_coords;
            alphas = self.id_to_coords;
            betas = self.id_to_coords;
            
            % calculate sum over beta A(x - beta) for each x
            %
            sums_A = zeros(N, 1);
            for x_id = 1:N
                x = self.id_to_coords(x_id, :);
                A_a_b = self.get_diff(self.A, x - betas); % vector of A(x - b) for all b
                sum_A = sum(A_a_b);
                sums_A(x_id) = sum_A;
            end
            
            % calculate sum over beta LS(x, beta) for each x
            %
            sums_LS = zeros(N, 1);
            for x_id = 1:N
                sum_LS = sum(LS_ON(x_id, :)) + sum(LS_OFF(x_id, :));
                sums_LS(x_id) = sum_LS;
            end
            
            dS_ON = zeros(N, N);
            dS_OFF = zeros(N, N);
            for x_id = 1:N
                x = self.id_to_coords(x_id, :);
                for alpha_id = 1:N
                    alpha = self.id_to_coords(alpha_id, :);
                    A_x_a = self.get_diff(self.A, x - alpha);
                    dS_ON(x_id, alpha_id) = LS_ON(x_id, alpha_id) - (A_x_a / (2 * sums_A(x_id))) * sums_LS(x_id);
                    dS_OFF(x_id, alpha_id) = LS_OFF(x_id, alpha_id) - (A_x_a / (2 * sums_A(x_id))) * sums_LS(x_id);
                end
            end
        end
        
        
        
        %
        % --------------- EVERYTHING BELOW IS FOR SANITY CHECKS ------------
        %

        

        % compare slow & fast dS, just in case
        %
        function compare_dS(self)
            maxX = self.maxX;
            N = self.N;
            [S_ON, S_OFF] = self.initialize();
            
            % fast dS
            %
            tic
            [dS_ON, dS_OFF] = self.dS_all(S_ON, S_OFF);
            toc
            
            % SUPER slow dS
            %
            tic
            [dS_ON_sanity, dS_OFF_sanity] = self.dS_all_sanity(S_ON, S_OFF);
            toc
            
            % compare
            %
            x = [3 4];
            x_id = self.coords_to_id(x(1) + self.coords_offset, x(2) + self.coords_offset);
            Z = reshape(dS_ON(x_id, :), [maxX + 1 maxX + 1]);
            Z_sanity = reshape(dS_ON_sanity(x_id, :), [maxX + 1 maxX + 1]);

            figure;
            subplot(1, 2, 1);
            surf(Z);
            title('fast');
            subplot(1, 2, 2);
            surf(Z_sanity);
            title('slow');
        end
        
        % compare slow & fast LS, just in case
        %
        function compare_LS(self)
            maxX = self.maxX;
            N = self.N;
            [S_ON, S_OFF] = self.initialize();
            
            tic % fast
            [LS_ON, LS_OFF] = self.LS_all(S_ON, S_OFF);
            toc
            
            tic % slow
            [LS_ON_slow, LS_OFF_slow] = self.LS_all_slow(S_ON, S_OFF);
            toc
            
            tic % SUPER slow
            %[LS_ON_sanity, LS_OFF_sanity] = self.LS_all_sanity(S_ON, S_OFF);
            LS_ON_sanity = LS_ON_slow;
            LS_OFF_sanity = LS_OFF_slow;
            toc
            
            % compare
            %
            x = [3 4];
            x_id = self.coords_to_id(x(1) + self.coords_offset, x(2) + self.coords_offset);
            Z = reshape(LS_ON(x_id, :), [maxX + 1 maxX + 1]);
            Z_slow = reshape(LS_ON_slow(x_id, :), [maxX + 1 maxX + 1]);
            Z_sanity = reshape(LS_ON_sanity(x_id, :), [maxX + 1 maxX + 1]);

            figure;
            subplot(1, 3, 1);
            surf(Z);
            title('fast');
            subplot(1, 3, 2);
            surf(Z_slow);
            title('slow');
            subplot(1, 3, 3);
            surf(Z_sanity);
            title('sanity');
        end
        
        % very slow version of LS -- for sanity checking
        % technically LS * dt
        % also note the helper functions are listed there explicitly
        %
        function [LS_ON, LS_OFF] = LS_sanity(self, x, alpha, S_ON, S_OFF)
            N = self.N;

            R = self.R;
            G = @(x, r) exp(- norm(x) / (r * R)^2);
            kronecker = @(x) real(norm(x) == 0);

            a = self.a;
            r1 = self.r1;
            k = self.k;
            x1 = self.x1;
            I = @(x) (a + (1 - a) * kronecker(x)) .* (G(x, r1) - k * G(x, 3 * r1)) .* (norm(x) <= x1);

            rc = self.rc; % 'excit': 0.28, 'excit/inhib': 0.2
            C_ON_ON = @(x) G(x, rc) - (1/9) * G(x, 3*rc);
            C_ON_OFF = @(x) -0.5; % ???
            C_OFF_OFF = C_ON_ON; % ???????????
            C_OFF_ON = C_ON_OFF; % ???????????

            A = @(x) circle_intersect(norm(x), 5, 2.5) * (norm(x) <= 5.5);

            lambda_dt = self.lambda_dt; % 'excit': 0.0012; 'excit/inhib': 0.0019

            id_to_coords = self.id_to_coords;
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

        % very slow version of dS
        % for sanity checks
        %
        function [dS_ON, dS_OFF] = dS_sanity(self, x, alpha, S_ON, S_OFF)
            N = self.N;

            A = @(x) circle_intersect(norm(x), 5, 2.5) * (norm(x) <= 5.5);

            sum_A = 0;
            for beta_id = 1:N
                beta = self.id_to_coords(beta_id, :);
                sum_A = sum_A + A(x - beta);
            end

            sum_LS = 0;
            shit = [];
            for beta_id = 1:N
                beta = self.id_to_coords(beta_id, :);
                [LS_ON, LS_OFF] = self.LS_sanity(x, beta, S_ON, S_OFF);
                shit = [shit LS_ON];
                sum_LS = sum_LS + LS_ON + LS_OFF;
            end

            [LS_ON, LS_OFF] = self.LS_sanity(x, alpha, S_ON, S_OFF);
            dS_ON = LS_ON - (A(x - alpha) / (2 * sum_A)) * sum_LS;
            dS_OFF = LS_OFF - (A(x - alpha) / (2 * sum_A)) * sum_LS;
        end
        
        % VERY slow version of LS_all
        % for sanity checks
        %
        function [LS_ON, LS_OFF] = LS_all_sanity(self, S_ON, S_OFF)
            N = self.N;
            maxX = self.maxX;
            
            LS_ON = zeros(size(S_ON));
            LS_OFF = zeros(size(S_OFF));
            for x_id = 1:N
                x = self.id_to_coords(x_id, :);
                for alpha_id = 1:N
                    alpha = self.id_to_coords(alpha_id, :);
                    [LS_ON(x_id, alpha_id), LS_OFF(x_id, alpha_id)] = self.LS_sanity(x, alpha, S_ON, S_OFF);
                end
            end
            save(['LS_all_sanity_', num2str(maxX), '.mat']);
        end

        % slow version of LS_all
        %
        function [LS_ON, LS_OFF] = LS_all_slow(self, S_ON, S_OFF)
            N = self.N;
            
            LS_ON = zeros(size(S_ON));
            LS_OFF = zeros(size(S_OFF));
            for y_id = 1:N
                y = self.id_to_coords(y_id, :);
                for beta_id = 1:N
                    beta = self.id_to_coords(beta_id, :);
                    [LS_ON(y_id, beta_id), LS_OFF(y_id, beta_id)] = self.LS(y, beta, S_ON, S_OFF);
                end
            end
        end
    
        % VERY VERY slow version of dS_all
        % used for sanity checks
        %
        function [dS_ON, dS_OFF] = dS_all_sanity(self, S_ON, S_OFF)
            N = self.N;
            maxX = self.maxX;
            
            dS_ON_sanity = zeros(size(S_ON));
            dS_OFF_sanity = zeros(size(S_OFF));
            for x_id = 1:N
                x = self.id_to_coords(x_id, :)
                for alpha_id = 1:N
                    alpha = self.id_to_coords(alpha_id, :)
                    [dS_ON_sanity(x_id, alpha_id), dS_OFF_sanity(x_id, alpha_id)] = self.dS_sanity(x, alpha, S_ON, S_OFF);
                end
            end
            save(['dS_all_sanity_', num2str(maxX), '.mat']);
        end
        
        % Slow version of initialize
        %
        function [S_ON, S_OFF] = initialize_slow(self)
            N = self.N;
            rng('default');
            
            S_ON = zeros(N, N);
            S_OFF = zeros(N, N);
            for x_id = 1:N
                x = self.id_to_coords(x_id, :);
                for alpha_id = 1:N
                    alpha = self.id_to_coords(alpha_id, :);
                    S_ON(x_id, alpha_id) = (rand * 0.4 + 0.8) * self.get_diff(self.A, x - alpha);
                    S_OFF(x_id, alpha_id) = (rand * 0.4 + 0.8) * self.get_diff(self.A, x - alpha);
                end
            end            
        end

    end
end

