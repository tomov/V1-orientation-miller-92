classdef Simulator
    
    properties (Constant = true)
        maxX = 5;
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
            %{
            assert(id_to_coords(1, 1) == 0 && id_to_coords(1, 2) == 0);
            id_to_diff = [  id_to_coords(:,1),   id_to_coords(:,2);
                          - id_to_coords(:,1),   id_to_coords(:,2);
                            id_to_coords(:,1), - id_to_coords(:,2);
                          - id_to_coords(:,1), - id_to_coords(:,2)];
            diff_to_id = zeros(size(X) * 2); % offset by maxX + 1
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
            res = A(self.diff_to_id(diffs(:,1) + self.diff_offset, diffs(:,2) + self.diff_offset));
        end

        function [LS_ON, LS_OFF] = LS(self, x, alpha, S_ON, S_OFF)
            rc = Simulator.rc;
            maxX = Simulator.maxX;
            
            y = self.id_to_coords;
            I_x_y = I(x - y);
            
            beta = self.id_to_coords;
            C_a_b = self.get_diff(self.C_ON_ON, alpha - beta);
            
            
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

    end
    
end

