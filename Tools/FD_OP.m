% ----------------------------------------------------------------------------------------------------------------------
% Constructs Finite Difference Operators u_x, u_xx in 1 and 2D
% ----------------------------------------------------------------------------------------------------------------------

classdef FD_OP
    
    methods(Static)
        
        % == START 1D OPERATORS ========================================================================================
        
        function LOP = UX_D(N) % U_x Homogenious Dirchlet
            
            function nnz = num_nonzero_elements(N)
                nnz = 2 * (N-2) + 2;  % number of nonzero entries: interior grid points (N-2 total points with  2 node stencils) + boundries (2 boundries of with 1 node stencils).
            end
            
            function [left, center, right] = weights(~, ~)
                center = 0;
                right  = 1;
                left   = -1;
            end
            
            LOP = FD_OP.construct_1d(@weights, @num_nonzero_elements, N);
            
        end
        
        function LOP = UX_N(N) % U_x Homogenious Neumann
            
            function nnz = num_nonzero_elements(N)
                nnz = 2 * (N-2);  % number of nonzero entries: interior grid points (N-2 total points with  2 node stencils) + boundries (2 boundries of with 1 node stencils).
            end
            
            function [left, center, right] = weights(i, N)
                % -- center --------------------------------------------------------------------------------------------
                center = 0;                
                % -- left ----------------------------------------------------------------------------------------------
                if(i == 1)
                    left = 0;
                else
                    left = -1;
                end
                % -- right ---------------------------------------------------------------------------------------------
                if(i == N)
                    right = 0;
                else
                    right = 1;
                end
            end
            
            LOP = FD_OP.construct_1d(@weights, @num_nonzero_elements, N);
            
        end
        
        function LOP = UXX_D(N) % U_xx Homogenious Dirchlet
            
            function nnz = num_nonzero_elements(N)
                nnz = 3 * (N-2) + 4;  % number of nonzero entries: interior grid points (N-2 total points with  3 node stencils) + boundries (2 boundries of with 2 node stencils).
            end
            
            function [left, center, right] = weights(~, ~)
                center = -2;
                right  = 1;
                left   = 1;
            end
            
            LOP = FD_OP.construct_1d(@weights, @num_nonzero_elements, N);
            
        end
        
        function LOP = UXX_N(N) % U_x Homogenious Neumann
            
            function nnz = num_nonzero_elements(N)
                nnz = 2 * (N-2) + 2;  % number of nonzero entries: interior grid points (N-2 total points with  2 node stencils) + boundries (2 boundries of with 1 node stencils).
            end
            
            function [left, center, right] = weights(i, N)
                % -- center --------------------------------------------------------------------------------------------
                center = -2;                
                % -- left ----------------------------------------------------------------------------------------------
                if(i == 1)
                    left = 2;
                else
                    left = 1;
                end
                % -- right ---------------------------------------------------------------------------------------------
                if(i == N)
                    right = 2;
                else
                    right = 1;
                end
            end
            
            LOP = FD_OP.construct_1d(@weights, @num_nonzero_elements, N);
            
        end
        
        % == START 2D OPERATORS ========================================================================================
        
        function LOP = UXX_UYY_D(N)
            % UXX_UYY_D: 2nd-order Finite Difference Operator for 2D Laplacian with all Dirchelt boundary conditions.
            % N represents number of interior grid points (not including boundary)
            
            % -- Set weights for 5-Point stencil -----------------------------------------------------------------------
            function [center, top, right, bottom, left] = weights(~, ~, ~)
                % -- center --------------------------------------------------------------------------------------------
                center = -4;
                top    = 1;
                right  = 1;
                bottom = 1;
                left   = 1;
            end
            
            function nnz = num_nonzero_elements(N)
                nnz = (N - 2) * (N - 2) * 5 + 4 * (N - 2) * 4 + 4 * 3; % number of nonzero entries: inner grid points ((N-2)(N-2) 5 node stencils) + grid boundries (4 sides of width N-2 with 4 node stencils) + corners (4 corners with 3 node stencils)
            end
            
            LOP = FD_OP.construct_2d(@weights, @num_nonzero_elements, N);
        end
        
        function LOP = UXX_UYY_N(N)
            % UXX_UYY_NMN: 2nd order FD operator for 2D Laplacian with all Neumann boundary conditions. N represents number of total grid
            % points (including boundary)
            %
            % Note: Implementation uses ghost nodes, and then reduce the new augemented system back into an N*N matrix.
            % See for example: R. J. LeVeque "Finite Difference Methods for Ordinary DIfferential Equations", Section 2.12 (p. 31).
            
            % -- Set weights for 5-Point stencil -----------------------------------------------------------------------
            function [center, top, right, bottom, left] = weights(i, j, N)
                % -- center --------------------------------------------------------------------------------------------
                center = -4;
                % -- top -----------------------------------------------------------------------------------------------
                if(i == N)
                    top = 2;
                else
                    top = 1;
                end
                % -- right ---------------------------------------------------------------------------------------------
                if(j == 1)
                    right = 2;
                else
                    right = 1;
                end
                % -- bottom --------------------------------------------------------------------------------------------
                if(i == 1)
                    bottom = 2;
                else
                    bottom = 1;
                end
                % -- left ----------------------------------------------------------------------------------------------
                if(j == N)
                    left = 2;
                else
                    left = 1;
                end
            end
            
            function nnz = num_nonzero_elements(N)
                nnz = (N - 2) * (N - 2) * 5 + 4 * (N - 2) * 4 + 4 * 3;  % number of nonzero entries: inner grid points ((N-2)(N-2) 5 node stencils) + grid boundries (4 sides of width N-2 with 4 node stencils) + corners (4 corners with 3 node stencils)
            end
            
            LOP = FD_OP.construct_2d(@weights, @num_nonzero_elements, N);
        end
        
        
        function LOP = UX_UY_N(N)
            % 2D FD Operator for U_x + U_y with all Neumann Boundary conditions. N represents number of total grid
            % points (including boundary)
            %
            % Note: Implementation uses ghost nodes, and then reduce the new augemented system back into an N*N matrix.
            % See for example: R. J. LeVeque "Finite Difference Methods for Ordinary DIfferential Equations", Section 2.12 (p. 31).
            
            % -- Set weights for 5-Point stencil -----------------------------------------------------------------------
            function [center, top, right, bottom, left] = weights(i, j, N)
                % -- center --------------------------------------------------------------------------------------------
                center = 0;
                % -- top -----------------------------------------------------------------------------------------------
                if(i == N)
                    top = 0;
                else
                    top = 1/2;
                end
                % -- right ---------------------------------------------------------------------------------------------
                if(j == 1)
                    right = 0;
                else
                    right = 1/2;
                end
                % -- bottom --------------------------------------------------------------------------------------------
                if(i == 1)
                    bottom = 0;
                else
                    bottom = -1/2;
                end
                % -- left ----------------------------------------------------------------------------------------------
                if(j == N)
                    left = 0;
                else
                    left = -1/2;
                end
            end
            
            function nnz = num_nonzero_elements(N)
                nnz = (N-2)^2 * 4 + 4 * (N-2) * 2;  % number of nonzero entries: interior grid points ((N-2)(N-2) total points with  4 node stencils) + inner boundries (4 boundries of N-2 points with 2 node stencils).
            end
            
            LOP = FD_OP.construct_2d(@weights, @num_nonzero_elements, N);
        end
        
        
        function LOP = UX_UY_D(N)
            % 2D FD Operator for U_x + U_y with all Dirchlet Boundary conditions. N represents number of interior grid
            % points (not including boundary)
            
            % -- Set weights for 5-Point stencil -----------------------------------------------------------------------
            function [center, top, right, bottom, left] = weights(~, ~, ~)
                % -- center --------------------------------------------------------------------------------------------
                center = 0;
                top = 1/2;
                right = 1/2;
                bottom = -1/2;
                left = -1/2;
            end
            
            function nnz = num_nonzero_elements(N)
                nnz = (N-2)^2 * 4 + 4 * (N-2) * 3 + 4 * 2;  % number of nonzero entries: interior grid points ((N-2)(N-2) total points with  4 node stencils) + inner boundries (4 boundries of N-2 points with 3 node stencils) + 4 corners with 2 node stencils.
            end
            
            LOP = FD_OP.construct_2d(@weights, @num_nonzero_elements, N);
        end
        
        
        % == START MATRIX CONSTRUCTION FUNCTIONS ===
        
        function LOP = construct_2d(weight_handle, nnz_handle, N)
            % CONSTRUCT_2D : Constructs linear operator for 5-point stencils in 2 dimensions.
            
            num_non_zero = nnz_handle(N);
            data = zeros(3, num_non_zero); % stores (i, j, v) pairs
            count = 1;
            
            for i = 1 : N
                for j = 1 : N
                    index = (i - 1) * N + j;
                    % --- CREATE STENCIL -------------------------------------------------------------------------------
                    
                    [w_center, w_top, w_right, w_bottom, w_left] = weight_handle(i, j, N);
                    
                    if(w_center ~= 0) % -- center connection -----------------------------------------------------------
                        data(:, count) = [index, index, w_center];
                        count = count + 1;
                    end
                    if (j < N && w_right ~= 0) % -- right connection ---------------------------------------------------
                        data(:, count) = [index, index + 1, w_right];
                        count = count + 1;
                    end
                    
                    if(j > 1 && w_left ~= 0) % -- left connection ------------------------------------------------------
                        data(:, count) = [index, index - 1, w_left];
                        count = count + 1;
                    end
                    
                    if(i > 1 && w_top ~= 0) % -- top connection --------------------------------------------------------
                        data(:, count) = [index, index - N, w_top];
                        count = count + 1;
                    end
                    
                    if(i < N && w_bottom ~= 0) % -- bottom connection --------------------------------------------------
                        data(:, count) = [index, index + N, w_bottom];
                        count = count + 1;
                    end
                    
                end
            end
            
            LOP = sparse(data(1,:), data(2, :), data(3, :), N * N, N * N);
            
        end
        
        function LOP = construct_1d(weight_handle, nnz_handle, N)
            % CONSTRUCT_1D : Constructs linear operator for 3-point stencils in 2 dimensions.
            
            num_non_zero = nnz_handle(N);
            data = zeros(3, num_non_zero); % stores (i, j, v) pairs
            count = 1;
            
            for i = 1 : N
                
                % --- CREATE STENCIL -----------------------------------------------------------------------------------
                [w_left, w_center, w_right] = weight_handle(i, N);
                
                if(w_center ~= 0) % -- center weight -------------------------------------------------------------------
                    data(:, count) = [i, i, w_center];
                    count = count + 1;
                end
                
                if(i > 1 && w_left ~= 0) % -- left weight ---------------------------------------------------------------
                    data(:, count) = [i, i - 1, w_left];
                    count = count + 1;
                end
                
                if(i < N && w_left ~= 0) % -- right weight -----------------------------------------------------------
                    data(:, count) = [i, i + 1, w_right];
                    count = count + 1;
                end
                
            end
            
            LOP = sparse(data(1,:), data(2, :), data(3, :), N, N);
            
        end
        
    end
    
end