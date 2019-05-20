% ======================================================================================================================
% ADR2D : Problem class for 2D ADR equation:
%
%       u_t = \epsilon * (u_{xx} + u_{yy}) + \alpha * (u_x + u_y) + \gamma * u .* (u - 1/2) .* (1 - u)
%
% NOTE: To learn about observables see: 
%   - https://www.mathworks.com/help/matlab/matlab_oop/listening-for-changes-to-property-values.html
%   - https://www.mathworks.com/help/matlab/matlab_oop/set-events-when-value-does-not-change.html
% ======================================================================================================================

classdef ADR2d < Problem

    properties(SetObservable)
        tspan     = [0, 0.05];
        params    = struct(     ...
            'N',        400,     ...
            'epsilon',  1/100,  ...
            'alpha',    -10,    ...
            'gamma',    100     ...      
        );        
    end
    
	properties(SetAccess = protected)
    	name = 'ADR2d';
        dimension;
        initial_condition;
        real_valued = true;
        description;
    end
    
    properties(Access = private)
    	linear_operator;
    end
        
    methods
        
        function desc = get.description(this)
            desc = sprintf('2D ADR, N = %d^2', this.params.N);
        end
        
        function up = RHS(this, u)
            N = this.params.gamma * u .* (u - 1/2) .* (1 - u);
            up = this.linear_operator * u + N;
        end

        function J = J(this, u)
            temp = this.params.gamma * (3*u - 1/2 - 3*u.^2);
            J = this.linear_operator + spdiags(temp, 0, this.dimension, this.dimension);
        end

        function N = N(this, u)
            N = this.params.gamma * u .* (u - 1/2) .* (1 - u);
        end
        
        function L = L(this)
            L = this.linear_operator;
        end       
        
    end
    
    methods(Access = protected)
        
        function reset(this, varargin)
            this.setInitialCondition();
            this.setLinearOperator();
        end
        
        function setDimension(this)
            this.dimension = this.params.N ^ 2;
        end
        
        function setInitialCondition(this)
            N = this.params.N;
            h = 1 / (N - 1);
            size = N * N;

            y0 = zeros(size,1);
            for j = 1:N
                y = (j-1) * h;
                for i = 1:N
                    x = (i-1) * h;
                    index = (j-1)*(N) + i;
                    y0(index) = 256 * ((1 - x) * x * (1 - y) * y)^2 + 0.3;
                end
            end
            
            this.initial_condition = y0;             
        end
        
        function setLinearOperator(this)
            epsilon = this.params.epsilon;
            alpha   = this.params.alpha;
            N       = this.params.N;
            h       = 1 / (N - 1);
            this.linear_operator = (epsilon / h^2) * FD_OP.UXX_UYY_N(N) + (alpha / h) * FD_OP.UX_UY_N(N);
        end
        
    end
    
end