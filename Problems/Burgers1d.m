% ======================================================================================================================
% BURGERS1D : Problem class for 1D Burgers' equation:
%
%       u_t = \nu * u_{xx} + (u^2)_x
%
% NOTE: To learn about observables see: 
%   - https://www.mathworks.com/help/matlab/matlab_oop/listening-for-changes-to-property-values.html
%   - https://www.mathworks.com/help/matlab/matlab_oop/set-events-when-value-does-not-change.html
% ======================================================================================================================

classdef Burgers1d < Problem
    
    properties(SetObservable)
        tspan     = [0, 1];
        params    = struct(     ...
            'N',        2000,   ...
            'Lx',       1,      ...
            'nu',       3e-4    ...      
        );        
    end

	properties(SetAccess = protected)
    	name = 'Burgers1d';
        dimension;
        initial_condition;
        real_valued = true;
        description;
    end
    
    properties(Access = private)
    	linear_operator;
        ux_operator; % stores linear operator for (-1/2) d/dx 
    end

    methods
        
        function desc = get.description(this)
            desc = sprintf('Burgers 1d, \nu = %2.2g, N = %d', this.params.nu, this.params.N);
        end

        function yp = RHS(this, u)
            yp = this.linear_operator * u + this.N(u);
        end

        function J = J(this, u)
            N = this.params.N;
            h = this.params.Lx / (N + 1);
            c = 1 / (2*h);

            J = c * (spdiags(u, -1, N, N) - spdiags(u, 1, N, N)) + this.linear_operator;
        end

        function N = N(this, u)
            N = this.ux_operator * (u.^2);
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
            this.dimension = this.params.N;
        end
        
        function setInitialCondition(this)
            N = this.params.N;
            x = linspace(0, this.params.Lx, this.params.N + 2);
            y0 = zeros(N, 1);
            for i = 2 : N + 1
                y0(i-1,1) = (sin(3 * pi * x(i))) .^ 2 .* (1 - x(i)) .^ (3/2);
            end
            this.initial_condition = y0;            
        end
        
        function setLinearOperator(this)
            N = this.params.N;
            h = this.params.Lx / (N + 1);
            this.linear_operator = (this.params.nu / h^2) * FD_OP.UXX_D(N);
            this.ux_operator = -1 / (4 * h) * FD_OP.UX_D(N); % store operator for computing (-1/2) * (d/dx) - this is used by nonlinear evaluation
        end
        
    end
    
end