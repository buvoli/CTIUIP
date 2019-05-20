% ======================================================================================================================
% BDFConst : BDF Methods of any order 
%
% -- Constructor Arguments --------------------------------------------------------------------------------------------- 
%
%   options (struct) - contains fields:
%       order -> positive integer pertaining to desired order of bdf method. Defaults to 2.
%       nonlinear_solver -> nonlinear solver to be used at each timestep. Defaults to Newton.
%
% ======================================================================================================================

classdef BDFConst < ImplicitIntegratorConst
    
    properties
        graph_line_style = {};
    end
    
    properties(SetAccess = protected)
        name  = '';
        description = '';
        order = [];
        starting_times   = []; 
        coefficients = [];
    end

    methods
        
        function this = BDFConst(options)
            if(nargin == 0)
            	options = struct();
            end
            default_options = {{'order', 2}};
            options = setDefaultOptions(options, default_options);
            this = this@ImplicitIntegratorConst(options);
            this.setOrder(options.order);            
        end
        
        function d = get.description(this)
            d = sprintf('Order %i BDF', this.order);
        end
        
        function setOrder(this, order)
            this.order            = order;
            this.name             = ['BDF', num2str(order)];
            this.coefficients     = reshape(double(polyCoefficients(0:order, 1:order, order+1, order)), [order + 1, 1]);
            this.starting_times   = 0:this.order-1;
        end     
        
    end
    
    methods (Access = protected)
        
        function setStepsize(this, problem)
            this.h = (problem.tspan(end) - problem.tspan(1))/(this.order - 1 + this.num_timesteps);
        end
        
        function [step_struct, y_in] = initStepStruct(~, ~, y_in, ~)
            step_struct = struct();
        end
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem, final_step)
            
            step_start_time = tic; % -- start step time clock ----------------------------------------------------------                 
            
            h = this.h; % stepsize
            q = this.order; % number of inputs
            % -- newton solve y = b + c * f(y) -------------------------------------------------------------------------
            c      = h * this.coefficients(q + 1); % nonlinear system leading coefficeint
            b      = y_in * this.coefficients(1 : q); % nonlinear rhs
            y_next = this.nonlinear_solver.solveBC(problem, b, c, y_in(:, end));
            % -- update state ------------------------------------------------------------------------------------------
            y_out = [y_in(:, 2 : end) y_next];
            t_out = t_in + h;
            % -- keep only final output at last step -------------------------------------------------------------------
            if(final_step) 
                t_out = t_out(q);
                y_out = y_out(:, q);
            end
            
            this.step_stats.recordStep(toc(step_start_time)); % -- stop step time clock --------------------------------
            
        end
        
    end

end