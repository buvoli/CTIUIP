% ======================================================================================================================
% AMConst : Adams-Moulton Methods of any order 
%
% -- Constructor Arguments --------------------------------------------------------------------------------------------- 
%
%   options (struct) - contains fields:
%       order -> positive integer pertaining to desired order of bdf method. Defaults to 2.
%       nonlinear_solver -> nonlinear solver to be used at each timestep. Defaults to Newton.
%
% ======================================================================================================================

classdef AMConst < ImplicitIntegratorConst
    
    properties
        graph_line_style = {};
        post_step_rhs_eval = false;
    end
    
    properties(SetAccess = protected)
        name  = '';
        description = '';
        order;
        starting_times;
        coefficients;
    end
    
    methods
        
        function this = AMConst(options)
            if(nargin == 0)
                options = struct();
            end
            default_options = {{'order', 2}};
            options = setDefaultOptions(options, default_options);
            this = this@ImplicitIntegratorConst(options);
            this.setOrder(options.order);          
        end
        
        function d = get.description(this)
            d = sprintf('Order %i Adams Moulton', this.order);
        end
        
        function setOrder(this, order)
            this.order = order;
            this.name  = ['AM', num2str(order)];
            this.starting_times = 0:max(0,this.order-2);
            this.coefficients = transpose(double(iPolyCoefficients(-(this.order - 2) : 1, [0 1])));
        end
            
    end
    
    methods (Access = protected)
        
        function setStepsize(this, problem)
            this.h = (problem.tspan(end) - problem.tspan(1))/(max(2,this.order) - 2 + this.num_timesteps);
        end
        
        function [step_struct, y_in] = initStepStruct(~, ~, y_in, problem)
        	% -- evaluate RHS ------------------------------------------------------------------------------------------
            F_in = zeros(size(y_in));
            for i = 1 : size(y_in, 2)
                F_in(:, i) = problem.RHS(y_in(:,i));
            end
            step_struct = struct('F_in', F_in);
        end
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem, final_step)
            
            step_start_time = tic; % -- start step time clock ----------------------------------------------------------                 
            
            h = this.h; % stepsize
            q = this.order; % number of inputs
            % -- newton solve y = b + c * f(y) -------------------------------------------------------------------------
           
            if(q == 1) % Backwards Euler
                b = y_in(:, end);
                c = h * this.coefficients;
            else % Adams Moulton
                b = y_in(:, end) + step_struct.F_in * (h * this.coefficients(1 : q-1));
                c = h * this.coefficients(q);
            end
            y_next = this.nonlinear_solver.solveBC(problem, b, c, y_in(:, end));
            step_struct.F_in(:, 1:end-1) = step_struct.F_in(:, 2:end);
            if(this.post_step_rhs_eval)
            	start_rhs_time = tic;
            	step_struct.F_in(:, end) = problem.RHS(y_next);
              	this.rhs_stats.recordRHSEval(toc(start_rhs_time));
            else
            	step_struct.F_in(:, end) = (y_next(:, end) - b) / c;
            end
            % -- update state ------------------------------------------------------------------------------------------
            y_out = [y_in(:, 2 : end) y_next];
            t_out = t_in + h;
            % -- keep only final output at last step -------------------------------------------------------------------
            if(final_step) 
                t_out = t_out(end);
                y_out = y_out(:, end);
            end
            
            this.step_stats.recordStep(toc(step_start_time)); % -- stop step time clock --------------------------------
            
        end
        
    end
    
end