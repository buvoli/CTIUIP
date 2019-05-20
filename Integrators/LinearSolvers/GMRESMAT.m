% ======================================================================================================================
% GMRESMAT : Linear solver that wraps Matlab's gmres
%
% -- Public Properties -------------------------------------------------------------------------------------------------
%
%  1. stats - generic stats object for linear solver 
%  2. tol - tolerance for GMRES solve
%
% -- Public Methods ----------------------------------------------------------------------------------------------------
%
%  1. solve(A, b, x0)    - Solve system of form A x = b using initial guess x0
%  2. solveBC(A,b,c, x0) - Solve system of form x = b + c * A * x using initial guess x0 
%
% ======================================================================================================================

classdef GMRESMAT < LinearSolver
    
    properties
        tol   = 1e-10;
        stats;
    end
    
    properties (Access = private)
        max_iterations;
    end
    
    methods
        
        function this = GMRESMAT(options)
            if(nargin == 0)
                options = struct();
            end
            options = this.DefaultOptions(options);
            this.max_iterations = options.max_iterations;
            this.tol = options.tolerance;
            this.stats = KrylovStats();
            this.stats.record = options.record_stats;
        end
        
        function [y, exit_flag, residual] = solve(this, A, b, x0)
            
            starting_time = tic;
        	if(isnumeric(A))
                real_valued = isreal(A) & isreal(b);
                if( ~ real_valued)
                    A_r = real(A);
                    A_i = imag(A);                    
                    A_aug  = [A_r -A_i; A_i A_r];
                    b_aug  = [real(b); imag(b)];
                    x0_aug = [real(x0); imag(x0)];
                    [y_aug, flag, residual, iter, resvec] = gmres(A_aug, b_aug, [], this.tol, min(length(x0_aug), this.max_iterations), [], [], x0_aug);
                    y = y_aug(1:end/2) + 1i*y_aug(end/2+1:end);
                else
                    [y, flag, residual, iter, resvec] = gmres(A, b, [], this.tol, min(length(x0), this.max_iterations), [], [], x0);
                end
            else
                [y, flag, residual, iter, resvec] = gmres(A, b, [], this.tol, min(length(x0), this.max_iterations), [], [], x0);
            end 
            if(flag == 0)
                exit_flag = true;
            else
                exit_flag = false;
            end
            this.stats.recordSolve(iter, resvec, toc(starting_time));
        end
        
        
        function [y, exit_flag, residual] = solveBC(this, A, b, c, x0)
            if(nargin == 4)
                x0 = zeros(size(b));
            end
            if(isnumeric(A))
                A_hat = speye(size(A)) - c * A;
                [y, exit_flag, residual] = solve(this, A_hat, b, x0);
            else
                A_hat = @(x) x - c * A(x);
                [y, exit_flag, residual] = solve(this, A_hat, b, x0);
            end
            
        end
        
    end
    
    methods(Access = private)
        function options = DefaultOptions(this, options)
            if(~isfield(options, 'tolerance'))
                options.tolerance = 1e-10;
            end
            if(~isfield(options, 'max_iterations'))
                options.max_iterations = 100;
            end
            if(~isfield(options, 'record_stats'))
                options.record_stats = false;
            end 
        end 
    end

end