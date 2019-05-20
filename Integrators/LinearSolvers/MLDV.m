% ======================================================================================================================
% MLDV : Linear solver that wraps Matlab's backslash (mldv)
%
% -- Public Properties -------------------------------------------------------------------------------------------------
%
%  1. stats - generic stats object for linear solver 
%
% -- Public Methods ----------------------------------------------------------------------------------------------------
%
%  1. solve(A, b, ~)    - Solve system of form A x = b  
%  2. solveBC(A,b,c, ~) - Solve system of form x = b + c * A * x 
%
% ======================================================================================================================

classdef MLDV < LinearSolver
    
    properties
        stats;
    end
    
    methods
        
        function this = MLDV(options)
            if(nargin == 0)
                options = struct();
            end
            options = this.DefaultOptions(options);
            this.stats = MLDVStats();
            this.stats.record = options.record_stats;
        end
        
        function [x, exit_flag, residual] = solve(this, A, b, varargin)            
            starting_time = tic;
            x = A\b;
            exit_flag = true;
            if(nargout > 2)
                residual = norm(A*x - b, 2);
            end            
            this.stats.recordSolve(toc(starting_time));
        end
        
        function [x, exit_flag, residual] = solveBC(this, A, b, c, varargin)
            % solves (I - c * A) x = b
            starting_time = tic;
            n = size(A,1);
            x = (speye(n) - c * A)\b;  % (I - c * A) = b 
            exit_flag = true;
            if(nargout > 2)
                residual = norm(A*x - b, 2);
            end
            this.stats.recordSolve(toc(starting_time));
        end
       
    end
    
    methods(Access = private)
        function options = DefaultOptions(this, options)
            if(~isfield(options, 'record_stats'))
                options.record_stats = false;
            end            
        end
    end
end