% ======================================================================================================================
% LINEARSOLVER : A Basic implementation of Newtons Method
%
% -- Public Properties -------------------------------------------------------------------------------------------------
%
%  1. max_iterations - newton iteration will terminate with exit_flag = false if iterations exceed this quantity
%  2. max_residual_tolerance - newton iteration will terminate with exit_flag = false if residual excedes this quantity
%  3. residual_tolerance - newton interation will terminate with exit_flag = true if residual is below this quantity.
%  4. delta_tolerance - newton interation will terminate with exit_flag = true if change in residual is below this quantity.
%  5. tol_norm - norm for determining tolerances 
%  6. linear_solver - linear solver. Must inherit from LinearSolver
%
% -- Public Methods ----------------------------------------------------------------------------------------------------
%
%  1. solve   - solve system of form F(x) = 0  
%  2. solveBC - solve system of form Solve system of form x = b + c * F(x) 
%
% ======================================================================================================================


classdef Newton < NonlinearSolver
    
    properties
        max_iterations
        max_residual_tolerance % iteration exits if residual exceded quantity
        residual_tolerance
        delta_tolerance
        tol_norm
        linear_solver
        stats
    end
    
    methods
    
        function this = Newton(options)
            if(nargin == 0)
                options = struct();
            end
            this.stats = NewtonStats();
            options = this.DefaultOptions(options);
            this.max_iterations = options.max_iterations;
            this.residual_tolerance = options.residual_tolerance;
            this.max_residual_tolerance = options.max_residual_tolerance;
            this.delta_tolerance = options.delta_tolerance;
            this.tol_norm = options.tol_norm;
            this.linear_solver = options.linear_solver;
            this.stats.record = options.record_stats;
        end
        
        function storeStats(this, flag)
            this.stats.store = flag;
        end
    
        function [x_k, clean_exit_flag, final_residual] = solve(this, problem, x0)
        % SOLVEBC solves system x = b + c * F(x)
        % PARAMETERS
        %   problem (class or struct)   : must have the parameters:
        %                                   varArgs (Cell)   - additional arguments for F and J
        %                                   F       (Handle) - function for F(x)
        %                               and one of the following
        %                                   Jxy (Handle) - function @(x,y) returning product of jacobian at x with vector y        
        %                                   J (Handle) - function @(x) returning matrix representing the jacobian of F at x
        %   b       (vector)            : constant b in nonlinear system
        %   c       (double)            : constant c in nonlinear system
        %   x0      (vector)            : initial guess
            
        % -- read problem properties --------
        dim = length(x);
        
        % -- determine if jacobian is provided explicitly or if multiplication function exists
        if(isfield(problem, 'J') || isprop(problem, 'J'))
            explicit_jacobian = true;
            e = ones(dim, 1);
        elseif(isfield(problem, 'Jxy') || isprop(problem, 'Jxy'))
            explicit_jacobian = false;
        else
            error('Jacobian or Jacobian mult function not provided');
        end
        
        % -- Newton Iteration ---------
        iterations = 0; 
        x_k        = x0;
        residuals  = zeros(this.max_iterations + 1, 1);
        deltas     = zeros(this.max_iterations, 1);
        seconds    = zeros(this.max_iteration, 1);
        clean_exit_flag = false;
        
        while(true) % DOWHILE loop
                iteration_start_time = tic;
                % -- Evaluate NL Function G(y) = b + c * f(y) - y ------------------------------------------------------
                G  = problem.RHS(x_k);
                % -- Calculate residual and test exit conditions -------------------------------------------------------
                residual_index = iterations + 1;
                current_residual = this.tol_norm(G);
                residuals(residual_index) = current_residual;               
                if (current_residual < this.residual_tolerance) % clean residual exit condition
                    clean_exit_flag = true;
                    break;
                end
                if(iterations >= this.max_iterations || current_residual >= this.max_residual_tolerance || isnan(current_residual)) % bad exit conditions
                   break
                end
                % -- Construct Jacobian and Solve ----------------------------------------------------------------------
                if(explicit_jacobian)
                    GP  = problem.J(x_k);
                    rhs = GP * x_k - G;
                else
                    GP = @(y) problem.Jxy(x_k, y);
                    rhs = GP(x_k) - G;
                end                
                x_km1 = x_k;
                x_k = this.linear_solver.solve(GP, rhs, x_k);
                % -- increment iteration count and test exit conditions ------------------------------------------------
                iterations = iterations + 1;
                deltas(iterations) = this.tol_norm(x_k - x_km1);
                seconds(iterations)  = toc(iteration_start_time);
                if(deltas(iterations) < this.delta_tolerance) % clean delta exit condition
                    clean_exit_flag = true;
                    break;
                end
                
        end
        final_residual = residuals(residual_index);
        this.stats.addSolve(iterations, residuals(1:iterations+1), deltas(1:iterations), seconds(1:iterations))
        end
        
        function [x_k, clean_exit_flag, final_residual] = solveBC(this, problem, b, c, x0)
        % SOLVEBC solves system x = b + c * F(x)
        % PARAMETERS
        %   problem (class or struct)   : must have the parameters:
        %                                   varArgs (Cell)   - additional arguments for F and J
        %                                   F       (Handle) - function for F(x)
        %                               and one of the following
        %                                   Jxy (Handle) - function @(x,y) returning product of jacobian at x with vector y        
        %                                   J (Handle) - function @(x) returning matrix representing the jacobian of F at x
        %   b       (vector)            : constant b in nonlinear system
        %   c       (double)            : constant c in nonlinear system
        %   x0      (vector)            : initial guess
            
        % -- read problem properties --------
        dim = length(x0);
        
        % -- determine if jacobian is provided explicitly or if multiplication function exists
        if(isfield(problem, 'J') || ismethod(problem, 'J') || isprop(problem, 'J'))
            explicit_jacobian = true;
            e = ones(dim, 1);
        elseif(isfield(problem, 'Jxy') || ismethod(problem, 'Jxy') || isprop(problem, 'Jxy'))
            explicit_jacobian = false;
        else
            error('Jacobian or Jacobian mult function not provided');
        end
        
        % -- Newton Iteration ---------
        iterations = 0; 
        x_k        = x0;
        residuals  = zeros(this.max_iterations + 1, 1);
        deltas     = zeros(this.max_iterations, 1);
        seconds    = zeros(this.max_iterations, 1);
        clean_exit_flag = false;
        
        while(true) % DOWHILE loop
                iteration_start_time = tic;
                % -- Evaluate NL Function G(y) = b + c * f(y) - y ------------------------------------------------------
                G  = b + c * problem.RHS(x_k) - x_k;
                % -- Calculate residual and test exit conditions -------------------------------------------------------
                residual_index = iterations + 1;
                current_residual = this.tol_norm(G);
                residuals(residual_index) = current_residual;               
                if (current_residual < this.residual_tolerance) % clean residual exit condition
                    clean_exit_flag = true;
                    break;
                end
                if(iterations >= this.max_iterations || current_residual >= this.max_residual_tolerance || isnan(current_residual)) % bad exit conditions
                   break
                end                
                % -- Construct Jacobian and Solve ----------------------------------------------------------------------
                if(explicit_jacobian)
                    GP  = c * problem.J(x_k) - spdiags(e, 0, dim, dim);
                    rhs = GP * x_k - G;
                else
                    GP = @(y) c * problem.Jxy(x_k, y) - y;
                    rhs = GP(x_k) - G;
                end
                x_km1 = x_k;
                x_k   = this.linear_solver.solve(GP, rhs, x_k);
                % -- increment iteration count and test exit conditions ------------------------------------------------
                iterations = iterations + 1;
                deltas(iterations) = this.tol_norm(x_k - x_km1);
                seconds(iterations)  = toc(iteration_start_time);
                if(deltas(iterations) < this.delta_tolerance) % clean delta exit condition
                    clean_exit_flag = true;
                    break;
                end
                
        end
        final_residual = residuals(residual_index);
        this.stats.addSolve(iterations, residuals(1:iterations+1), deltas(1:iterations), seconds(1:iterations))
        end        
     
    end
    
    methods(Access = private)
        
        function options = DefaultOptions(this, options)
            if(~isfield(options, 'max_iterations'))
                options.max_iterations = 100;
            end
            if(~isfield(options, 'residual_tolerance'))
                options.residual_tolerance = 1e-12;
            end
            if(~isfield(options, 'max_residual_tolerance'))
                options.max_residual_tolerance = Inf;
            end
            if(~isfield(options, 'delta_tolerance'))
                options.delta_tolerance = 1e-12;
            end
            if(~isfield(options, 'tol_norm'))
                options.tol_norm = @(x) norm(x, Inf);
            end
            if(~isfield(options, 'linear_solver'))
                options.linear_solver = MLDV();
            end
            if(~isfield(options, 'record_stats'))
                options.record_stats = false;
            end  
        end    
        
    end
    
end

