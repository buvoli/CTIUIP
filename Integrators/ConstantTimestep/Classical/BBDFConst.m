% ======================================================================================================================
% BBDFConst : Serial implementation of block BDF methods (BBDF) of any order.
%
% -- Constructor Arguments --------------------------------------------------------------------------------------------- 
%
%   options (struct) - contains fields:
%
%       z      -> (vector) vector of nodes z_j.
%       alpha  -> (double) nonlinear solver to be used at each timestep. Defaults to Newton.
%       z_name -> (string) optional string describing nodes
%       nonlinear_solver -> nonlinear solver to be used at each timestep. Defaults to Newton.
%
% ======================================================================================================================

classdef BBDFConst < ImplicitIntegratorConst
    
    properties
        graph_line_style = {};
    end
    
    properties(SetAccess = protected)
        name  = '';
        description = '';
        order = [];
        starting_times      = [];
        pbdf_coefficients   = [];
        output_coefficients = [];
        output_index        = [];
        z        = [];
        alpha    = [];
        r        = [];
        z_output = 0;
    end
    
    properties(Access = private)
        z_name = '';
    end
    
    methods
        
        function this = BBDFConst(options)
            if(nargin == 0)
                options = struct();
            end
            default_options = {{'z', 0},{'alpha', 1},{'z_name', ''}};
            options = setDefaultOptions(options, default_options);
            
            this = this@ImplicitIntegratorConst(options);
            this.setParameters(options.z, options.alpha, options.z_name);
        end
        
        function d = get.description(this)
            [an, ad] = rat(this.alpha);
            if(ad == 1)
                alpha_str = sprintf('alpha = %i', an);
            else
                alpha_str = sprintf('alpha = %i/%i', an, ad);
            end           
            if(isempty(this.z_name))
                d = sprintf('Order %i Block BDF (%s)', this.order, alpha_str);
            else
                d = sprintf('Order %i Block BDF (%s, z = %s)', this.order, alpha_str, this.z_name);
            end
        end
        
        function setParameters(this, z, alpha, z_name)
            if(nargin < 4)
                z_name = '';
            end
            this.z = z(:);
            this.z_name = z_name;
            this.alpha = alpha;
            this.starting_times = double(z(:) / alpha);
            this.order = length(z);
            calculateCoefficients(this);            
            % -- adjust stat objects -----------------------------------------------------------------------------------
            this.nonlinear_solver.stats.reset(length(z));
            this.nonlinear_solver.linear_solver.stats.reset(length(z)); 
            this.step_stats.reset(length(z)); 
        end
    
    end
    
    methods (Access = private)
        
        function setName(this)
            if(~isempty(this.z_name))
                this.name = ['BBDF', num2str(length(this.z)), '_{', this.z_name, '}'];
            else
                this.name = ['BBDF', num2str(length(this.z))];
            end
        end
        
        function [solve_flags, conj_get, conj_set] = solveIndices(this, problem)
            % SOLVEINDICES: computes indices of z that require a nonlinear solve and specifies which indices of z can
            % be obtained via swartz reflection. y^{[n+1]}(conj_set) = nl_solve(conj_set)
            % Returns:
            % solve_flags   (array) : 	array of size q specifying nonlinear solve for jth output:
            %                               value of 2 -> complex solve, 
            %                               value of 1 -> real solve
            %                               value of 0 -> no solve (reflected value)
            % conjugate_get (array) :   indices of solve_index that can be computed via swartz reflection
            % conjugate_set (array) :   indices of z that can be computed via via swartz reflection
            
            zd = double(this.z);
            q = length(zd);
            conj_get = [];
            conj_set = {};
            abs_tol  = 1e-14; % neceessary or solve indices will miss classify solution points due to rouding error
            
            if(problem.real_valued)                
                real_index   = find(abs(imag(zd)) < abs_tol);
                solve_flags  = zeros(1, q);
                    solve_flags(real_index) = 1;                
                search_queue = 1:q;                
                    search_queue(real_index) = [];
                while(~isempty(search_queue))
                    % pop first value
                    z_index = search_queue(1);
                    search_queue(1) = [];
                    % add to complex solve
                    solve_flags(z_index) = 2;
                    % find all reflected nodes
                    z_new   = zd(z_index);
                    reflected_z_index = find(abs(zd(search_queue) - conj(z_new)) <= abs_tol);% find(z(search_queue) == conj(z_new)); REPLACED due to rounding errors
                    if(~isempty(reflected_z_index))
                        conj_get(end+1) = z_index;
                        conj_set{end+1} = search_queue(reflected_z_index);
                        search_queue(reflected_z_index) = [];
                    end
                end
            else 
                solve_flags = 2 * ones(1, q);
            end
            
        end
        
        function calculateCoefficients(this)
            if(~isempty(this.alpha) && ~isempty(this.z))
                q   = length(this.z);
                tau_out = this.z + this.alpha;
                C   = zeros(q, q + 1);
                for j=1:q
                    zeta_j = [this.z; tau_out(j)];
                    A_j    = 1:q;
                    B_j    = q+1;
                    C(j,:) = double(polyCoefficients(zeta_j, A_j, B_j, tau_out(j)));
                end
                this.pbdf_coefficients = C;
                % -- set output coefficients --
                [~, index] = ismember(this.z_output, this.z);
                if(index)
                    this.output_index = index;
                    this.output_coefficients = [];  % interpolation formula
                else
                    this.output_index = [];
                    zeta_out = [this.z; this.z_output + this.alpha];
                    tau_out  = this.z_output + this.alpha;
                    this.output_coefficients = double(polyCoefficients(zeta_out, 1:q, q + 1, tau_out));  % nonlinear formula
                end
                this.setName();
            end
        end
        
 
        
    end
    
    methods (Access = protected)
        
        
        function [step_struct, y_in] = initStepStruct(this, t_in, y_in, problem)
            
            [solve_flags, conj_get, conj_set] = solveIndices(this, problem);
            
            step_struct = struct(                                    ...
                'coeff',        transpose(this.pbdf_coefficients),   ...
                'output_coeff', transpose(this.output_coefficients), ...
                'solve_flags',  solve_flags,                         ...
                'conj_get',     conj_get,                            ...
                'conj_set',     {conj_set}                           ...
                );
            
        end
        
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem, final_step)
            
            % -- store class variables ---------------------------------------------------------------------------------
            h         = this.h;
            alpha     = double(this.alpha);
            r         = h / alpha;
            z         = double(this.z);
            q         = length(this.z);
            
            % -- store step_struct variables ---------------------------------------------------------------------------
            coeff = step_struct.coeff;
            solve_flags  = step_struct.solve_flags;
            
            % -- create empty time record for current step -------------------------------------------------------------
            for j = 1 : q
                this.step_stats.setIndex(j);
                this.step_stats.recordStep(0);
            end
            
            if(final_step)
                
                this.step_stats.setIndex(1);
                output_coeff = step_struct.output_coeff;
                output_index = this.output_index;
                [~, output_y0_index] = min(abs(this.z - this.z_output)); % closest initial condition
                
                if(isempty(output_index)) % output node is not a member of node set
                    y0 = y_in(:, output_y0_index);           % initial gues
                    b  = y_in * output_coeff(1 : q);         % vector b
                    c  = r * output_coeff(q + 1);            % constant
                else % output node is a member of node set
                    y0 = y_in(:, output_index);              % initial guess
                    b  = y_in * coeff(1 : q, output_index);  % vector b
                    c  = r * coeff(q + 1, output_index);     % constants
                end
                
                t_out = t_in + h;
                y_out = this.nonlinear_solver.solveBC(problem, b, c, y0);
                if(problem.real_valued)
                    y_out = real(y_out);
                end
                
            else
                
                y_out = zeros(size(y_in));
                % -- solve nonlinear systems ---------------------------------------------------------------------------
                for j = 1 : q
                    if(solve_flags(j) > 0)
                        output_start_time = tic;
                        this.step_stats.setIndex(j);
                        this.nonlinear_solver.stats.setIndex(j);
                        this.nonlinear_solver.linear_solver.stats.setIndex(j);
                        b = y_in * coeff(1 : q, j);
                        c = r * coeff(q + 1, j);
                        if(solve_flags(j) == 1) % real valued ----------------------------------------------------------
                            [y_out(:,j), clean_exit] = this.nonlinear_solver.solveBC(problem, real(b), real(c), real(y_in(:,j)));
                        else % complex valued --------------------------------------------------------------------------
                            [y_out(:,j), clean_exit] = this.nonlinear_solver.solveBC(problem, b, c, y_in(:,j));
                        end
                        this.step_stats.recordSubstep(toc(output_start_time));
                        % -- test for emergency exit condition ---------------------------------------------------------
                        emergency_exit = ((~clean_exit) || norm(y_out(:,j), Inf) >= Inf || any(isnan(y_out(:,j))));
                        if(emergency_exit)
                            y_out = NaN; t_out = t_in;
                            return;
                        end
                    end
                end
                
                % -- count remaining step time on first output coutner -------------------------------------------------
                this.step_stats.setIndex(1);
                step_start_time = tic; % remaining step time
                
                % -- reflected values ----------------------------------------------------------------------------------
                conj_get = step_struct.conj_get;
                conj_set = step_struct.conj_set;
                for j = 1:length(conj_get)
                    y_out(:,conj_set{j}) = repmat( conj(y_out(:,conj_get(j))), 1, length(conj_set{j}));
                end
                
                % -- update state --------------------------------------------------------------------------------------
                t_out = t_in + h;
                this.step_stats.recordSubstep(toc(step_start_time));
                
            end
            
            
        end
          
    end
    
end