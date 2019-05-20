% ======================================================================================================================
% BBDFConst_spmd : A parallel implementation of block BDF Methods (BBDF) of any order, using Matlab's parallel toolkit.
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


classdef BBDFConst_spmd < ImplicitIntegratorConst
    
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
        
        function this = BBDFConst_spmd(options)
            if(nargin == 0)
                options = struct();
            end
            default_options = {{'z', 0},{'alpha', 1},{'z_name', ''}};
            options = setDefaultOptions(options, default_options);
            
            this = this@ImplicitIntegratorConst(options);
            this.setParameters(options.z, options.alpha, options.z_name);
        end
        
        % redefine solve to allow for spmd
        function [t_out, y_out] = solve(this, problem)
            if(isempty(this.num_timesteps))
                error('Integrator: number of timesteps has not been set.');
            end
            % -- set stepsize ------------------------------------------------------------------------------------------
            this.setStepsize(problem);
            % -- obtain initial conditions -----------------------------------------------------------------------------
            [t_in, y_in, clean_exit] = this.initialConditions(problem);
            if(~clean_exit)
                t_out = NaN;
                y_out = NaN;
                warning('Integrator: Emergency exit during initial condition computation.');
                return;
            end
            % -- solve problem -----------------------------------------------------------------------------------------
            start_step_time = tic;
            [struct_in, y_in] = this.initStepStruct(t_in, y_in, problem);
            num_steps = this.num_timesteps;
            num_proc  = struct_in.num_proc;
            spmd(num_proc)
                if(labindex < num_proc)
                    y_in = []; t_in = []; % helper processors do not propagate solution.
                end               
                for i = 1 : num_steps
                    final_step_flag = (i == num_steps);              
                    [t_in, y_in, struct_in] = this.step(t_in, y_in, struct_in, problem, final_step_flag);
                    % -- emergency exit conditions ---------------------------------------------------------------------
                    if(any(isinf(y_in(:))) || any(isnan(y_in(:))))
                        t_in = NaN;
                        y_in = NaN;
                        warning('Integrator: Emergency exit at step %i', i);
                        break;
                    end
                end
            end
            t_out = t_in; y_out = y_in{num_proc};
            this.step_stats.recordStep(toc(start_step_time));
            % -- output solution if only one output argument -----------------------------------------------------------
            if(nargout <= 1)
                t_out = y_out;
            end
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
            this.nonlinear_solver.stats.reset();
            this.nonlinear_solver.linear_solver.stats.reset(); 
            this.step_stats.reset(); 
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
            solve_inds = find(solve_flags > 0);
            num_solves = sum(solve_flags ~= 0);
            
            step_struct = struct(                                    ...
                'coeff',        transpose(this.pbdf_coefficients),   ...
                'output_coeff', transpose(this.output_coefficients), ...
                'solve_flags',  solve_flags,                         ...
                'solve_inds',   solve_inds,                          ...
                'conj_get',     conj_get,                            ...
                'conj_set',     {conj_set},                          ...
                'num_proc',     1 + num_solves                       ...
            );
            
        end
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem, final_step)
            
            % -- store class variables ---------------------------------------------------------------------------------
            h         = this.h;
            alpha     = double(this.alpha);
            r         = h / alpha;
            z         = double(this.z);
            q         = length(this.z);
            output_index = this.output_index;
            
            % -- store step_struct variables ---------------------------------------------------------------------------
            coeff        = step_struct.coeff;
            output_coeff = step_struct.output_coeff;
            solve_flags  = step_struct.solve_flags;
            solve_inds   = step_struct.solve_inds;
            conj_get     = step_struct.conj_get;
            conj_set     = step_struct.conj_set;
            num_proc     = step_struct.num_proc;
            num_solves   = num_proc - 1;
            
            if(final_step)
                if(labindex == num_proc)  % -- compute output on leader ------------------------------------------------              
                    if(isempty(output_index)) % output node is not a member of node set
                        [~, output_y0_index] = min(abs(this.z - this.z_output)); % index of closest initial condition
                        y0 = y_in(:, output_y0_index);            % initial gues
                        b  = y_in * output_coeff(1 : q);          % vector b
                        c  = r * output_coeff(q + 1);             % constant
                    else
                        y0 = y_in(:, output_index);               % initial guess
                        b  = y_in * coeff(1 : q, output_index);   % vector b
                        c  = r * coeff(q + 1, output_index);      % constants
                    end

                    t_out = t_in + h;
                    y_out = this.nonlinear_solver.solveBC(problem, b, c, y0);
                    if(problem.real_valued)
                        y_out = real(y_out);
                    end                    
                else
                    y_out = [];
                    t_out = [];
                end
            else                
                if(labindex == num_proc) % LEADER NODE
                    
                    y_out = zeros(size(y_in));
                    % -- compute and send information for nonlinear system y = b + c * F(y) ----------------------------
                    for j = 1 : num_solves
                        sol_index = solve_inds(j);
                        y0 = y_in(:, sol_index);                % initial guess
                        b  = y_in * coeff(1 : q, sol_index);    % vector b
                        c  = r * coeff(q + 1, sol_index);       % constants
                        labSend({b, c, y0}, j, 0);
                    end
                    
                    % -- recieve solutions to nonlinear systems y = b + c * F(y) ---------------------------------------
                    for j = 1 : num_solves
                        [data, srcWkrIdx] = labReceive('any', 1);
                        y_out(:, solve_inds(srcWkrIdx)) = data;
                    end
                    
                    % -- set reflected values ----------------------------------------------------------------------
                    for j = 1:length(conj_get)
                        y_out(:,conj_set{j}) = repmat( conj(y_out(:,conj_get(j))), 1, length(conj_set{j}));
                    end
                    
                    t_out = t_in + h;
                    
                else % WORKER NODE
                    
                    sol_index = solve_inds(labindex);
                    
                    % -- recieve information for nonlinear system y = b + c * F(y) -------------------------------------
                    data = labReceive(num_proc, 0);
                    b  = data{1};
                    c  = data{2};
                    y0 = data{3};
                    
                    % -- solve nonlinear system ------------------------------------------------------------------------
                    if(solve_flags(sol_index) == 1) % real valued
                        [y, clean_exit] = this.nonlinear_solver.solveBC(problem, real(b), real(c), real(y0));
                    else % complex valied
                        [y, clean_exit] = this.nonlinear_solver.solveBC(problem, b, c, y0);
                    end
                    
                    % -- send solution ---------------------------------------------------------------------------------
                    labSend(y, num_proc, 1);
                    y_out = [];
                    t_out = [];
                    
                end
                
            end
            
            labBarrier;
            
        end
          
    end
    
end