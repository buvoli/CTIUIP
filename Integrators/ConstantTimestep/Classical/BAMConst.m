% ======================================================================================================================
% BBDFConst : Serial implementation of block Adams-Moulton Methods (BAM) of any order 
%
% -- Constructor Arguments --------------------------------------------------------------------------------------------- 
%
%   options (struct) - contains fields:
%
%       z      -> (vector) vector of nodes z_j.
%       b      -> (vector) vector of nodes b_j.
%       alpha  -> (double) nonlinear solver to be used at each timestep. Defaults to Newton.
%       z_name -> (string) optional string describing nodes
%       nonlinear_solver -> nonlinear solver to be used at each timestep. Defaults to Newton.
%
% ======================================================================================================================

classdef BAMConst < ImplicitIntegratorConst
    
    properties
        graph_line_style = {};
        post_step_rhs_eval = false;
    end
    
    properties(SetAccess = protected)
        name  = '';
        description = '';
        order = [];
        starting_times = [];
        y_coefficients = [];
        f_coefficients = [];
        y_output_coefficients = [];
        f_output_coefficients = [];
        output_index = [];
        z     = []; % quadrature points
        alpha = []; % alpha
        b     = []; % left integration endpoints
        r     = []; % radius
        z_output = 0;
    end
    
    properties(Access = private)
        z_name = '';
    end
    
    methods
        
        function this = BAMConst(options)
            
            if(nargin == 0)
                options = struct();
            end
            
            default_options = { ...
                {'z', 0}        ...
                {'b', 0}        ...
                {'alpha', 1}    ...
                {'z_name', ''}  ...
            };
            options = setDefaultOptions(options, default_options);
            
            this = this@ImplicitIntegratorConst(options);
            this.setParameters(options.z, options.b, options.alpha, options.z_name);
        end
        
        function d = get.description(this)
            [an, ad] = rat(this.alpha);
            if(ad == 1)
                alpha_str = sprintf('alpha = %i', an);
            else
                alpha_str = sprintf('alpha = %i/%i', an, ad);
            end           
            if(isempty(this.z_name))
                d = sprintf('Order %i Block Adams Moulton (%s)', this.order, alpha_str);
            else
                d = sprintf('Order %i Block Adams Moulton (%s, z = %s)', this.order, alpha_str, this.z_name);
            end
        end
        
        function setParameters(this, z, b, alpha, z_name)
            if(nargin < 5)
                z_name = '';
            end
            if(length(z) ~= length(b))
                error('quadrature points z and endpoints b are not of equal length');
            end
            
            this.z = z(:);
            this.z_name = z_name;
            this.b = b(:);
            this.alpha = alpha;
            this.starting_times = double(z(:) / this.alpha);
            this.order = length(z) + 1;
            calculateCoefficients(this);
             % -- adjust stat objects -----------------------------------------------------------------------------------
            this.nonlinear_solver.stats.reset(length(z));
            this.nonlinear_solver.linear_solver.stats.reset(length(z));
            this.step_stats.reset(length(z)); 
            this.rhs_stats.reset(length(z)); 
        end
        
    end
    
    methods (Access = private)
        
        function setName(this)
            if(~isempty(this.z_name))
                this.name = ['BAM', num2str(this.order), '_{', this.z_name, '}'];
            else
                this.name = ['BAM', num2str(this.order)];
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
            bd = double(this.b);
            q  = length(zd);
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
                    z_new  = zd(z_index);
                    b_new  = bd(z_index);
                    z_reflected = abs(zd(search_queue) - conj(z_new)) <= abs_tol;
                    b_reflected = abs(bd(search_queue) - conj(b_new)) <= abs_tol;
                    reflected_z_index = find(z_reflected & b_reflected);
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
                q = length(this.z);
                tau_out = this.z + this.alpha;
                CY = zeros(q, q);
                CF = zeros(q, q + 1);
                for j = 1:q
                    CY(j,:) = double(polyCoefficients(this.z, 1:q, [], this.b(j)));
                    CF(j,:) = double(iPolyCoefficients([this.z; tau_out(j)], [this.b(j) tau_out(j)]));
                end
                this.y_coefficients = CY;
                this.f_coefficients = CF;
                % -- set output coefficients --
                [~, index] = ismember(this.z_output, this.z);
                if(index)
                    this.output_index = index;
                    this.y_output_coefficients = [];  % interpolation formula
                    this.f_output_coefficients = [];  % interpolation formula
                else
                    this.output_index = [];
                    [~, nearest_index] = min(abs(this.z - this.z_output));
                    this.y_output_coefficients = zeros(1, q);
                    this.y_output_coefficients(nearest_index) = 1;
                    zeta_out = [this.z; this.z_output + this.alpha];
                    tau_out  = this.z_output + this.alpha;
                    this.f_output_coefficients = double(iPolyCoefficients(zeta_out, [this.z(nearest_index) tau_out]));
                end
                this.setName();
            end
        end
        
    end
    
    methods (Access = protected)
        
        
        function [step_struct, y_in] = initStepStruct(this, t_in, y_in, problem)
            
            % -- evaluate RHS ------------------------------------------------------------------------------------------
            q   = length(this.z);
            F_in = zeros(size(y_in));
            for j = 1 : q
                F_in(:,j) = problem.RHS(y_in(:,j));                   
            end
            
            % -- real and complex solve indices ------------------------------------------------------------------------
            [solve_flags, conj_get, conj_set] = solveIndices(this, problem);
            
            step_struct = struct(                                        ...
                'F_in',           F_in,                                  ...
                'y_coeff',        transpose(this.y_coefficients),        ...
                'f_coeff',        transpose(this.f_coefficients),        ...
                'y_output_coeff', transpose(this.y_output_coefficients), ...
                'f_output_coeff', transpose(this.f_output_coefficients), ...
                'solve_flags',    solve_flags,                           ...
                'conj_get',       conj_get,                              ...
                'conj_set',       {conj_set}                             ...
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
            y_coeff = step_struct.y_coeff;
            f_coeff = step_struct.f_coeff;
            solve_flags = step_struct.solve_flags;
            
            % -- create empty time record for current step -------------------------------------------------------------
            for j = 1 : q
                this.step_stats.setIndex(j);
                this.step_stats.recordStep(0);
            end
            
            if(final_step)
                
                this.step_stats.setIndex(1);
                output_index = this.output_index;
                [~, output_y0_index] = min(abs(this.z - this.z_output)); % closest initial condition
                
                if(isempty(output_index)) % output node is not a member of node set
                    y0 = y_in(:, output_y0_index);           % initial guess
                    b  = y_in * step_struct.y_output_coeff + step_struct.F_in * (r * step_struct.f_output_coeff(1 : q));
                    c  = r * step_struct.f_output_coeff(q + 1);            % constant
                else % output node is a member of node set
                    y0 = y_in(:, output_index);              % initial guess
                    b  = y_in * y_coeff(:, output_index) + step_struct.F_in * (r * f_coeff(1 : q, output_index));  % vector b
                    c  = r * f_coeff(q + 1, output_index);     % constants
                end
                
                t_out = t_in + h;
                y_out = this.nonlinear_solver.solveBC(problem, b, c, y0);
                if(problem.real_valued)
                    y_out = real(y_out);
                end
                       
            else
                
                y_out = zeros(size(y_in));
                F_out = zeros(size(y_in));
                % -- solve nonlinear systems ---------------------------------------------------------------------------
                for j = 1 : q
                    if(solve_flags(j) > 0)
                        output_start_time = tic;
                        this.step_stats.setIndex(j);
                        this.rhs_stats.setIndex(j);
                        this.nonlinear_solver.stats.setIndex(j);
                        this.nonlinear_solver.linear_solver.stats.setIndex(j);
                        bv = y_in * y_coeff(:, j) + step_struct.F_in * (r * f_coeff(1 : q, j));
                        cf = r * f_coeff(q + 1, j);
                        if(solve_flags(j) == 1) % real valued
                            [y_out(:,j), clean_exit] = this.nonlinear_solver.solveBC(problem, real(bv), real(cf), real(y_in(:,j)));
                        else % complex valued
                            [y_out(:,j), clean_exit] = this.nonlinear_solver.solveBC(problem, bv, cf, y_in(:,j));
                        end
                        % -- test for emergency exit condition ---------------------------------------------------------
                        emergency_exit = ((~clean_exit) || norm(y_out(:,j), Inf) >= Inf || any(isnan(y_out(:,j))));
                        if(emergency_exit)
                            clean_exit
                            y_out = NaN; t_out = t_in;
                            return;
                        end
                        % -- evaluate RHS ------------------------------------------------------------------------------
                        if(this.post_step_rhs_eval)
                            start_rhs_time = tic;
                            F_out(:,j) = problem.RHS(y_out(:,j));
                            this.rhs_stats.recordRHSEval(toc(start_rhs_time));
                        else
                            if(solve_flags(j) == 1) % real valued
                                F_out(:,j) = real(y_out(:,j) - bv)/real(cf);
                            else
                                F_out(:,j) = (y_out(:,j) - bv)/cf;
                            end
                        end
                        this.step_stats.recordSubstep(toc(output_start_time));
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
                    F_out(:,conj_set{j}) = repmat( conj(F_out(:,conj_get(j))), 1, length(conj_set{j}));
                end
                
                % -- update state --------------------------------------------------------------------------------------
                t_out = t_in + h;
                step_struct.F_in = F_out;
                this.step_stats.recordSubstep(toc(step_start_time));
                
            end
             
        end
               
    end
    
end