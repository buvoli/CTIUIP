% ======================================================================================================================
% PROBLEM : The base class for all initial value problems.
%
% -- Public Properties -------------------------------------------------------------------------------------------------
%
% 1. tspan  : (2x1 vector) integration window
%
% 2. params : (struct) struct storing any additional parameters associated with the IVP. 
%                    
% -- Private Properties ------------------------------------------------------------------------------------------------
%
% 1. dimension : (integer) dimension
%
% 2. real_valued : (boolean) true if IVP is stricly real_valued, false if not.
%
% 3. initial_condition : (vector) initial condition for IVP.
%
% 4. reference_solution : (string) automaticallt genereated final solution for IVP.
%
% 5. name : (string) name of problem.
%
% 6. description : (string) description of problem.
%
% NOTE: to learn about Copyable classes see:
% - https://www.mathworks.com/help/matlab/ref/matlab.mixin.copyable-class.html
% ======================================================================================================================


classdef Problem < matlab.mixin.Copyable
    
	properties (Abstract = true, SetObservable)
        tspan;
        params;        
    end
    
    properties(SetAccess = protected)
    	reference_solution = [];
    end
    
	properties(Abstract, SetAccess = protected)
    	dimension;
        name;
        description;
        initial_condition;
        real_valued;
    end

    methods(Abstract = true, Access = protected)
        reset(this);
        setDimension(this);
    end
    
    methods
        
        function this = Problem()
            addlistener(this,'tspan','PostSet',@this.ObservablesHaveChanged);
            addlistener(this,'params','PostSet',@this.ObservablesHaveChanged);
            this.setDimension();
            this.reset();
        end
        
        function flag = has(this, name)
            if(ismethod(this, name))
                flag = true;
            elseif(isprop(this, name) && ~isempty(this.(name)))
                flag = true;
            else
                flag = false;
            end
        end
        
        function sol = get.reference_solution(this)
            if(isempty(this.reference_solution))
                this.setRefSolution();
            end            
            sol = this.reference_solution;
        end
        
    end
    
    methods(Access = protected)
    
        function ObservablesHaveChanged(this, varargin)
            this.reference_solution = [];
            this.setDimension();
            this.reset();
        end
        
        function setRefSolution(this)
            if(isfield(this.params, 'ref_solve_tol'))
                ref_solve_tol = this.params.ref_solve_tol;
            else
                ref_solve_tol = 3e-14;
            end
            options = odeset('Jacobian', @(varargin) this.J(varargin{2:end}), 'RelTol', ref_solve_tol, 'AbsTol', ref_solve_tol);
            timeSpan = [this.tspan(1), (this.tspan(1) + this.tspan(2))/2, this.tspan(2)];
            [~,exsol] = ode15s(@(varargin) this.RHS(varargin{2:end}), timeSpan, this.initial_condition, options);
            this.reference_solution = transpose(exsol(end,:));
        end
    
    end
    
end