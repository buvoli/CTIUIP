% ======================================================================================================================
% EXPONENTIALINTEGRATORCONST : The base class for exponential integrators
%
% -- Public Properties -------------------------------------------------------------------------------------------------
%
% 1. phi_evaluator - Class for initializing phi functions. Must have NonlinearSolver as superclass
%
% ======================================================================================================================

classdef ExponentialIntegratorConst < IntegratorConst
    
    properties
        phi_evaluator
    end
    
    methods
        
        function this = ExponentialIntegratorConst(options)
            % -- set basic integrator properties -----------------------------------------------------------------------
            this@IntegratorConst(options);
            % -- set local properties ----------------------------------------------------------------------------------
            default_field_value_pairs = {
                {'phi_evaluator', KIOPS()}
            };
            options = setDefaultOptions(options, default_field_value_pairs);
            % -- set props ---------------------------------------------------------------------------------------------
            props = {'phi_evaluator'};
            this.setEmptyClassProps(props, options);
        end
        
    end
    
end