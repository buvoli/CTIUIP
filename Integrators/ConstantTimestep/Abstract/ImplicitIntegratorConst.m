% ======================================================================================================================
% IMPLICITINTEGRATORCONST : The base class for all implicit integrators. 
%
% -- Public Properties -------------------------------------------------------------------------------------------------
%
% 1. nonlinear solver - Class for solving nonlinear systems. Must have NonlinearSolver as superclass
%
% ======================================================================================================================

classdef ImplicitIntegratorConst < IntegratorConst
    
    properties
        nonlinear_solver
    end
    
    methods
        
        function this = ImplicitIntegratorConst(options)
            % -- set basic integrator properties -----------------------------------------------------------------------
            this@IntegratorConst(options);
            % -- set local properties ----------------------------------------------------------------------------------
            default_field_value_pairs = {
                {'nonlinear_solver', Newton()}
            };
            options = setDefaultOptions(options, default_field_value_pairs);
            % -- set props ---------------------------------------------------------------------------------------------
            props = {'nonlinear_solver'};
            this.setEmptyClassProps(props, options);
        end
        
    end
    
end