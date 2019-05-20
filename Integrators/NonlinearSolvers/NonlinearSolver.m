% ======================================================================================================================
% NONLINEARSOLVER : Abstract class for any linear solver
%
% -- Public Properties -------------------------------------------------------------------------------------------------
%
%  1. stats - generic stats object for nonlinear solver 
%
% -- Public Methods ----------------------------------------------------------------------------------------------------
%
%  1. solve   - solve system of form F(x) = 0  
%  2. solveBC - solve system of form Solve system of form x = b + c * F(x) 
%
% ======================================================================================================================

classdef NonlinearSolver < handle
    
    properties(Abstract = true)
        stats;
    end
    
    methods(Abstract = true)
        [y, exit_flag, residual] = solve(this)       % Solve system F(x) = 0
        [y, exit_flag, residual] = solveBC(this)     % Solve system of form x = b + c * F(x) 
    end
    
end

