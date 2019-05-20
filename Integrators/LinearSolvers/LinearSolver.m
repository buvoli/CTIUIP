% ======================================================================================================================
% LINEARSOLVER : Abstract class for any linear solver
%
% -- Public Properties -------------------------------------------------------------------------------------------------
%
%  1. stats - generic stats object for linear solver 
%
% -- Public Methods ----------------------------------------------------------------------------------------------------
%
%  1. solve   - Solve system A * x = b 
%  2. solveBC - Solve system x = b + c * A * x 
%
% ======================================================================================================================

classdef LinearSolver < handle
    
    properties(Abstract = true)
        stats;
    end
    
    methods(Abstract = true)
        [y, flag, residual] = solve(this);       % Solve system A * x = b
        [y, flag, residual] = solveBC(this);     % Solve system of form x = b + c * A * x 
    end
    
end