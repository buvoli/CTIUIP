% ======================================================================================================================
% PHIEVALUATOR : Abstract class for any method for computing vector products with matrix phi functions
%
% -- Public Properties -------------------------------------------------------------------------------------------------
%
%  1. stats - generic stats object for linear solver 
%
% -- Public Methods ----------------------------------------------------------------------------------------------------
%
%  1. phi   - initialize phi_n(t*A)  
%  2. phib  - initialize phi_n(t*A) b 
%  3. phibs - initialize phi_0(t*A) b(:, 1) + phi_1(t * A) b(:, 2) + ... + phi_{length(b)-1}(t * A) b(:,end)
%  4. solve - phi_0(t*A) b(:, 1) + t phi_1(t * A) b(:, 2) + ... + phi_{length(b)-1}(t * A) b(:,end)
%
% ======================================================================================================================

classdef PhiEvaluator < handle
    
    properties(Abstract = true)
        stats;
    end
    
    methods(Abstract = true)
        [y, exit_flag, residual] = phi(this, n, t, A)       % Initialize phi_n(t*A)
        [y, exit_flag, residual] = phib(this, n, t, A, b)   % Initialize phi_n(t*A) b
        [y, exit_flag, residual] = phibs(this, t, A, b)     % Initialize phi_0(t*A) b(:, 1) + phi_1(t * A) b(:, 2) + ... + phi_{length(b)-1}(t * A) b(:,end)
        [y, exit_flag, residual] = solve(this, t, A, b)     % Initialize phi_0(t*A) b(:, 1) + t phi_1(t * A) b(:, 2) + ... + phi_{length(b)-1}(t * A) b(:,end)
    end
    
end