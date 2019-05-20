classdef spmd_burnin
    %SPMD_BURNIN Summary of this class goes here
    %   Detailed explanation goes here
    
    
    methods
        function obj = spmd_burnin(problem)
            %SPMD_BURNIN Construct an instance of this class
            %   Detailed explanation goes here
            obj.burnin(problem);
        end
        
        function output = burnin(~, problem)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            spmd(2)
                if(labindex == 1)
                    val = problem.initial_condition;
                else
                    val = 2;
                end
            end
            
            output = val{2};
            
        end
    end
end

