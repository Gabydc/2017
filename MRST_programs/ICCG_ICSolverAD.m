classdef PCG_ICSolverAD_cn < LinearSolverAD
%Preconditioned CG solver.
%
% SYNOPSIS:
%   solver = PCG_ICSolverAD()
%
% DESCRIPTION:
%   Solve a linearized problem using a CG solver with IC
%   preconditioner.
%
% REQUIRED PARAMETERS:
%   None
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   See class properties.
%
%
% SEE ALSO:
%   BackslashSolverAD, CPRSolverAD, LinearizedProblem, LinearSolverAD

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
    %}
    properties
        % Type: {'nofill'}, 'crout', 'ictp'
        ictype
        % Drop tolerance for elements (default: 0)
        dropTolerance
        % Modified ic type: 'row', 'col', {'off'}
        michol
        % Boolean indicating replacement of zero diagonal entries
        diagcomp 
        % Pivoting threshold for ilupt. Default 1.
        shape 
        %  Determines which triangle is referenced and returned        
           x0
         % Initial condition
         cn 
         %Compute the condition number
   
    end
    methods
        function solver = PCG_ICSolverAD_cn(varargin)
            solver = solver@LinearSolverAD();
            %      type     --- Type of factorization.
%      droptol  --- Drop tolerance when type is 'ict'.
%      michol   --- Indicates whether to perform Modified incomplete Cholesky.
%      diagcomp --- Perform Compensated incomplete Cholesky with the specified
%                   coefficient.
%      shape    --- Determines which triangle is referenced and returned. The default value is 'lower'
%
            solver.ictype               = 'nofill';
            solver.dropTolerance         = 0;
            solver.michol = 'off';
            solver.diagcomp     = 0;
            solver.shape        = 'lower';
            solver.x0 = [];  
            solver.cn = 0;
            
            solver = merge_options(solver, varargin{:});
        end
        
        function [result, report] = solveLinearSystem(solver, A, b)
 
            nel = size(A, 1);
%             if solver.reorderEquations
%                 [A, b] = reorderForILU(A, b);
%             end
%tic
            [L] = ichol(A, solver.getOptsIC());
           % tl=toc
            %pcg(A,b,tol,maxit,M1,M2,x0,varargin)
           % prec = @(x) L'\(L\x);
            
           
          if solver.cn == 0 
% %tic

%               [result,resvec,iter,flag] = pcg_m(A,solver.x0, b,L,  ...
%                 min(solver.maxIterations, nel),solver.tolerance);
%               [result,resvec,iter,flag] = pcg_cn(A,solver.x0, b,L,  ...
%                 min(solver.maxIterations, nel),solver.tolerance);

              [result,flag,relres,iter,resvec] = pcg_17_10_9(A,b,...
                solver.tolerance,min(solver.maxIterations,nel),L,L',solver.x0);
                    iter   
  A(nel+1:end,:)
 b(nel+1:end,:)
x = A\b;
diff = norm(result - x)  
                    
            % [x, error, iter, flag] = pcg_m(A, x, b, M, max_it, tol)
          
%             [result,resvec,iter,flag] = pcg(A, solver.x0, b, M,min(solver.maxIterations, nel), solver.tolerance);
%             relres = max(error);
%             resvec = error;

%             [result,flag,relres,iter,resvec] = pcg_Matlab(A,b,...
%                 solver.tolerance,min(solver.maxIterations, nel),...
%                 L,L',solver.x0);
% %             
            
%tpcg=toc
                     report = struct('CGflag',  flag, ...
                            'residual',   resvec,...
                            'iterations', iter);
            else   
         %Solver computing estimated condition number
%tic
           [result,flag,res,its,resvec,cond] = pcg_cn(A,b,solver.tolerance,min(solver.maxIterations, nel),L,L',solver.x0);          
 % tpcg=toc       
  report = struct('CGflag',  flag, ...
                            'residual',   resvec,...
                            'iterations', its,...
                            'cond', cond);
          end
                       
        end
        
        function opts = getOptsIC(solver)
            opts = struct('type',    solver.ictype, ...
                'droptol', solver.dropTolerance, ...
                'michol',    solver.michol , ...
                'diagcomp',   solver.diagcomp, ...
                'shape',  solver.shape);
        end
        
       function [dx, result, report] = solveLinearProblem(solver, problem, model)
           % Reduce to cell variables before solving
           [dx, result, report]= solver.solveCellReducedLinearProblem(problem, model);
       end
    end
end
