

p0 = rSol.pressure;
n=nx*ny*nz;
for i=1:numel(W)
p0(n+i) = 0;
end
    switch lsolver
        case 1
            mrstModule add agmg
            solver = AGMGSolverAD('tolerance', tol);
        case 2
            solver = GMRES_ILUSolverAD('tolerance', tol);
        case 3
            solver = ICCGSolverAD('tolerance', tol_p,'maxIterations', 2000,'x0',p0,'W', W);
        case 4
            solver = DICCGSolverAD('tolerance', tol_p,'maxIterations', 2000,'Z',Z,'x0',p0,'W', W);
        otherwise
            solver = BackslashSolverAD();
    end