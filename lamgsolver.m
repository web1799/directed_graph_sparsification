% solve equation Ax=b using LAMG

function x=lamgsolver(A, b)
	lamg    = Solvers.newSolver('lamg', 'randomSeed', 1);
    tStart = tic;
    setup = lamg.setup('sdd', A);
	tSetup = toc(tStart);
    setRandomSeed(now);
	% Turn on debugging printouts during the run
	core.logging.Logger.setLevel('lin.api.AcfComputer', core.logging.LogLevel.DEBUG);
	tStart = tic;
	%[x, ~, ~, details] = lamg.solve(setup, b, 'errorReductionTol', 1e-12);
    	[x, ~, ~, details] = lamg.solve(setup, b, 'finalErrorNorm', 1e-3);
	tSolve = toc(tStart);
	%disp(setup.setupLaplacian);
	tMvm    = mvmTime(A, 5);
	nnz     = numel(nonzeros(A));

	%fprintf('MVM time [sec]       elapsed %.3f, per nonzero %.2e\n', tMvm, tMvm / nnz);
    %fprintf('Setup time [sec]   elapsed %.3f, Solve time [sec]     elapsed %.3f',tSetup, tSolve);
	%fprintf('|A*x-b|/|b|: %.2e: ', norm(A*x-b)/norm(b));

    if(norm(A*x-b)/norm(b) > 1e-3)
        fprintf('|A*x-b|/|b|: %.2e: Not reaching desired solution for LAMG\n', norm(A*x-b)/norm(b));
        %x=A\b;
        %fprintf('|A*x-b|/|b|: %.2e\n', norm(A*x-b)/norm(b));
    end
	if (isfield(details, 'acf'))
    	%fprintf('Convergence factor   %.3f\n', details.acf);
	end
end
