function result = linearSolvers(jobName, numRuns, runId, action, keyRegexp, minEdges, maxEdges, solvers)
%LINEARSOLVERS LAMG projec graph Laplacian linear solver batch runs - main entrypoint.
%   RESULTBUNDLE=LINEARSOLVERS(JOBNAME,NUMRUNS,RUNID,'run',KEYREGEXP,
%   MINEDGES,MAXEDGES,SOLVERS) runs the list of specified solvers in the
%   SOLVERS cell array on all graph instances in the input directory data
%   index GLOBAL_VARS.DATA_DIR/lap_index.mat whose key matches the regular
%   expression KEYREGEXP and have at least MINEDGES edges and at most
%   MAXEDGES edges. Use NUMRUNS>1 for NUMRUNS parallel runs; test instances
%   are divided roughly equally between run ids RUNID=1..NUMRUNS. This is a
%   function wrapper of Solvers.linearSolvers(). JOBNAME is a
%   uniquely-identifying job name common to all NUMRUNS runs.
%
%   RESULTBUNDLE holds run statistics. The results are also saved under the
%   GLOBAL_VARS.OUT_DIR/<current_date>_<RUNID> directory (the RUNID suffix
%   is omitted if NUMRUNS=1). PRINTER is the printer used to print the
%   result table.
%
%   LINEARSOLVERS([],[],'download',[],MINEDGES,MAXEDGES) downloads instances
%   to the 'lap' data directory instead.
%
%   LINEARSOLVERS(JOBNAME,NUMRUNS,[],'reduce',[],[],[]) merges results from
%   parallel runs and saves the merged output under the output directory.
%
%   See also: SOLVERS, OPTIONS, PLOTRESULTBUNDLE.

%solvers
if (isdeployed)
    % For UChicago's Beagle Crey XE6 deployment
    %str2cell = @(s)(cellfun(@(x)(x(2:end-1)), regexp(s(2:end-1), '\s*,\s*', 'split'), 'UniformOutput', false));
    %str2cell = @(s)(regexp(s(2:end-1),':','split'));
    str2cell = @(s)(regexp(s,':','split'));
    numRuns  = str2double(numRuns);
    runId    = str2double(runId);
    minEdges = str2double(minEdges);
    maxEdges = str2double(maxEdges);
    solvers  = str2cell(solvers);
end
solvers

switch (action)
    case 'run'
        result = Solvers.runSolvers('dieOnException', true, ...
            'numRuns', numRuns, 'runId', runId, ...
            'keyRegexp', keyRegexp, 'minEdges', minEdges, 'maxEdges', maxEdges, ...
            'solvers', solvers, ...
            'outputDir', jobName, 'outputFile', 'solvers.html', ...
            'format', 'html', 'save', true, ...
            'svnCommit', false, 'sendMail', false);
    case 'download',
        Problems.downloadCollection('lap', minEdges, maxEdges);
    case 'reduce',
        % Merge results from parallel runs and save final output
        result = reduceResults(jobName, numRuns, true);
    otherwise,
        error('Unknown action ''%s''', action);
end
