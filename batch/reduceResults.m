function r = reduceResults(jobName, numRuns, varargin)
%REDUCERESULTS Merge result bundles of a distributed job.
%   R=REDUCERESULTS(JOBNAME,NUMRUNS,SAVE) loads linear solver result
%   bundles from the output files of runs 1..NUMRUNS of the job JOBNAME,
%   and merges them into a single bundle R. If SAVE=1, the merged file is
%   saved under the output directory.
%
%   See also: RESULTBUNDLE, LINEARSOLVERS.

if (numel(varargin) < 1)
    saveReduced = false;
else
    saveReduced = varargin{1};
end

r = load(getOutputFile(jobName, 1));
r = r.result;
for runId = 2:numRuns
    s = load(getOutputFile(jobName, runId));
    s = s.result;
    r.appendAll(s);
end

% Save reduced results
if (saveReduced)
    mergedFile = getOutputFile(jobName, []);
    create_dir(mergedFile, 'file');
    save(mergedFile ,'r');
end

%-----------------------------------------------------
function file = getOutputFile(jobName, runId)
% Return the output file corresponding to run id # runId. If runId = [],
% returns the output file corresponding to the merged run.
global GLOBAL_VARS;
if (isempty(runId))
    file = sprintf('%s/%s/cycle_results.mat', GLOBAL_VARS.out_dir, jobName);
else
    file = sprintf('%s/%s_%03d/cycle_results.mat', GLOBAL_VARS.out_dir, jobName, runId);
end
