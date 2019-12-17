function r = reduceResultsBatch(jobName, numRuns)
%REDUCERESULTSBATCH Merge result bundles of a distributed job (on Beagle).
%   R=REDUCERESULTSBATCH(JOBNAME,NUMRUNS,SAVE) loads linear solver result
%   bundles from the output files of runs 1..NUMRUNS of the job JOBNAME,
%   merges them into a single bundle R. and saves R under the output
%   directory.
%
%   See also: REDUCERESULTS, RESULTBUNDLE, LINEARSOLVERS.

r = reduceResults(jobName, numRuns, 1);
end