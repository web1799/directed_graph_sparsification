function mccCompile(f)
%MCCCOMPILE Compile to a Beagle executable.
%   Compile a LAMG entrypoint function to an executable on the Beagle Cray
%   XE6 machine.
%
%   See also: MCC.

config;

% Change to the directory under which we keep batch runnables
currentDir = pwd;
global GLOBAL_VARS
homeDir = [GLOBAL_VARS.home_dir '/code'];
eval(sprintf('cd %s/../bin', homeDir));

% MCC-compile
cmd = sprintf('mcc -v -a %s -m %s.m -R -nojvm -R -singleCompThread -R -nodisplay -o %s', homeDir, f, f);
fprintf('Running %s\n', cmd);
eval(cmd);

% Go back to original directory
eval(sprintf('cd %s', currentDir));
end