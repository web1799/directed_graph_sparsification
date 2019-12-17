function p = rmpathsub(dir)
%RMPATHSUB Remove directory and sub-directories from search path.
%   P=RMPATHSUB(DIR) removes the specified directory and all its
%   sub-directories found on the matlabpath and returns the updated path
%   string. This is case-sensitive, so RMPATHSUB('C:\dir') will not remove
%   the directory 'c:\dir' from the matlabpath on Windows matchines.
%
%   See also: PATH, RMPATH.

warning('off', 'MATLAB:rmpath:DirNotFound');

% Tokenize MATLAB path to a cell array
a = regexp(pathdef,['([^' pathsep '])*'], 'tokens');
a = [a{:}];

% Identify directories on the matlabpath that start with the string DIR
s = regexptranslate('escape', dir);
r = sprintf('%s(.+)', s);
b = a(cellfun(@(x)(~isempty(regexp(x, r, 'start'))), a));

% Remove those directories from the matlabpath
if (~isempty(b))
    rmpath(b{:});
end

% Make sure the current directory is removed. This doesn't always seem to
% happen.
rmpath(dir);

p = path;

warning('on', 'MATLAB:rmpath:DirNotFound');
end
