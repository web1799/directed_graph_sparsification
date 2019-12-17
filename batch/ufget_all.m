%UFGET_ALL Download the entire UF collection.
%   This example script gets the index file of the UF sparse matrix
%   collection, and then downloads in all matrices, in increasing order of
%   number of rows in the matrix.
%
%   Example:
%       type ufget_example ; % to see an example of how to use ufget
%
%   See also ufget, ufweb, ufgrep.

index = ufget ;
f = 1:numel(index.nrows);
%f = find(strcmp(index.Group, 'SNAP'))'; f = find (index.numerical_symmetry
%== 1 & ~index.isBinary) ;
[y, j] = sort (index.nrows (f)) ;
f = f (j) ;

j = 0;
n = numel(f);
for i = f %(1:end)
    j = j+1;
    fprintf ('[%d/%d] Loading %s%s%s, please wait ...\n', ...
        j, n, index.Group {i}, filesep, index.Name {i}) ;
    ufget (i,index,'download') ;
end

