fprintf('Read matrix\n');
file = 'LL.mtx';
fp = fopen(file, 'r');
B = textscan(fp, '%d %d %f', 'headerlines', 1); % read data instead of first line
row = cell2mat(B(1));
col = cell2mat(B(2));
val = cell2mat(B(3));
fclose(fp);
nz = length(row); %number of edge
A = sparse(double(row), double(col), val);
n = size(A,1);
b = rand(n, 1);

[x, setup] = run_solver('sdd', 'lamg', A, b);
