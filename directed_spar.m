%add 5% 10% 15% 20% 25% and see the lambda ratio drop 
%change the way to simulate the edge similarity
clear;clc;
fprintf('Read matrix\n');
iteration_times=3;
file = 'gre_115.mtx';
tic;
fp = fopen(file, 'r');
B = textscan(fp, '%d %d %f ', 'headerlines',1); % read data instead of firt lin
%B = textscan(fp, '%d %d', 'headerlines',1); % read data instead of first lin
row = cell2mat(B(1));
col = cell2mat(B(2));
if(min(row)<1)
	row = row+1;
end
if(min(col)<1)
	col = col+1;
end
val = abs(cell2mat(B(3)));
%val = ones(length(row),1);
n= max(max(row), max(col));
fclose(fp);
toc;
%%

% n = max(max(row), max(col));%number of node 
A_G = sparse(double(row), double(col), double(val), double(n), double(n));
n=length(A_G);%number of node n=nn in square matrix;
A_G=A_G-spdiags(diag(A_G), 0, n, n);
L_G= spdiags(sum(A_G,2),0, n, n)-A_G'+1e-6*spdiags(ones(n,1), 0, n, n);
[row,col,val]=find(A_G);


%% maxmium spanning tree
tic;
A_un=sparse(n,n);A_S=sparse(n,n);A_S_old=sparse(n,n);
A_un=spdiags(1./sum(A_G,2), 0, n, n)*A_G';
A_un=triu((A_un+A_un')/2); 
[row_m,col_m,val_m]=find(A_un);

Graph_modified=graph(row_m,col_m,1./val_m);
T= minspantree(Graph_modified) 
%************************************************
ss= table2array(T.Edges);
spanning_row=ss(:,1);
spanning_col=ss(:,2);
flag2=sparse(double(spanning_row),double(spanning_col),1,n,n);
A_S=flag2.*A_G;
%% Adding outgoing edges for nodes without one
row_og=find(sum(A_S,2)==0); %outgoing edges
if(~isempty(row_og))
    for i=1:length(row_og)
        [max_num,max_loc]=max(A_G(row_og(i),:));
        flag2(row_og(i),max_loc)=1;
    end
    A_S=flag2.*A_G;
end
toc;

% file = 'initial guess.mtx';
% fp = fopen(file, 'r');
% C = textscan(fp, '%d %d %f', 'headerlines', 1); % read data instead of first line
% aa = cell2mat(C(1));
% bb = cell2mat(C(2));
% cc = cell2mat(C(3));
% fclose(fp);
% A_S = sparse(double(aa), double(bb), cc);
%%*************************************************************
[A_S,rate,lambda_ratio,lambda_record]=sparsification(A_G,A_S,50,7,8);
%[A_S,rate,lambda_ratio,lambda_record]=sparsification(A_G,A_S,times,outdegree,lambda_limit)

fprintf('------------------------------------------------------------------\n');
for i=1:length(lambda_ratio)
   fprintf('(ratio=%f); Added off-tree edge %f\n\n', lambda_ratio(i),rate(i));
end
lambda_record



