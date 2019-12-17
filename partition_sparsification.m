clc;
clear;
%add 5% 10% 15% 20% 25% and see the lambda ratio drop 
%change the way to simulate the edge similarity
clear;clc;
fprintf('Read matrix\n');
iteration_times=3;
file = 'ibm32.mtx';
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
%val = abs(cell2mat(B(3)));
val = ones(length(row),1);
n= max(max(row), max(col));
fclose(fp);
toc;
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
%%*************************************************************
[A_S,rate,lambda_ratio,lambda_record]=sparsification(A_G,A_S,25,50,8)
%[A_S,rate,lambda_ratio,lambda_record]=sparsification(A_G,A_S,times,outdegree,lambda_limit)

fprintf('------------------------------------------------------------------\n');
for i=1:length(lambda_ratio)
    fprintf('(ratio=%f); Added off-tree edge %f\n\n', lambda_ratio(i),rate(i));
end

L_S= spdiags(sum(A_S,2),0, n, n)-A_S'+1e-6*spdiags(ones(n,1), 0, n, n);
LGG=sparse(L_G*L_G');
LSS=sparse(L_S*L_S');
eigs(LGG,10,'sm')
eigs(LSS,10,'sm')
%[vg,dg]=eigs(LGG,num_vect, 'sm');
%[vs,ds]=eigs(LSS,num_vect, 'sm');
k=2;
tol=1e-3;
[flagG, flagS] = directed_partitioningFEB_21(LGG, LSS, k, tol);
figure('Name', 'original graph partition', 'NumberTitle', 'off');
GG = digraph(A_G, 'OmitSelfLoops');
groupg1=find(flagG==1);
groupg2=find(flagG==2);
h=plot(GG, 'LineWidth', 2, 'EdgeColor',[ 0 1 1],'Marker', 'o', 'MarkerSize', 6);
%title('original graph partition');
highlight(h,groupg1,'NodeColor','r', 'Marker', 'o', 'MarkerSize', 6);
highlight(h,groupg2,'NodeColor','g',  'Marker', 'o', 'MarkerSize', 6);
savefig('GGGG.fig');

figure('Name', ' reduced graph partition', 'NumberTitle', 'off');
SS = digraph(A_S, 'OmitSelfLoops');
groups1=find(flagS==1);
groups2=find(flagS==2);
h=plot(GG, 'LineWidth', 2, 'EdgeColor',[1 1 1],'Marker', 'o', 'MarkerSize', 6);
%title('sparse graph partition')
highlight(h,groups1,'NodeColor','r',  'Marker', 'o','MarkerSize', 6);
highlight(h,groups2,'NodeColor','g', 'Marker', 'o', 'MarkerSize', 6);
highlight(h,SS,'EdgeColor',[0 1 1],  'LineWidth', 2);
savefig('SSS.fig');
h=plot(GG, 'LineWidth', 2, 'EdgeColor',[1 1 1],'Marker', 'o', 'MarkerSize', 6);
%title('sparse graph partition')
highlight(h,groups1,'NodeColor','g',  'Marker', 'o','MarkerSize', 6);
highlight(h,groups2,'NodeColor','r', 'Marker', 'o', 'MarkerSize', 6);
highlight(h,SS,'EdgeColor',[0 1 1],  'LineWidth', 2);
savefig('SSS_opposite.fig');
groupg1'
groups1'
