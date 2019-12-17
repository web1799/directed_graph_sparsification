function [AG,AS]=generated_direct(n)

nn=2*n;% how many egdges to be picked up 
iter=11;
AG=zeros(2*n,2*n);%A adjancy matrix of oringinal  directed grap
AS=zeros(2*n,2*n);%A adjancy matrix of  sparsifier
LG=zeros (2*n,2*n); %G is the oringinal  directed grap
LS=zeros (2*n,2*n); % sparsifier
A_per=zeros (2*n,2*n);
score=zeros(1,2*n);
for i=1:2*n-1
    AG(i,i+1)=1; 
    AS(i,i+1)=1;
end
AG(2*n,1)=1;
AS(2*n,1)=1;
for i=2:n-1
    AG(i,2*n+1-i)=1;
end
% AG=(AG+AG')/2;
% AS=(AS+AS')/2;
%%
 L_G = digraph(AG); % G is the directed graph 
 L_S = digraph(AS); %%%%%% original sparsifier
 figure; plot(L_G)
 figure; plot(L_S)
fp = fopen('directed_test.mtx', 'w');
fprintf(fp, '%d %d %d\n', length(AG), length(AG), nnz(AG));
[aa, bb, cc] = find(AG);
for i = 1:length(aa)
    fprintf(fp, '%d %d %f\n', aa(i), bb(i), cc(i));
end
fclose(fp);
ep = fopen('directed_initial.mtx', 'w');
fprintf(fp, '%d %d %d\n', length(AS), length(AS), nnz(AS));
[a, b, c] = find(AS);

for i = 1:length(a)
    fprintf(fp, '%d %d %f\n', a(i), b(i), c(i));
end
fclose(ep);
end