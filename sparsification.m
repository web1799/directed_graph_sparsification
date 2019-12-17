function[A_S,rate,lambda_ratio,lambda_record]=sparsification(A_G,A_S,times,outdegree,lambda_limit)
delta=0.15;
n=length(A_G);
number=int16(n*0.05);
iteration_times=2;
nnz_A_Spanning=nnz(A_S);
nnz_A_G=nnz(A_G); %number of edge
nnz_off_edge=nnz_A_G-nnz(A_S);
[aa,bb,cc]=find(A_S);
flag2=sparse(double(aa),double(bb),1,n,n);
[row,col,val]=find(A_G);
flag=sparse(double(row),double(col),-1,double(n), double(n));

L_G= spdiags(sum(A_G,2), 0, n, n)-A_G'+1e-6*spdiags(ones(n,1), 0, n, n); % directed laplacian
L_S= spdiags(sum(A_S,2),0, n, n)-A_S'+1e-6*spdiags(ones(n,1), 0, n, n);
%[ht,h0,lambda_new]= Ht(L_G, L_S, iteration_times);   
[ht,h0,lambda_new]= Ht_nov7(L_G, L_S, iteration_times);   
lambda_inital=lambda_new;
A_S_old=[];
step=1;
while((lambda_new>lambda_limit)&(step<times))

    lambda_old =lambda_new;
    A_S_old=A_S; 
    [flagoff_row,flagoff_col]=find((flag+flag2)==-1);
    score=zeros(length(flagoff_row),1);
    top=zeros(length(flagoff_row),1);
    top_row=zeros(length(flagoff_row),1);
%*************************************the similarity of off tree edges
    s = Score(A_G, L_S, flagoff_row, flagoff_col, ht);
    score = s*lambda_new;   
    %%
    if 0
        for i=1:length(flagoff_row)  %edge index    
            p =flagoff_row(i);
            q = flagoff_col(i);                                             %%the score of each edge 
            lsep = L_S(:, p); 
            score(i)=(ht(p)-ht(q))*lsep'*ht+ht'*lsep*(ht(p)-ht(q)); 
            score(i)=score(i)*A_G(p,q);%+A_G(p,q)^2*(ht(p)-ht(q))^2;  
        end
        
    end
     %%

%%*************************add the first important edges into sparsifier

    [top,top_row]=sort(score,'descend');
    n_positive=int16(length(find(score>0)));
    n_ch=min(n_positive, number);
    if 1
        score_similarity=zeros(number,1);   
        n_positive=min(number,n_positive);
        for i=1:n_positive
            p = flagoff_row(top_row(i));
            q = flagoff_col(top_row(i));                                             %%the score of each edge
            temps=zeros(3,1);
            lsep = L_S(:, p);
            for j=1:5
                lsep_ht=lsep'*h0(:,j);
                epq_ht=(h0(p,j)-h0(q,j));
                score_similarity(i,j)=2*epq_ht*lsep_ht+(ht(p)-ht(q))^2;
                score_similarity(i,j)=lambda_new*score_similarity(i,j);
            end
        end
    end
	if n_ch>=1;
    flag2(flagoff_row(top_row(1)),flagoff_col(top_row(1)))=1;
    for i=2:n_ch
        temp=zeros(i-1,1);
        for j=1:i-1
           temp(j)=norm(score_similarity(i,:)-score_similarity(j,:));
           temp(j)=temp(j)/norm(score_similarity(i,:));
        end
        outd=length(find(flag2(flagoff_row(top_row(i)),:)));% out dgrees number
        if ((outd<outdegree)&(max(temp)<delta))
       
             flag2(flagoff_row(top_row(i)),flagoff_col(top_row(i)))=1;
		
	
        end
	
    end
end

   
    A_S=flag2.*A_G;
    L_S= spdiags(sum(A_S,2), 0, n, n)-A_S'+1e-6*spdiags(ones(n,1), 0, n, n);
    %[ht,h0,lambda_new]= Ht(L_G, L_S, iteration_times);   
    [ht,h0,lambda_new]= Ht_nov7(L_G, L_S, iteration_times);   


    if lambda_new>lambda_old;
       A_S=A_S_old;
       L_S= spdiags(sum(A_S,2), 0, n, n)-A_S'+1e-6*spdiags(ones(n,1), 0, n, n);
       [ht,h0,lambda_new123]= Ht_nov7(L_G, L_S, iteration_times);   
       [aa,bb,cc]=find(A_S);
	   [lambda_new, lambda_new123]
       flag2=sparse(double(aa),double(bb),1,n,n);
    end
 
%%***************************************************************    
    lambda_record(step)=lambda_new;
    lambda_ratio(step)= lambda_inital/lambda_new;
    rate(step)=nnz(flag2)/nnz_A_G;
    step=step+1;
end

end
