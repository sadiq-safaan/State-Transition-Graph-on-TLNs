function [adj_matrix]=stg_adjacency(c)
    %This function produces a geometric adjacency matrix of regions for a feasible linear chamber/nullcline hyperplane arrangement of a CTLN
    
    %Input: Hyperplane arrangement chamber partition
    %Output: Geometric adjacency matrix
    
    k=length(c);
    adj_matrix=zeros(k);
    for i=1:k
        for j=(i+1):k
            lc1=c(i).lc;
            lc2=c(j).lc;
            nc1=c(i).nc;
            nc2=c(j).nc;
            diff=[setxor(lc1,lc2) setxor(nc1,nc2)];
            if length(diff)==1
                adj_matrix(i,j)=1;
                adj_matrix(j,i)=1;
            end
        end
    end
end