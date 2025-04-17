function [stg]=stg_edges(c,geo_adj,W,theta)
    
    % This function produces the state transition graph for the hyperplane arrangement given by the linear chamber boundaries and nullclines
    
    %Inputs: "geo_adj" Geometric adjacency matrix (i.e. undirected graph of which chambers are adjacent to one another in phase space), "c" list of feasible chambers, "W" synaptic weight matrix, "theta" external input value
    %Output: Adjacency Matrix of state transition graph
    
    k=length(c);
    stg=zeros(k);
    
    chambers=c;
    
    for i=1:k
        for j=1:k
            if geo_adj(i,j)==1
                
                %"Edge Receiving" Chamber
                lc_r=c(i).lc;
                nc_r=c(i).nc;
                
                %"Candidate Edge Giving" Chamber
                lc_g=c(j).lc;
                nc_g=c(j).nc;
                
                %Checking if it is possible to move from j to i
                stg(i,j)=edge_check(lc_g,lc_r,nc_g,nc_r,W,theta);
            end
        end
    end
end