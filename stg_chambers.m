function [chambers]=stg_chambers(W,theta)
    
    %This function produces the chambers of the linear boundary and nullcline hyperplane arrangement which lie in the positive orthant of a CTLN
    
    %Inputs: W (Synaptic Weight Matrix), epsilon, delta, theta
    %Outputs: Array of feasible regions with nc,lc fields which index the nullclines (nc) and linear chamber boundaries (lc) that the region lies beneath
    
    tol=1e-6; %parameter to enforce strict inequality
    op = optimoptions('linprog','Display','none');

    n=size(W,1);
    chambers=[];
    lower_bound=zeros(1,n);
    f=ones(1,n);
    for i=0:n
        sets=nchoosek(1:n,i);
        for j=1:size(sets,1)
            lin_cham=sets(j,:);
            b1=-theta*ones(n,1);
            b1(lin_cham)=theta;
            for k=0:i
                null_sets=nchoosek(lin_cham,k);
                
                for z=1:size(null_sets,1)
                   null_cham=null_sets(z,:);
                   
                   if length(lin_cham)==1 && k==0
                       null_cham=1:0;
                   end
                   
                   b2=-theta*ones(n,1);
                   b2(null_cham)=theta;
                   b2=b2(lin_cham);
                   b=[b1-tol;b2-tol]';
                   
                   A1=W;
                   A1(lin_cham,:)=-A1(lin_cham,:);
                   
                   A2=-eye(n)+W;
                   A2(null_cham,:)=-A2(null_cham,:);
                   A2=A2(lin_cham,:);
                   
                   A=[A1;A2];
                   
                   solu=linprog(f,A,b,[],[],lower_bound,[],op);
                   if isempty(solu)==0
                       cham=struct('nc',null_cham,'lc',lin_cham);
                       chambers=[chambers, cham];
                   end
                end
            end
        end
    end
end