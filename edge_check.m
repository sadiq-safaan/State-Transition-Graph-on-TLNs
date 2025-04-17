function [yn]=edge_check(lc_g,lc_r,nc_g,nc_r,W,theta)
    
    % This function determines whether an edge between two regions produced by the linear boundary/nullcline hyperplane arrangement of a CTLN exists in the state transition graph
    
    %Inputs: linear chamber and nullcline codewords of the receiving region (lc_r,nc_r) and of the candidate edge giving region (lc_g,nc_g), "W" synaptic weight matrix, "theta" external input
    %Output: 1 if edge exists, 0 if edge does not exist
    
    op = optimoptions('linprog','Display','none');
    
    tol=1e-6; %parameter to enforce strict inequality
    solution=[];
    count=0;
    if isempty(lc_g)==1
        lc_g=[];
    end

    if isempty(lc_r)==1
        lc_r=[];
    end

    if isempty(nc_g)==1
        nc_g=[];
    end
    
    if isempty(nc_r)==1
        nc_r=[];
    end
    
    yn=0;

    n=size(W,1);
    
    z=ones(1,n); %optimization function
    lower_bound=zeros(1,n); %positive orthant restriction
    ub=1000*ones(1,n); %dummy upper bound
    
    %Chamber Boundary Constraints
    A1=W;
    A1(lc_g,:)=-A1(lc_g,:);
    b1=-theta*ones(n,1);
    b1(lc_g)=theta;
    
    %Nullcline Constraints
    A2=-eye(n)+W;
    A2(nc_g,:)=-A2(nc_g,:);
    b2=-theta*ones(n,1);
    b2(nc_g)=theta;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Linear Chamber Boundary Face
    if length(nc_g)==length(nc_r)
        
        A2=A2(lc_g,:);
        b2=b2(lc_g);
        
        f=setxor(lc_g,lc_r); %Transition face
        
        f_vec=zeros(1,n);
        
        
        %dxdt dot normal
        for k=1:n
            f_vec(k)=-W(f,k)+sum(W(f,lc_g).*W(lc_g,k)');
        end
        f_c=sum(theta.*W(f,lc_g));
        
        %Restrict to face
        A1_r=A1(setdiff(1:n,f),:);
        b1_r=b1(setdiff(1:n,f));
        Aeq=A1(f,:);
        beq=b1(f);
        
        %Receiver is 'minus' side of face
        if length(lc_g)>length(lc_r)
            A=[A1_r;A2;f_vec];
            b=[b1_r-tol;b2-tol;-f_c-tol]';
            
            solution=linprog(z,A,b,Aeq,beq,lower_bound,ub,op);
            count=count+1;
        end
        
        %Receiver is 'plus' side of face
        if length(lc_g)<length(lc_r)
            A=[A1_r;A2;-f_vec];
            b=[b1_r-tol;b2-tol;f_c-tol]';
            
            solution=linprog(z,A,b,Aeq,beq,lower_bound,ub,op);
            count=count+1;
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Nullcline Face
    if length(lc_g)==length(lc_r)
        f=setxor(nc_g,nc_r); %Transition face
        
        f_vec=zeros(1,n);
        
        %dxdt dot normal
        for k=1:n
            f_vec(k)=-W(f,k)+sum(W(f,lc_g).*W(lc_g,k)');
        end
        f_c=sum(theta.*W(f,lc_g));
        

        %Restrict to face
        A2_r=A2(setdiff(lc_g,f),:);
        b2_r=b2(setdiff(lc_g,f));
        Aeq=A2(f,:);
        beq=b2(f);

        %Receiver is 'minus' side of face
        if length(nc_g)>length(nc_r)
            A=[A1;A2_r;f_vec];
            b=[b1-tol;b2_r-tol;-f_c-tol]';
            
            solution=linprog(z,A,b,Aeq,beq,lower_bound,ub,op);
            count=count+1;
        end
        
        %Receiver is 'plus' side of face
        if length(nc_g)<length(nc_r)
            A=[A1;A2_r;-f_vec];
            b=[b1-tol;b2_r-tol;f_c-tol]';
            
            solution=linprog(z,A,b,Aeq,beq,lower_bound,ub,op);
            count=count+1;
        end
    end
    
    if isempty(solution)==0
        yn=1;
    end
    count;
end