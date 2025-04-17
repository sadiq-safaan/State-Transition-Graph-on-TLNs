%Script producing the state transtion graph for the hyperplane arrangement given by lin. chamber/nullcline partition of a CTLN.  Chambers are labelled in STG in the form (NC,LC) where NC indexes the nullclines the region lies beneath and LC the linear chamber the region lies in.

%Safaan Sadiq October 2023

%Parameters
delta=0.5;
eps=0.25;
theta=1;

%Adjacency Matrix of CTLN
%a_mat=[0 0 1 0; 1 0 1 0; 1 0 0 0; 1 1 1 0];
%a_mat=[0 1 1; 1 0 1; 1 1 0];
a_mat=[0 0 1; 0 0 0; 0 0 0];
%a_mat=[0 1 0; 1 0 0; 1 1 0];

%a_mat=[0 1 0 1 1;0 0 1 0 0;0 0 0 1 1; 0 0 0 0 1; 0 0 0 0 0];

%a_mat=[0 0 0 0 1; 1 0 1 0 0; 0 0 0 0 1; 1 0 1 0 0; 0 1 0 1 0];

n=size(a_mat,1);

%Synaptic Weight Matrix
W=(-1-delta)*ones(n)+(1+delta)*eye(n)+(delta+eps)*a_mat;

%Calculate all feasible chambers in positive orthant
chambers=stg_chambers(W,theta);

%Find geometrically adjacent chambers
geo_adj=stg_adjacency(chambers);

%Calculate state transition graph
stg=stg_edges(chambers,geo_adj,W,theta);


labels=string(1:length(chambers));

for i=1:length(chambers)
    nc=chambers(i).nc;
    lc=chambers(i).lc;
    
    if isempty(nc)==1
        nc=0;
    end
    
    if isempty(lc)==1
        lc=0;
    end
    
    %labels(i)=join(string(nc))+","+join(string(lc));
    
    %Create Region Labels for STG
    l_lc="";
    l_nc="";
    for j=1:n
        if sum(ismember(lc,j))>0.5
            l_lc=l_lc+"+";
        else
            l_lc=l_lc+"-";
        end
        
        if sum(ismember(nc,j))>0.5
            l_nc=l_nc+"+";
        else
            l_nc=l_nc+"-";
        end
    end
    
    labels(i)=l_nc+","+l_lc;
    
end

%Find Strongly Connected Components
G=digraph(stg');
bins=conncomp(G);

%Plot State Transition Graph
figure(1)
subplot(2,1,1)
p=plot(G,'nodeLabel',labels);
p.NodeCData=bins;
colormap(hsv(max(bins)))
rotate(p,[0 0 1],90)


G2=condensation(G);

%Plot SCC DAG
subplot(2,1,2)
p2=plot(G2);
p2.NodeCData=1:max(bins);
colormap(hsv(max(bins)))

%plot(graph(geo_adj),'nodeLabel',labels)