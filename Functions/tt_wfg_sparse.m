%% this wfg_sparse version is to reduce recycling of code
%
%%
% X=sptensor(randn(8,8,8,8));
% S=size(X);N=ndims(X);
% r=16*ones(1,numel(S)-1); 
% missing_rate=0.9;              % Missing rate: [0 ~ 1]
% Omega = randperm(prod(S)); 
% Omega = Omega(1:round((1-missing_rate)*prod(S)));
% W = zeros(S); 
% W(Omega) = 1;
% W=sptensor(W);
% G=initcoreten(S,r);
% for i=1:N
%     G{i}=12*G{i};
% end
% G0_v=Gm2Gv(Gt2Gm(G));

%%
function [f,G_gradient]=tt_wfg_sparse(X,W,G0_v,r)

% X and W is the sparse data type by tensor toolbox
% G0_v is the vector form of G

% Set observed entries to a vector form, ignore all the missing entries
% Use W to size Zvals because there are maybe 0 values in X
%% Get a vector  named Yvals only contains observed entries
X=sptensor(X);W=sptensor(W);
Y=sptensor(double(W).*double(X));
[~,loc] = ismember(Y.subs,W.subs,'rows');
Yvals=zeros(size(W.vals));
Yvals(loc) = Y.vals;

Q=numel(Yvals);
N=ndims(X);
S=size(X);
Wsubs = W.subs;
G=Gm2Gt(Gv2Gm(G0_v,S,r),r); % transit G from vector to cell mode tensor

%% Get the vector named Zvals only contains observed entries from core tensors
% get g_matrix of size NvalsxN containing every slice factor g accroding to
% observed entries

% get Zvals
Zvals=zeros(Q,1);
for q=1:Q
    Z=reshape(G{1}(:,Wsubs(q,1),:),[1,size(G{1},3)]);
    for n=2:N
        Z=Z*squeeze(G{n}(:,Wsubs(q,n),:));
    end
    Zvals(q)=Z;
end

%% Calculate function value & gradient
% get function value
f=0.5*sum(Yvals.^2)-Yvals'*Zvals+0.5*sum(Zvals.^2);

% get gradiet_matrix the same size as g_matrix containing every factor
% gradient according to observed entries
Tvals=Zvals-Yvals;

% accumulate gradient of every slice g with same index
G_gradient=zeros(size(G0_v));
G_gradient=Gm2Gt(Gv2Gm(G_gradient,S,r),r); 

for q=1:Q
    for n=1:N
        if n==1
             g_nl=1;
           g_nr=squeeze(G{2}(:,Wsubs(q,2),:));
           for i=3:N
              g_nr=g_nr*squeeze(G{i}(:,Wsubs(q,i),:));
           end    
        G_gradient{1}=squeeze(G_gradient{1});
        G_gradient{1}(Wsubs(q,1),:)=G_gradient{1}(Wsubs(q,1),:)+Tvals(q)*g_nl'*g_nr';
        G_gradient{1}=reshape(G_gradient{1},[1,size(G_gradient{1},1),size(G_gradient{1},2)]);
        
        elseif n==2
            g_nl=reshape(G{1}(:,Wsubs(q,1),:),[1,size(G{1},3)]);
           g_nr=squeeze(G{3}(:,Wsubs(q,3),:));
           for i=4:N
               g_nr=g_nr*squeeze(G{i}(:,Wsubs(q,i),:));
           end
           g=Tvals(q)*g_nl'*g_nr';
            G_gradient{2}=permute(G_gradient{2},[2,1,3]);
           G_gradient{2}(Wsubs(q,2),:,:)=G_gradient{n}(Wsubs(q,2),:,:)+reshape(g,[1,size(g,1),size(g,2)]);
             G_gradient{2}=permute(G_gradient{2},[2,1,3]);
           
        elseif n==N
            g_nr=1;
            g_nl=reshape(G{1}(:,Wsubs(q,1),:),[1,size(G{1},3)]);
           for i=2:N-1
               g_nl=g_nl*squeeze(G{i}(:,Wsubs(q,i),:));
           end
             G_gradient{N}(:,Wsubs(q,N))=G_gradient{N}(:,Wsubs(q,N))+Tvals(q)*g_nl'*g_nr';
        elseif n==N-1
            g_nr=squeeze(G{N}(:,Wsubs(q,N),:));
            g_nl=reshape(G{1}(:,Wsubs(q,1),:),[1,size(G{1},3)]);
            for i=2:N-2
                g_nl=g_nl*squeeze(G{i}(:,Wsubs(q,i),:));
            end
           g=Tvals(q)*g_nl'*g_nr';
            G_gradient{N-1}=permute(G_gradient{N-1},[2,1,3]);
            G_gradient{N-1}(Wsubs(q,N-1),:,:)=G_gradient{n}(Wsubs(q,N-1),:,:)+reshape(g,[1,size(g,1),size(g,2)]);
            G_gradient{N-1}=permute(G_gradient{N-1},[2,1,3]);
    else % n=2:N-2
                    g_nl=reshape(G{1}(:,Wsubs(q,1),:),[1,size(G{1},3)]);
            for i=2:n-1
                g_nl=g_nl*squeeze(G{i}(:,Wsubs(q,i),:));
            end
            g_nr=squeeze(G{n+1}(:,Wsubs(q,n+1),:));
            for i=n+2:N
                g_nr=g_nr*squeeze(G{i}(:,Wsubs(q,i),:));
            end
        g=Tvals(q)*g_nl'*g_nr';
        G_gradient{n}=permute(G_gradient{n},[2,1,3]);
        G_gradient{n}(Wsubs(q,n),:,:)=G_gradient{n}(Wsubs(q,n),:,:)+reshape(g,[1,size(g,1),size(g,2)]);
        G_gradient{n}=permute(G_gradient{n},[2,1,3]);
         end
    end
end
G_gradient=Gm2Gv(Gt2Gm(G_gradient));
end