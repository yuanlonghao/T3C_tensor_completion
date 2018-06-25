
function G_opted=tt_sgd(varargin)
%% variable set
% The 1st and 2nd input must be X and W
if (nargin == 0||nargin == 1)
    G_opted=0;
    fprintf('Not enough input.\n');
else
%% 
ip = inputParser;
ip.addParamValue('Tol', 1e-20, @isscalar);
ip.addParamValue('Rank', 16*ones(1,ndims(varargin{1})-1), @ismatrix);
ip.addParamValue('MaxIter', 5e5, @isscalar);
%ip.addParamValue('Alg', 'TT', @ischar);
ip.addParamValue('lr', 5e-4, @isscalar);
ip.parse(varargin{3:end});

Tol = ip.Results.Tol;
Rank = ip.Results.Rank;
lr= ip.Results.lr;
MaxIter=ip.Results.MaxIter;

X=varargin{1};
W=varargin{2};W=sptensor(W);
S=size(X);
r=Rank(:);
N=numel(S);
G=initcoreten(S,r);
G0_v=1e-1*Gm2Gv(Gt2Gm(G));
G=Gm2Gt(Gv2Gm(G0_v,S,r),r);
grad=cell(1,N);

for i=1:N
            if i==1
            grad{1}=zeros(1,S(i),r(i));
           elseif i==N
            grad{N}=zeros(r(i-1),S(i),1);
            else
            grad{i}=zeros(r(i-1),S(i),r(i));
            end
end
            
%% TT-SGD

        %         rng('default');
        sgdmethod = 'Adam';
        lrate =1e-3;    % learning rate
        %lr=[1e-3,1e-3,1e-3];
        beta1 = 0.9;
        beta2 = 0.999;
        
        momentum = 0.5;  % forget rate 
        mu = 0.1;        % constant for diagnoal hessian matrix   
        m_t = cell(1,N);
        v_t = cell(1,N);
        for i=1:N
            m_t{i} = zeros(size(G{i}));
            v_t{i} = zeros(size(G{i}));
        end
        it=0;
        Wsubs=W.subs;
         index=randperm(size(Wsubs,1));
         Q=numel(index);
        while  it<MaxIter
            it = it + 1;    
            
       % sample one tensor entry
       q=round(Q*rand);
       if mod(q,1)~=0 || q<=0
           q=1;
       end
        I=Wsubs(q,:);
       % compute z-y
       z=reshape(G{1}(:,I(1),:),[1,numel(G{1}(:,I(1),:))]);
       for n=2:N-1
           z=z*squeeze(G{n}(:,I(n),:));
       end
       z=z*G{N}(:,I(N));
       pos=num2cell(I);
       y=X(pos{:});
       T=z-y;
            % compute stochastic gradient
            for n=1:N
        if n==1
           g_nl=1;
           g_nr=squeeze(G{2}(:,I(2),:));
           for i=3:N
           g_nr=g_nr*squeeze(G{i}(:,I(i),:));
           end
           g1=T*g_nl'*g_nr';
        elseif n==2
           g_nl=reshape(G{1}(:,I(1),:),[1,size(G{1},3)]);
           g_nr=squeeze(G{3}(:,I(3),:));
           for i=4:N
               g_nr=g_nr*squeeze(G{i}(:,I(i),:));
           end
           g1=T*g_nl'*g_nr';
        elseif n==N
            g_nr=1;
            g_nl=reshape(G{1}(:,I(1),:),[1,size(G{1},3)]);
           for i=2:N-1
               g_nl=g_nl*squeeze(G{i}(:,I(i),:));
           end
             g1=T*g_nl'*g_nr';
        elseif n==N-1
            g_nr=squeeze(G{N}(:,I(N),:));
            g_nl=reshape(G{1}(:,I(1),:),[1,size(G{1},3)]);
            for i=2:N-2
                g_nl=g_nl*squeeze(G{i}(:,I(i),:));
            end
            g1=T*g_nl'*g_nr';
    else % n=2:N-2
                g_nl=reshape(G{1}(:,I(1),:),[1,size(G{1},3)]);
            for i=2:n-1
                g_nl=g_nl*squeeze(G{i}(:,I(i),:));
            end
                g_nr=squeeze(G{n+1}(:,I(n+1),:));
            for i=n+2:N
                g_nr=g_nr*squeeze(G{i}(:,I(i),:));
            end
            g1=T*g_nl'*g_nr';
        end
 
%                 if strcmp(sgdmethod, 'Hessian')
%                     g2 = g1./(b'.^2 + mu);
%                     g2 = reshape(g2, [r(n),1,r(n+1)]);
%                     grad{n}(:,I(n),:) = momentum * grad{n}(:,I(n),:)+ lrate * g2;
%                 end
%                 if strcmp(sgdmethod, 'Adadelta')
%                     g1 = reshape(g1,[r(n),1,r(n+1)]);
%                     v_t{n}(:,I(n),:) = beta2 * v_t{n}(:,I(n),:) + (1-beta2) * (g1.^2);
%                     grad{n}(:,I(n),:) =  (sqrt(m_t{n}(:,I(n),:)) + 1e-8) .* g1 ./ (sqrt(v_t{n}(:,I(n),:)) + 1e-8);
%                     m_t{n}(:,I(n),:) = beta2 * m_t{n}(:,I(n),:) + (1-beta2) * (grad{n}(:,I(n),:).^2);
%                 end
%                 if strcmp(sgdmethod, 'RMSprop')
%                     g1 = reshape(g1,[r(n),1,r(n+1)]);
%                     v_t{n}(:,I(n),:) =  beta2 * v_t{n}(:,I(n),:) + (1-beta2) * (g1.^2);
%                     grad{n}(:,I(n),:) = lrate * g1 ./ (sqrt(v_t{n}(:,I(n),:)) + 1e-8);
%                 end
                if strcmp(sgdmethod,'Adam')
                    if n==1
                    g1 = reshape(g1,[1,1,r(n)]);
                    elseif n==N
                    g1 = reshape(g1,[r(N-1),1]);
                    else
                    g1=reshape(g1,[r(n-1),1,r(n)]);
                    end
                    m_t{n}(:,I(n),:) =  beta1 * m_t{n}(:,I(n),:) + (1-beta1) * g1;
                    v_t{n}(:,I(n),:) =  beta2 * v_t{n}(:,I(n),:) + (1-beta2) * (g1.^2);
                   % m_t{j}(:,I(j),:) = m_t{j}(:,I(j),:)./(1-beta1.^(it+1));
                   % v_t{j}(:,I(j),:) = v_t{j}(:,I(j),:)./(1-beta2.^(it+1));
                    grad{n}(:,I(n),:) = lrate * m_t{n}(:,I(n),:) ./ (sqrt(v_t{n}(:,I(n),:)) + 1e-8);
                end
                
            end
            
% %             % gradient descent
            for n=1:N
                G{n}(:,I(n),:) = G{n}(:,I(n),:) - grad{n}(:,I(n),:);
            end
            
            if (mod(it+5000,1e5)==0 ) 
                LastG=G{1}; % for calculating err
            end
            
            if (mod(it,1e5)==0 ) 
                X_hat_w=double(W).*coreten2tt(G);
                err=abs(norm(G{1}(:)- LastG(:)))/norm(G{1}(:));
                RSE_all=RSE_fun(X,X_hat_w,double(W));
                RSE=RSE_all(2);
                 fprintf('it:%d, RSE=%.4f, tol=%.2e \n',it,RSE,err);
                
                if err<Tol
                        fprintf('Tol reached!\n');
                        break;
                end

            end
            
           
         
            
        end
            
        end
        G_opted=G;
    
end