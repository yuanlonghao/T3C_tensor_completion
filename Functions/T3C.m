
function [X_hat,G_opted]=T3C(varargin)

% The 1st and 2nd input must be X and W
if (nargin == 0||nargin == 1)
    X_hat=0;
    G_opted=0;
    fprintf('Not enough input.\n');
else
    
    ip = inputParser;
    ip.addParamValue('Tol', 0, @isscalar);
    ip.addParamValue('Rank', 16*ones(1,ndims(varargin{1})-1), @ismatrix);
    ip.addParamValue('MaxIter', 5e2, @isscalar);
    ip.addParamValue('Alg', 'TTWOPT', @ischar);
    ip.parse(varargin{3:end});
    
    Tol = ip.Results.Tol;
    Rank = ip.Results.Rank;
    MaxIter = ip.Results.MaxIter;
    Alg = ip.Results.Alg;
    X=varargin{1};
    W=varargin{2};
    
    S=size(X);
    N=numel(S);
    r=Rank;
%     r_init=Rank;
%     [G0,r]=TT_G_init(X,W,r_init); % initial r and G by truncate-SVD
    G0=Abs_initcoreten(S,r);
    
            % TTWOPT
    if strcmp(Alg,'TTWOPT')
        fprintf('------ TTWOPT ------ \n\n');
        G0_v=1/N/N*Gm2Gv(Gt2Gm(G0));
        options = lbfgs('defaults');
        options.MaxIters = MaxIter;
        options.StopTol = Tol;
        options.RelFuncTol = 1.0e-100;
        options.DisplayIters =100;
        options.MaxFuncEvals=1.0e+100;
        out=lbfgs(@(G) tt_wfg(W.*X,W,G,r),G0_v, options);
        
        G_opted=Gm2Gt(Gv2Gm(out.X,S,r),r); % core tensors
        X_hat=coreten2tt(G_opted);
        X_hat(W==1)=X(W==1);
        
           % TTSGD     
    elseif strcmp(Alg,'TTSGD')
        fprintf('------ TTSGD ------ \n\n');
        G_opted=tt_sgd(X,W,'Rank',r,'MaxIter',MaxIter,'Tol',Tol);
        X_hat=coreten2tt(G_opted);
        X_hat(W==1)=X(W==1);
    end
end
end