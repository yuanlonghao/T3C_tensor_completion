
% Given tensor X without missing entries, initial rank r, iteration k and
% error e, to do truncated SVD to descent tt-ranks to a optimal one
function [G_opted, r_opted]=TT_BALS(X,r,k,e)

S=size(X);
N=ndims(X);


%% Block ALS to optimize G and tt-ranks
G=initcoreten(S,r);
 for i=1:k
    for n=1:N-1
        [G_l,~]=G_nmode_devide(G,n); 
        [~,G_r1]=G_nmode_devide(G,n+1); 
        k_l=double(tenmat(G_r1,1));
        k_r=double(tenmat(G_l,n));
        G_m0=kron(k_l,k_r);     %G minus zero
        S1=S;% new dimension with big core
        S1(n)=S1(n)*S1(n+1);
        S1(n+1)=[];% S1 generated!!!
        X_n_np1=double(tenmat(tensor(reshape(X,S1)),n)); % X_(n,n+1) matrix
        big_G=X_n_np1*pinv(G_m0);  % ALS
        a=size(G{n});
        b=size(G{n+1});
        if n==N-1
        big_G_opted=reshape(big_G,[a(2)*b(2),a(1),1]);
        big_G_opted=permute(big_G_opted,[2,1,3]);
        big_G_opted=reshape(big_G_opted,[a(1)*a(2),b(2)]); % X_(n,n+1) matrix
        else
        big_G_opted=reshape(big_G,[a(2)*b(2),a(1),b(3)]);
        big_G_opted=permute(big_G_opted,[2,1,3]);
        big_G_opted=reshape(big_G_opted,[a(1)*a(2),b(2)*b(3)]); % X_(n,n+1) matrix
        end
        % truncated SVD
        [G_n,G_np1]=tranSVD(big_G_opted,G,n,e);
        G{n}=G_n;
        G{n+1}=G_np1;
    end
 end

 G_opted=G;
 r_opted=zeros(N-1,1);
for i=1:N-1
    r_opted(i)=size(G{i},3);
end
end