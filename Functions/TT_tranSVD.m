function [G_n,G_np1]=TT_tranSVD(big_G,G,n,e)
% G is cell mode core tensors
% n is from 1 to N-1, the mark of combination of n, n+1 core tensors
% big_G is the big core tensor with size of Rn-1In  x  In+1Rn+1
% e is the accuracy
% output G_n, G_np1 are 3-D tensors

a=size(G{n});
b=size(G{n+1});
N=numel(G);
 [u,s,v]=svd(big_G,'econ');
 s=diag(s);
 r=my_chop2(s,e);
 if r<=2
     r=2;
 elseif 2<r<=a(3)
     r=r;
 else
     r=a(3);
 end
 %r=8;
 u=u(:,1:r);
 s=s(1:r);
 v=v(:,1:r);
 G_n=reshape(u,[a(1),a(2),r]);
 if n~=N-1
 G_np1=reshape(diag(s)*v',[r,b(2),b(3)]);
 else
 G_np1=reshape(diag(s)*v',[r,b(2)]);
 end
end