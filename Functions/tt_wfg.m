
function [f,Gradient] = tt_wfg(X,W,G,r)
% G is vector form core tensors
% Compute Y&Z
S=size(X);
Y = W.*X;
Gt=Gm2Gt(Gv2Gm(G,S,r),r);
Z = W.*coreten2tt(Gt);
Yt=tensor(Y);
Zt=tensor(Z);
% function value
f = 0.5 * norm(Yt)^2- innerprod(Zt,Yt) + 0.5 * norm(Zt)^2;
% gradient computation
N = ndims(X);
Gradient = cell(N,1);
T = Zt - Yt;
for n = 1:N
    [G_nl,G_nr]=G_nmode_devide(Gt,n); 
    G_m0=kron(double(tenmat(G_nr,1)),double(tenmat(G_nl,n)));     %G minus zero
    Gradient{n} = double(tenmat(T,n))*G_m0';
end          

Gradient=Gm2Gv(Gradient);
end
 
