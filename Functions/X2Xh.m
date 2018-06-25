function Xh=X2Xh(X,vs)
% vs{1} original image size
% vs{2} tensorized size
% vs{3} long array reshape
% vs{4} permute order
% vs{5} repermute order
Xh=reshape(permute(reshape(X,vs{3}),vs{4}),vs{2});
end