function X=Xh2X(Xh,vs)
% vs{1} original image size
% vs{2} tensorized size
% vs{3} long array reshape
% vs{4} permute order
% vs{5} repermute order
X=reshape(permute(reshape(Xh,vs{3}),vs{5}),vs{1});
end