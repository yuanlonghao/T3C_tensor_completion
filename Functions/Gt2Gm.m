function Gm=Gt2Gm(Gt)
% Gt is a Nx1 cell of core tensors
N=size(Gt,1);
Gm=cell(N,1);
for i=1:N
    Gm{i}=tenmat(Gt{i},2);
    Gm{i}=double(Gm{i});
end
    
    