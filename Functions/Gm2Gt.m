function Gt=Gm2Gt(Gm,r)
% from mode 2 matricization tensor core to original tensor core
% all are Nx1 cell mode
N=size(Gm,1);
Gt=cell(N,1);
for i=2:N-1
    Gt{i}=reshape(Gm{i},[size(Gm{i},1),r(i-1),r(i)]);
    Gt{i}=permute(Gt{i},[2,1,3]);
end
Gt{1}=reshape(Gm{1},[size(Gm{1},1),1,r(1)]);
Gt{1}=permute(Gt{1},[2,1,3]);
Gt{N}=reshape(Gm{N},[size(Gm{N},1),r(N-1),1]);
Gt{N}=permute(Gt{N},[2,1,3]);
end

    
