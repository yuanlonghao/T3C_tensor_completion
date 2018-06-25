function x=Gm2Gv(Gm)
N=size(Gm,1);
Gv=cell(N,1);
x=[];
for i=1:N
    Gv{i}=reshape(Gm{i},[numel(Gm{i}),1]);
    xp=Gv{i};
    x=[x;xp];
end

end