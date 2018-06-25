function tt=coreten2tt(G)
% combine tt-cores back to tt-tensor
% input must be cell mode core tensors with size Nx1
P=G{1}; %initialize P
for n=1:size(G,1)-1
    L=reshape(P,[numel(P)/size(G{n},3),size(G{n},3)]);
    R=reshape(G{n+1},[size(G{n+1},1),numel(G{n+1})/size(G{n+1},1)]);
    P=L*R;
end
%reshape P to original tensor size
S=[];
for j=1:size(G,1)
    S=[S size(G{j},2)];
end
 tt=reshape(P,S);
end