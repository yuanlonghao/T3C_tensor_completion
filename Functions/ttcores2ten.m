function T=ttcores2ten(tt_core,X,rank)
% from cell mode tt cores and original data X and (N+1)x1 mode core ranks to new data tensor T
P=tt_core{1}; %initialize P
for n=1:size(tt_core,1)-1
    L=reshape(P,[numel(P)/rank(n+1),rank(n+1)]);
    R=reshape(tt_core{n+1},[rank(n+1),numel(tt_core{n+1})/rank(n+1)]);
    P=L*R;
end
%reshape P to original tensor size
T=reshape(P,size(X));
end