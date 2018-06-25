function [X_high,vs]=VDT4(X)
% only for 2^n \times 2^n type image to 4 \times 4\times \cdots \times 4.
S=size(X);
n=log2(S(1));
vs=cell(1,4);

% vs{1} original image size
% vs{2} tensorized size
% vs{3} long array reshape
% vs{4} permute order
% vs{5} repermute order
vs{1}=size(X);                         % original image shape
vs{2}=[4*ones(1,n) S(3:end)];    % index_reshape
vs{3}=[2*ones(1,2*n) S(3:end)];    %index_arrangement
vs{4}=zeros(size(vs(2)));    % index_permute
vs{5}=vs(3);% index_repermute

v1=1:n; 
v2=(n+1):(2*n);
vs{4}=[v1;v2]; 
vs{4}=[vs{4}(:)' (2*n+1):numel(vs{3})];

s1=[];
s2=[];
for i=1:n
    s1=[s1 2*i-1];
    s2=[s2 2*i];
end
vs{5}=[s1 s2 (2*n+1):numel(vs{3})];

X_high=reshape(permute(reshape(X,vs{3}),vs{4}),vs{2});



end