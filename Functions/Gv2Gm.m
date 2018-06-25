function Gm=Gv2Gm(Gv,S,r)
N=numel(S);
Gm=cell(N,1);
I=cell(N,1); % I is the number of elements of every G cores
for i=2:N-1 %  2 to N-1
    I{i}=S(i)*r(i-1)*r(i);
end
I{1}=S(1)*r(1);
I{N}=S(N)*r(N-1);
%
index=[];
add=I{1};
for j=1:N-1
    index=[index add];
    add=add+I{j+1};
end
index=[0 index add];% index is a N+1 element index
%
for k=1:N
    Gm{k}=Gv(index(k)+1:index(k+1),1);
end
%
for d=2:N-1
Gm{d}=reshape(Gm{d},[S(d),r(d-1)*r(d)]);
end
Gm{1}=reshape(Gm{1},[S(1),r(1)]);
Gm{N}=reshape(Gm{N},[S(N),r(N-1)]);
end
                 






% Gm{1}=Gv(1:r*r*size(X,1),1);
% Gm{1}=reshape(Gm{1},[size(X,1),r,r]);
% Gm{1}=permute(Gm{1},[2,1,3]);
% for i=1:N-1
%     num1=sum(Sz(1:i));
%     num2=sum(Sz(1:(i+1)));
%     Gm{i+1}=Gv(r^2*num1+1:r^2*num2,1);
%     Gm{i+1}=reshape(Gm{i+1},[size(X,i+1),r,r]);
%     Gm{i+1}=permute(Gm{i+1},[2,1,3]);
% end