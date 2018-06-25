function [G_nl,G_nr]=G_nmode_devide(Gt,n)

N=size(Gt,1);
% n=1
if n==1
 G_nl=1;
 P=Gt{2};
 for i=2:N-1
     L=reshape(P,[numel(P)/size(Gt{i},3),size(Gt{i},3)]);
     R=reshape(Gt{i+1},[size(Gt{i+1},1),numel(Gt{i+1})/size(Gt{i+1},1)]);
     P=L*R;
 end
 index=[];
 for i=2:N
    index=[index size(Gt{i},2)];
 end
    index=[size(Gt{2},1) index];
    G_nr=reshape(P,index);
%n=N
elseif n==N
        G_nr=1;
        P=Gt{1};
    for i=1:N-2
        L=reshape(P,[numel(P)/size(Gt{i},3),size(Gt{i},3)]);
        R=reshape(Gt{i+1},[size(Gt{i+1},1),numel(Gt{i+1})/size(Gt{i+1},1)]);
        P=L*R;
    end
    index=[];
    for i=1:N-1
    index=[index size(Gt{i},2)];
    end
    index=[index size(Gt{N},1)];
    G_nl=reshape(P,index);
   else
       % other 
    P=Gt{1};
    for i=1:n-2
        L=reshape(P,[numel(P)/size(Gt{i},3),size(Gt{i},3)]);
        R=reshape(Gt{i+1},[size(Gt{i+1},1),numel(Gt{i+1})/size(Gt{i+1},1)]);
        P=L*R;
    end
    index=[];
    for i=1:n-1
    index=[index size(Gt{i},2)];
    end
    index=[index size(Gt{n-1},3)];
    G_nl=reshape(P,index);
    
    P=Gt{n+1};
    for i=n+1:N-1
        L=reshape(P,[numel(P)/size(Gt{i},3),size(Gt{i},3)]);
        R=reshape(Gt{i+1},[size(Gt{i+1},1),numel(Gt{i+1})/size(Gt{i+1},1)]);
        P=L*R;
    end
    index=[];
    for i=n+1:N
        index=[index size(Gt{i},2)];
    end
    index=[size(Gt{n+1},1) index];
    G_nr=reshape(P,index);
end
%  fprintf('size of G_nl is  ');
%  fprintf('%d  ' ,size(G_nl));
%  fprintf('\n');
%  fprintf('size of G_nr is  ');
%  fprintf('%d  ',size(G_nr));
%  fprintf('\n');
%  fprintf('\n');
end


