function[avgDCVI,bins,sortmst,weight2]=testWeight(varargin)
k=varargin{1};
mst=varargin{2};
avgDCVI=0;
DCVI=zeros(1,k);
edge=mst(3,:);
[Y,I]=sort(edge,2,'descend');
sortmst=mst(1:2,I);
sortmst(3,:)=Y;
cutnumber=k-1;
weight=zeros(1,k);
weight2=zeros(1,k);
cutmst=sortmst;
wmst=sortmst;

G=graph(cutmst(1,:),cutmst(2,:));
for ii=1:cutnumber
    a=cutmst(1,ii);
    b=cutmst(2,ii);
    G=rmedge(G,a,b);
end

bins=conncomp(G);
maxedge=sortmst(:,1:cutnumber);
cutpoint=zeros(3,cutnumber);
wmst(:,1:cutnumber)=[];
for ii=1:2
    for jj=1:cutnumber
      cutpoint(ii,jj)=bins(maxedge(ii,jj));
     
    end
end
cutpoint(3,:)=maxedge(3,:);

sep(1:k)=cutpoint(3,1);


for ii=1:k
   for jj=1:cutnumber
       if(cutpoint(1,jj)==ii||cutpoint(2,jj)==ii)
           sep(ii)=min(cutpoint(3,jj),sep(ii));
       end
   end
   
   disp(sep(ii));
end
disp('*******************************');

for ii=1:k
    temp2=0;
    temp=0;
    count=0;
 for jj=1:size(wmst,2)
  if bins(wmst(1,jj))==ii
      
      temp2=max(temp2,wmst(3,jj));
      temp=temp+wmst(3,jj);
      count=count+1;
%       disp(temp2);
  end
 end
  if((count-1)==0)
      weight(ii)=0;
      
  else
  weight(ii)=temp/(count-1);
  end
  weight2(ii)=temp2;
  
end
temp=0;
for ii=1:k
if (weight2(ii)==0)
   
   temp=min(weight2(weight2~=0));
   weight2(ii)=temp;
   
end

disp(weight2(ii));
 
end

% for ii=1:k
% DCVI(ii)=weight(ii)/sep(ii);
% end
% avgDCVI=sum(DCVI)/k;

% for ii=1:k
% DCVI(ii)=(sep(ii)-weight(ii))/(sep(ii)+weight(ii));
% end
% avgDCVI=sum(DCVI)/k;

% avgDCVI=(sum(weight)/k)/(sum(sep)/k);
% avgDCVI=max(weight)/min(sep);

for ii=1:k
DCVI(ii)=weight2(ii)/sep(ii);
end
avgDCVI=sum(DCVI)/k;

end