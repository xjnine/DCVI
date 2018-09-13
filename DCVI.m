function [result,D1,LPS,CP,CL]=DCVI(varargin)
D=varargin{1};
[CL,LPS,CP,CPV] = findCore(D);%寻找密度核
% showCPV(D,CPV,LPS,CL);
D1 = D(LPS,1:end);         %核心点
[G,result] = Kruskal(D1);
% MST(D1,result);%画出最小生成树

% cut=result;%切边
% %G=graph(cut(1,:),cut(2,:));
% % figure(3)
% % plot(G);
% [a,b]=max(cut(3,:));
% % dist(''+a);
% cut(:,b)=[];
% G1=graph(cut(1,:),cut(2,:));
% % figure(4)
% % plot(G1);
% bins = conncomp(G1);
end
function[CL,LPS,CP,CPV]=findCore(varargin)
D=varargin{1};

[Sup,NN,RNN,NNN,nb,A]=NaNSearching(D);
[CPV,CP,Pr1,Pr2,r1,rf,r2,Nei1,CL]=findCenter2(D,NN,NNN,A,nb);
[LPS,FLP,T2]=findDensityPeak(CP,D,r1,rf,Pr1,Pr2,nb,Sup,CL);
end%寻找密度核

function[Sup,NN,RNN,NNN,nb,A]=NaNSearching(varargin)
D=varargin{1};
r=1;
nb=zeros(size(D,1),1);
C=cell(size(D,1),1);
NN=cell(size(D,1),1);%初始化每个点的KNN邻居
RNN=cell(size(D,1),1);%初始化每个点的RKNN邻居
NNN=cell(size(D,1),1);%是NN和RNN的交集，也就每个点的
A=pdist2(D,D);
Numb1=0;
Numb2=0;
for ii=1:size(D,1)
   [sa,index]=sort(A(:,ii));
   C{ii}=[sa,index];
end
while(r<size(D,1))
 for kk=1:size(D,1)
     x=kk;
     y=C{x}(r+1,2);
     nb(y)= nb(y)+1;
     NN{x}=[NN{x},y];
     RNN{y}=[RNN{y},x];
 end
    Numb1=sum(nb==0);
    if Numb2~=Numb1
        Numb2=Numb1;
    else
       break; 
    end
    r=r+1;
end
for jj=1:size(D,1)
NNN{jj}=intersect(NN{jj},RNN{jj});
end
Sup=r;
end



function[CPV,CP,Pr1,Pr2,r1,rf,r2,Nei1,CL]=findCenter2(varargin)
D=varargin{1};
NN=varargin{2};
NNN=varargin{3};
A=varargin{4};
nb=varargin{5};


r1=[];
r2=[];
rf=[];
Pr1=[];
Pr2=[];
CPV=zeros(size(D,1),1);%用于噪声点的检查***

CP=[];


Nei1=cell(size(D,1),1);                                              
Nei2=cell(size(D,1),1);
 
for kk=1:size(D,1)
 CL(kk)=0;
 end 
    for ii=1:size(D,1)
       if ~isempty(NN{ii})
%         r2(ii)=mean(pdist2(D(ii,:),D(NN{ii},:)));%以点p与周围的自然邻的距离均值作为寻求convergent point的半径 ********************************
         r1(ii)=1*max(pdist2(D(ii,:),D(NN{ii},:)));
        r2(ii)=max(pdist2(D(ii,:),D(NN{ii},:)));%以点p与周围的自然邻的距离极大值作为寻求density core的半径
        rf=r1*0.95;
        Nei1{ii}=find(A(:,ii)<r1(ii)); 
        Nei2{ii}=find(A(:,ii)<rf(ii));         
        Pr1(ii)=size(Nei1{ii},1);
        Pr2(ii)=size(Nei2{ii},1);
       else
           r1(ii)=0;
           r2(ii)=0;
           rf=r1*0.95;
       end
    end
%      
   B=mean(r2)+2*std(r2); %标记噪声点，不让其因为动态半径过大，而被误选进入RCP,当数据集中没有噪声点时，不需要。** 
        for ii=1:size(D,1)
            if r2(ii)>B
                CL(ii)=-1;
            end
            if r2(ii)==0
                CL(ii)=-1;
            end
            if nb(ii)<2
                CL(ii)=-1;
            end
        end

    for jj=1:size(D,1)%将每个点的扫面半径内的点中的异常清除
     Nei1{ii}(find(CL(Nei1{ii})==-1))=[];
    end
    
       for ii=1:size(D,1)
           if ~isempty(Nei1{ii}) 
           if CL(ii)~=-1
            [~,y]=min(pdist2(D(Nei1{ii},:),mean(D(Nei1{ii},:),1)));   
             if (CPV(Nei1{ii}(y))==ii)           %对于没有邻居的点，以及互为收敛的两个点进行特殊处理
             CPV(ii)=ii;
              else
                CPV(ii)=Nei1{ii}(y);   
              end
             else
               CPV(ii)=ii;
           end
           else
               CPV(ii)=ii;
           end
        end
        
   for ii=1:size(D,1)                                                        %find the convergent point for point pi CP(i)       
     if CL(ii)~=-1
        CP(ii)=ii;
      while(CP(ii)~=CPV(CP(ii)))
         disp(CP(ii));
          CP(ii)=CPV(CP(ii));
      end
     else
         CP(ii)=ii;
     end
   end
   
end

function[LPS,FLP,T2]=findDensityPeak(varargin)
CP=varargin{1};
D=varargin{2};
r1=varargin{3};
rf=varargin{4};
Pr1=varargin{5};
Pr2=varargin{6};
nb=varargin{7};
Sup=varargin{8};
CL=varargin{9};
LPS=[];
FLP=[];
T2=[];
  for ii=1:size(D,1)
        if(CL(ii)~=-1)
        if(CP(ii)==ii)
            LPS=[LPS,ii];
% %             T2=size(D,2)*log(rf(ii)/r1(ii))+log(Pr1(ii)); 
%             if nb(ii)< Sup/2
%                   FLP=[FLP,ii];
%              end
        end
        end
       
 end
end

function[]=showCPV(varargin)
D=varargin{1};
CPV=varargin{2};
LPS=varargin{3};
CL=varargin{4};
figure(1)
plot(D(CL~=-1,1),D(CL~=-1,2),'g.','markersize',5);
hold on
for ii=1:size(D,1)
   if (CPV(ii)~=ii&&CL(ii)~=-1)
   x=[D(ii,1),D(CPV(ii),1)];
   y=[D(ii,2),D(CPV(ii),2)];
   plot(x,y,'-','color',[0.4 0.8 0.9],'markersize',20)
   end
end

plot(D(LPS,1),D(LPS,2),'r.','markersize',10);
end%画出数据集

function [T,result]=Kruskal(varargin)
%构造最小生成树
%算法：Kruskal
D=varargin{1};
%%做连接图（邻接矩阵）
a = zeros(size(D,1));
for i = 1:size(D,1)
    for j = 1:size(D,1)
        if i<j
            a(i,j)=sqrt((D(i,1)-D(j,1))*(D(i,1)-D(j,1))+(D(i,2)-D(j,2))*(D(i,2)-D(j,2)));
        end
    end
end

G = graph(a,'upper');  
[T,pred] = minspantree(G,'Method','sparse');
result = [ T.Edges.EndNodes(:,1), T.Edges.EndNodes(:,2),T.Edges.Weight]';
end%构造最小生成树

function[]=MST(varargin)
D=varargin{1};
result=varargin{2};
figure(2)
plot(D(:,1),D(:,2),'r.','markersize',10);
hold on


for i=1:size(result,2)
    x=[D(result(1,i),1),D(result(2,i),1)];
    y=[D(result(1,i),2),D(result(2,i),2)];
    plot(x,y,'b-');
    hold on
    
end
    

end%画出最小生成树

