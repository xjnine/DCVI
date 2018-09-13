function[CL,on,onCVI]=DCVIClustering(varargin)
D=varargin{1};
[result,D1,LPS,CP,CL]=DCVI(D);
[on,onCVI,cl]=TestCVI(result);
CL(LPS)=cl;
for ii=1:size(D,1)
if(CL(ii)==0)
    CL(ii)=CL(CP(ii));
end
end
end