function[on,onCVI,cl]=TestCVI(varargin)
result=varargin{1};
kk=size(result,2)+1;
on=1;
onCVI=1;
for ii=2:ceil(sqrt(kk))
[avgDCVI,bins]=testWeight(ii,result);    
temp=min(onCVI,avgDCVI);
if onCVI~=temp
onCVI=temp;
on=ii;
cl=bins;
end
end
end