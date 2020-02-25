function vec=showneibfirst(x,y,Z,G_start)
% this function is used to count neighbours of a specific voxel. But it
% only applies to the voxel which is partial conditional on others.
% Therefore, in 8-neighbourhood structure, its only has 6 neighbourhood.
% if indicator==0
%     vec=zeros(1,G_start);
%     if x-1>0
%         vec(Z(x-1,y))=vec(Z(x-1,y))+1;
%     end
%     if x+1<=size(Z,1)
%         vec(Z(x+1,y))=vec(Z(x+1,y))+1;
%     end
%     if x-1>0 && y-1>0
%         vec(Z(x-1,y-1))=vec(Z(x-1,y-1))+1;
%     end
%     if x-1>0 && y+1<=size(Z,2)
%         vec(Z(x-1,y+1))=vec(Z(x-1,y+1))+1;
%     end
%     if x+1<=size(Z,1) && y+1<=size(Z,2)
%         vec(Z(x+1,y+1))=vec(Z(x+1,y+1))+1;
%     end
%     if x+1<=size(Z,1) && y-1>0
%         vec(Z(x+1,y-1))=vec(Z(x+1,y-1))+1;
%     end
% elseif indicator==1
mask=(Z>0);
vec=zeros(1,G_start);
if y-1>0 && mask(x,y-1)==1
    vec(Z(x,y-1))=vec(Z(x,y-1))+1;
end
if y+1<=size(Z,2) && mask(x,y+1)==1
    vec(Z(x,y+1))=vec(Z(x,y+1))+1;
end
if x-1>0 && y-1>0 && mask(x-1,y-1)==1
    vec(Z(x-1,y-1))=vec(Z(x-1,y-1))+1;
end
if x-1>0 && y+1<=size(Z,2) && mask(x-1,y+1)==1
    vec(Z(x-1,y+1))=vec(Z(x-1,y+1))+1;
end
if x+1<=size(Z,1) && y+1<=size(Z,2) && mask(x+1,y+1)==1
    vec(Z(x+1,y+1))=vec(Z(x+1,y+1))+1;
end
if x+1<=size(Z,1) && y-1>0 && mask(x+1,y-1)==1
    vec(Z(x+1,y-1))=vec(Z(x+1,y-1))+1;
end
end

        