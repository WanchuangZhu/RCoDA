function mat=conschess(z) 
% this is used to construct the chess matrix for any field, where z is
% random field. The size of z should be power of 2, therefore it can be
% divided by 2 always
%% input
% z: Potts model
%% output
% mat: the same size with z, it divides the Potts model into several parts
% according to coding method.
N=size(z);
% chess=zeros(N);
rowarray1=1:2:(N(1)-1);
rowarray2=2:2:N(1);
colarray1=1:2:(N(2)-1);
colarray2=2:2:N(2);
% switch n
%     case 0
%         for i=rowarray1
%             for j=colarray1
%                 mat(i,j)=1;
%             end
%         end
%         
%         for i=rowarray1
%             for j=colarray2
%                 mat(i,j)=2;
%             end
%         end
%         for i=rowarray2
%             for j=colarray1
%                 mat(i,j)=3;
%             end
%         end
%         for i=rowarray2
%             for j=colarray2
%                 mat(i,j)=4;
%             end
%         end
%     case 1
        for i=rowarray1
            for j=colarray1
                mat(i,j)=1;
            end
        end
        for i=rowarray1
            for j=colarray2
                mat(i,j)=3;
            end
        end
        for i=rowarray2
            for j=colarray1
                mat(i,j)=2;
            end
        end
        for i=rowarray2
            for j=colarray2
                mat(i,j)=4;
            end
        end
end
