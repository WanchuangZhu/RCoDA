function [y,Z_res] = composedecom( Z,times,G_start )
%It returns the matrix blocks for calculation of likelihood using RCoDA
%% this returns matrix blocks for calculation of likelihood using RCoDA
%% input
% Z: data, Potts model
% times: number of split
% G_start: number of component in Potts model
%% output
% y: matrix blocks for calculation of likelihood using RCoDA
% Z_res: the final small lattice in RCoDA
N=size(Z);
mask=(Z>0);
% chess=repmat(eye(2),N/2);
for t=1:times
    chess=conschess(Z);
    n=0;v=[];
    %     if mod(t,2)==1
    for i=1:N(1)
        for j=1:N(2)
            if chess(i,j)==4 && mask(i,j)==1
                n=n+1;
                nei=showneibsecond(i,j,Z,G_start);
                
                v(n,:)=[Z(i,j),nei];
            end
        end
    end
    m=0;w=[];
    for i=1:N(1)
        for j=1:N(2)
            if chess(i,j)==3 && mask(i,j)==1
                m=m+1;
                nei=showneibfirst(i,j,Z,G_start);
                
                w(m,:)=[Z(i,j),nei];
            end
        end
    end
    Z=Z(chess==1|chess==2);
    Z=reshape(Z,N(1),N(2)/2);
    if mod(t,2)==1
        Z=Z';
    end
    if mod(size(Z,1),2)==1
        Z=[Z;zeros(1,size(Z,2))];
    end
    if mod(size(Z,2),2)==1
        Z=[Z,zeros(size(Z,1),1)];
    end
    mask=(Z>0);
    y{t}=[v;w];
    N=size(Z);
end
Z_res=Z;
end