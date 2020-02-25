function [y,Z_res] = composedecom_newversion( Z,times,G_start)
%composelike computes the likelihood by composing potts into several
%independent blocks. kappa is hyperparameter in potts model. Z is the potts
%model realization. Times is how many times we split the potts model.  
%% second conditional algorithm
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
                nei=showneibsecond(i,j,Z,G_start);
                
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