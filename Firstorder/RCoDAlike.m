function loglike=RCoDAlike(neibcell,kappa_serial,times,q)
%% input
% neibcell: this is obtained in function composedecom.m
% kappa_serial: this is a vector. Each component in this vector corresponds
% to spatial correlation in each split
% times: number of splits
% q: number of components in Potts model.
%% output
 % loglike: loglikelihood of Potts model

loglike=0;
for i=1:times
    tem=neibcell{i};
    ttem=tem;
    ttem(:,1)=[];
    ttem=exp(ttem*kappa_serial(i));
    s=sum(ttem,2);
    S=repmat(s,[1,q]);
    scaletem=ttem./S;
    like=[];
    for j=1:size(scaletem,1)
        like(j)=scaletem(j,tem(j,1));
    end
    loglike=sum(log(like))+loglike;
end