function loglike=RCoDAlike(neibcell,kappa_serial,times,q)
%% this returns loglikelihood of Potts model
%% input
% neibcell: obtained in composedecom.m
% kappa_serial: spatial correlation in each split
% times: number of times of split
% q: number of component in Potts model.
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