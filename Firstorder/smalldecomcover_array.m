clear
N=[32 32];
q=3;
kappa=0.1:0.1:0.8;
for i=1:length(kappa)
    kappa_cell{i}=kappa(i);
end
clear('kappa');
% jobid=getenv('PBS_ARRAYID');
parpool('local',8)



spmd
    % rng(str2num(jobid));
    kappa_true=kappa_cell{labindex};
    times=2;
    Nhat=N(1)/(2^(times/2));
%     alpha_start=0.3;
    file=strcat('NCwhole_',num2str(Nhat),...
        'times',num2str(Nhat),num2str(q),'.txt');
    NCserial=importdata(file);
    NCserial=NCserial.data;
    for reps=1:100
        Z=randi(q,N);
        for i=1:1000
            Z=potts_prop(Z,kappa_true,q);
        end
        G_start=q;
        %% prior
        % we only infer \kappa in this script. The prior for \kappa is uniform
        % distribution between 0 and 1.
        %% starting points
        kappa_start=0.5;
        iteration=6000;
        alpha_start=0.3;
        %% MCMC simulation
        kappa_new=0.2;
        [~,Z_res]=composelike(kappa_new,Z,times,G_start,alpha_start);
        pairs=0;
        d=(Z_res==circshift(Z_res, [0 1]));
        d(:,1)=0;
        pairs=pairs+sum(sum(d));
        d=(Z_res==circshift(Z_res, [0 -1]));
        d(:,size(Z_res,2))=0;
        pairs=pairs+sum(sum(d));
        d=(Z_res==circshift(Z_res, [1 0]));
        d(1,:)=0;
        pairs=pairs+sum(sum(d));
        d=(Z_res==circshift(Z_res, [-1 0]));
        d(size(Z_res,1),:)=0;
        pairs=pairs+sum(sum(d));
        pairs=pairs/2;
        for iter=1:iteration
            kappa_new=normrnd(kappa_start,0.05);
            if kappa_new>0 && kappa_new<0.9
                like_new=composelike(kappa_new,Z,times,G_start,alpha_start);
                like_start=composelike(kappa_start,Z,times,G_start,alpha_start);
                kappa_prob=like_new-like_start + pairs*alpha_start^(times)*...
                    (kappa_new-kappa_start)...
                    +ncintegnew(alpha_start^(times)*kappa_start,NCserial) - ...
                    ncintegnew(alpha_start^(times)*kappa_new,NCserial);
                if kappa_prob>log(rand(1))
                    kappa_start=kappa_new;
                end
            end
            %% update alpha
            alpha_new=normrnd(alpha_start,0.06);
            if alpha_new>0 && alpha_new<1
                like_new=composelike(kappa_start,Z,times,G_start,alpha_new);
                like_start=composelike(kappa_start,Z,times,G_start,alpha_start);
                alpha_prob=like_new-like_start + pairs* ...
                    kappa_start*(alpha_new^(times)-alpha_start^(times))+...
                    ncintegnew(alpha_start^(times)*kappa_start,NCserial) - ...
                    ncintegnew(alpha_new^(times)*kappa_start,NCserial);
                if alpha_prob>log(rand(1))
                    alpha_start=alpha_new;
                end
            end
            kappa_mat(iter)=kappa_start;
            alpha_mat(iter)=alpha_start;
        end
        kappa_mat_var=var(kappa_mat(2000:end));
        alpha_mat_var=var(alpha_mat(2000:end));
        kappa_rep(:,reps)=mean(kappa_mat(2000:end));
        alpha_rep(:,reps)=mean(alpha_mat(2000:end));
        kappa_var_rep(:,reps)=kappa_mat_var;
        alpha_var_rep(:,reps)=alpha_mat_var;
        kappa_interval(reps,:)=[quantile(kappa_mat(2000:end),0.025),quantile(kappa_mat(2000:end),0.975)];
    end
end
q=3;
N=[32,32];
dsave(strcat('smalldecomcover_',num2str(N(1)),num2str(q),'jobid=',num2str(jobid),'.mat'));

%%
cover=zeros(1,8);
for id=1:20
    load(strcat('smalldecomcover_323jobid=',num2str(id),'.mat'));
    for i=1:8
        kappa_int=kappa_interval{i};
        kappa_point=kappa_true{i};
        for j=1:size(kappa_int,1)
            if kappa_point>kappa_int(j,1) && kappa_point<kappa_int(j,2)
                cover(i)=cover(i)+1;
            end
        end
    end
end
latex(cover/200,'nomath','%.3f')

%% collect all the results
exit
for j=1:8
    kappa_tep(j,:)=kappa_rep{j};
end
kappa_result=kappa_tep;
clearvars kappa_rep

for j=1:8
    kappa_tep(j,:)=kappa_var_rep{j};
end
kappa_var_result=kappa_tep;

% clearvars -except kappa_rep kappa_hat kappa_var_hat kappa_result kappa_var
kappa_var=cell2mat(kappa_var);
kappa_est=cell2mat(kappa_est);
kappa_est=mean(kappa_result,2);

for i=1:8
    kappa=kappa_true(i);
    kappa=cell2mat(kappa);
    kappa_point=kappa_result(i,:);
    kappa_var_point=kappa_var_result(i,:);
    s=0;
    for j=1:length(kappa_point)
        varstd=sqrt(kappa_var_point(j));
        
        kappa_upper=kappa_point(j)+1.96*(varstd);
        kappa_lower=kappa_point(j)-1.96*(varstd);
        if kappa>kappa_lower && kappa<kappa_upper
            s=s+1;
        end
    end
    kappa_perc(i)=s;
end

% for i=1:8
%     kappa_int=kappa_result(i,:);
%     varstd=sqrt(kappa_var(i));
%     kappa_mean=kappa_est(i);
%     kappa_upper=kappa_mean+1.96*(varstd);
%     kappa_lower=kappa_mean-1.96*(varstd);
%     kappa_percent=(kappa_int<kappa_upper) .* (kappa_int>kappa_lower);
%     kappa_perc(i)=sum(kappa_percent);
% end
latex(kappa_perc,'nomath','%.0f')

for i=1:8
    kappa_int=kappa_result(i,:);
    % varstd=sqrt(kappa_var(i));
    kappa_mean=kappa_est(i);
    kappa_upper=kappa_interval{i}(2);
    kappa_lower=kappa_interval{i}(1);
    kappa_percent=(kappa_int<kappa_upper) .* (kappa_int>kappa_lower);
    kappa_perc(i)=sum(kappa_percent);
end
latex(kappa_perc,'nomath','%.0f')
