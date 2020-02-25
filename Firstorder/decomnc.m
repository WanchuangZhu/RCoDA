% use simulated data from data.m 
clear
kappa=0.1:0.1:0.8;
for i=1:length(kappa)
    kappa_cell{i}=kappa(i);
end
m=256/8;n=m;
times=6;
Nhat=m/(2^(times/2));
q=2;
G_start=q;

iteration=6000;
parpool('local',8)
spmd
    for rep=1:200
        kappa_true=kappa_cell{labindex};
        fadd=fullfile('firstorderdata',strcat(num2str(m),num2str(q)),...
            strcat('Zst',num2str(kappa_true),num2str(rep),'.txt'));
        Z=dlmread(fadd);
        %         Z=randi(q,[m,n]);
        %         for i=1:5000
        %             Z=potts_prop(Z,kappa_true,q);
        %         end
        %% prior
        % we only infer \kappa in this script. The prior for \kappa is uniform
        % distribution between 0 and 1.
        %% starting points
        kappa_start=0.3;
        alpha_start=0.3;
        %% MCMC simulation
        file=strcat('NCwhole_',num2str(Nhat),'times',num2str(Nhat),num2str(q),'.txt');
        NCserial=importdata(file);
        NCserial=NCserial.data;
        % [breaks,func_cell]=fitncspline(NCserial);
        % [~,Z_res]=composelike(kappa_start,Z,times,G_start,alpha_start);
        [neibcell,Z_res]=composedecom(Z,times,q); % find the matrix blocks using RCoDA
        kappa_serial=kappa_start*alpha_start.^(0:(times-1));
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
            kappa_new=normrnd(kappa_start,0.03); % for 128*128, 0.06 is OK.
            if kappa_new>0 && kappa_new<0.9
                kappa_serial_start=kappa_start*alpha_start.^(0:(times-1));
                kappa_serial_new=kappa_new*alpha_start.^(0:(times-1));
                like_start=RCoDAlike(neibcell,kappa_serial_start,times,q);
                like_new=RCoDAlike(neibcell,kappa_serial_new,times,q);
                %             like_new=composelike(kappa_new,Z,times,G_start,alpha_start);
                %             like_start=composelike(kappa_start,Z,times,G_start,alpha_start);
                kappa_prob=like_new-like_start+pairs*(kappa_new-kappa_start)*alpha_start^(times)...
                    +ncintegnew(alpha_start^(times)*kappa_start,NCserial) -...
                    ncintegnew(alpha_start^(times)*kappa_new,NCserial);
                if kappa_prob>log(rand(1))
                    kappa_start=kappa_new;
                end
            end
            %% update alpha
            if labindex<4
                alpha_new=normrnd(alpha_start,0.09);
            else alpha_new=normrnd(alpha_start,0.04);
            end
            if alpha_new>0 && alpha_new<1
                kappa_serial_start=kappa_start*alpha_start.^(0:(times-1));
                kappa_serial_new=kappa_start*alpha_new.^(0:(times-1));
                like_alpha_start=RCoDAlike(neibcell,kappa_serial_start,times,q);
                like_alpha_new=RCoDAlike(neibcell,kappa_serial_new,times,q);
                alpha_prob=like_alpha_new-like_alpha_start + ...
                    pairs*kappa_start*(alpha_new^(times)-alpha_start^(times))+...
                    ncintegnew(alpha_start^(times)*kappa_start,NCserial) - ...
                    ncintegnew(alpha_new^(times)*kappa_start,NCserial);
                if alpha_prob>log(rand(1))
                    alpha_start=alpha_new;
                end
            end
            kappa_mat(iter)=kappa_start;
            alpha_mat(iter)=alpha_start;
        end
        kappa_est_var(rep)=var(kappa_mat(2000:end));
        alpha_est_var(rep)=var(alpha_mat(2000:end));
        kappa_est=mean(kappa_mat(2000:end));
        alpha_est=mean(alpha_mat(2000:end));
        kappa_rep(rep)= kappa_est;
        alpha_rep(rep)= alpha_est;
    end
end

N=[256,256]/2; q=2;
dsave(strcat('decomnc',num2str(N(1)),'_',num2str(q),'_',num2str(jobid),'.mat'));
exit
%% collect results
N=[256,256]; q=2;
a=[];b=a;aa=[];bb=[];
for jobid=1:20
    load(strcat('decomnc',num2str(N(1)),'_',num2str(q),'_',num2str(jobid),'.mat'));
    a=[a,cell2mat(kappa_rep')];
    b=[b,cell2mat(kappa_est_var')];
    aa=[aa,cell2mat(alpha_rep')];
    bb=[bb,cell2mat(alpha_est_var')];
end
for i=1:8
    kappa_rmse(:,i)=sqrt(mean((a(i,:)-0.1*i).^2));
end
latex(kappa_rmse,'nomath','%.3f')

latex(mean(a),'nomath','%.3f')
latex(sqrt(mean(b,2)),'nomath','%.3f')
latex(mean(aa),'nomath','%.3f')
latex(sqrt(mean(bb)),'nomath','%.3f')
%%

for i=1:8
    figure;
    plot(kappa_mat{i})
end
kappa_rep


for i=1:length(kappa_est)
    kappa_est_mat{i}=kappa_est{i};
end
kappa_re=cell2mat(kappa_est_mat); kappa_re
latex(kappa_re,'nomath','%.3f')


for i=1:length(alpha_est)
    alpha_est_mat{i}=alpha_est{i};
end
alphaaa=cell2mat(alpha_est_mat);
alphaaa
latex(alphaaa,'nomath','%.3f')


aa=cell2mat(kappa_est_mat);
kappa=0.05:0.05:0.8;
plot(aa,kappa)
bb=[aa;kappa];bb=bb(:,2:2:16);
dlmwrite('result_decomnc.txt',bb)
latex(bb,'nomath','%.3f')
%%%%%%%%%%%%%%%%%
for i=1:length(kappa_est_var)
    kappa_est_mat_var{i}=kappa_est_var{i};
end
kappa_re_var=cell2mat(kappa_est_mat_var); sqrt(kappa_re_var)
latex(sqrt(kappa_re_var),'nomath','%.3f')

for i=1:length(alpha_est_var)
    alpha_est_mat_var{i}=alpha_est_var{i};
end
alpha_re_var=cell2mat(alpha_est_mat_var); sqrt(alpha_re_var)
latex(sqrt(alpha_re_var),'nomath','%.3f')