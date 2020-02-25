clear
kappa=0.1:0.1:0.3;
for i=1:length(kappa)
    kappa_cell{i}=kappa(i);
end
m=256/8;n=m;
times=6;
q=2;
G_start=q;
chess=zeros(m,n);
chess=conschess(chess);
clear('kappa');
% jobid=getenv('PBS_ARRAYID');
% jobid=str2double(jobid);

parpool('local',3)
spmd
    chromatic=4;
    kappa_true=kappa_cell{labindex};
    for rep=1:200
        fadd=fullfile('secondorderdata',strcat(num2str(m),num2str(q)),...
            strcat('Znd',num2str(kappa_true),num2str(rep),'.txt'));
        Z=dlmread(fadd);
%         Z=randi(q,[m,n]);
%         for i=1:5000
%             Z=potts_prop(Z,kappa_true,q,chess);
%         end
        %% starting points
        kappa_start=0.2;
        iteration=6000;
        Nhat=m/(2^(times/2));
        alpha_start=0.3;
        %% MCMC simulation
        file=strcat('NCwhole_',num2str(Nhat), 'times',num2str(Nhat),num2str(q),'.txt');
        NCserial=importdata(file);
        NCserial=NCserial.data;
        % [~,Z_res]=composelike_newversion(kappa_start,Z,times,G_start,alpha_start);
        [neibcell,Z_res]=composedecom(Z,times,q); % find the matrix blocks using RCoDA
        
        pairs=0;
        movemat=[0,1;0,-1;1,0;-1,0;1,1;1,-1;-1,-1;-1,1];
        % need to compare with 8 neighbouhoods
        for i =1:8
            d=(Z_res==circshift(Z_res, movemat(i,:)));
            if movemat(i,1)==1
                d(1,:)=0;
            elseif movemat(i,1)==-1
                d(end,:)=0;
            end
            if movemat(i,2)==1
                d(:,1)=0;
            elseif movemat(i,2)==-1
                d(:,end)=0;
            end
            pairs=pairs+sum(sum(d));
        end
        pairs=pairs/2;
        for iter=1:iteration
            kappa_new=normrnd(kappa_start,0.01);
            if kappa_new>0 && kappa_new<0.5
                kappa_serial_start=kappa_start*alpha_start.^(0:(times-1));
                kappa_serial_new=kappa_new*alpha_start.^(0:(times-1));
                like_start=RCoDAlike(neibcell,kappa_serial_start,times,q);
                like_new=RCoDAlike(neibcell,kappa_serial_new,times,q);
                %             like_new=composelike_newversion(kappa_new,Z,times,G_start,alpha_start);
                %             like_start=composelike_newversion(kappa_start,Z,times,G_start,alpha_start);
                kappa_prob=like_new-like_start+pairs*(kappa_new-kappa_start)*alpha_start^times...
                    +ncintegnew(kappa_start*alpha_start^times,NCserial) - ...
                    ncintegnew(kappa_new*alpha_start^times,NCserial);
                if kappa_prob>log(rand(1))
                    kappa_start=kappa_new;
                end
            end
            %% update alpha
            alpha_new=normrnd(alpha_start,0.01);
            if alpha_new>0 && alpha_new<1
                kappa_serial_start=kappa_start*alpha_start.^(0:(times-1));
                kappa_serial_new=kappa_start*alpha_new.^(0:(times-1));
                like_alpha_start=RCoDAlike(neibcell,kappa_serial_start,times,q);
                like_alpha_new=RCoDAlike(neibcell,kappa_serial_new,times,q);
                alpha_prob=like_alpha_new-like_alpha_start + ...
                    pairs*(kappa_start*alpha_new^times-kappa_start*alpha_start^times)...
                    +ncintegnew(kappa_start*alpha_start^times,NCserial)- ...
                    ncintegnew(kappa_start*alpha_new^times,NCserial);
                if alpha_prob>log(rand(1))
                    alpha_start=alpha_new;
                end
            end
            kappa_mat(iter)=kappa_start;
            alpha_mat(iter)=alpha_start;
        end
        kappa_est_var=var(kappa_mat(2000:end));
        alpha_est_var=var(alpha_mat(2000:end));
        kappa_est=mean(kappa_mat(2000:end));
        alpha_est=mean(alpha_mat(2000:end));
        kappa_rep(rep)=kappa_est;
        alpha_rep(rep)=alpha_est;
    end
end

N=[256,256]/8;q=2;
dsave(strcat('decomncsecond',num2str(N(1)),num2str(q),num2str(jobid),'_M','.mat'));

exit
%%
N=[256 256]/2;q=3;
a=[];b=a;aa=[];bb=[];
for jobid=1:20
    load(strcat('decomncsecond',num2str(N(1)),num2str(q),num2str(jobid),'_M','.mat'));
    a=[a,cell2mat(kappa_rep')];
    b=[b,cell2mat(kappa_est_var')];
    aa=[aa,cell2mat(alpha_rep')];
    bb=[bb,cell2mat(alpha_est_var')];
end
for i=1:3
    kappa_rmse(:,i)=sqrt(mean((a(i,:)-0.1*i).^2));
end
latex(kappa_rmse,'nomath','%.3f')

latex(mean(a),'nomath','%.3f')
latex(sqrt(mean(b,2)),'nomath','%.3f')
latex(mean(aa),'nomath','%.3f')
latex(sqrt(mean(bb)),'nomath','%.3f')

%%
for i=1:length(kappa_est)
    kappa_est_mat{i}=kappa_est{i};
end
kappa_re=cell2mat(kappa_est_mat); 
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