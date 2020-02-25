clear
kappa=0.1:0.1:0.3;
for i=1:length(kappa)
    kappa_cell{i}=kappa(i);
end
chromatic=4;
N=[32 32];

q=3;
G_start=q;

clear('kappa');
parpool('local',3)

spmd
   
    kappa_true=kappa_cell{labindex};
    m=32;n=m;
    chess=zeros(m,n);
    chess=conschess(chess);

    %% MCMC simulation

    %% coverage
    for reps=1:200
        Z=randi(q,N);
        for i=1:5000
            Z=potts_prop(Z,kappa_true,q,chess);
        end
        kappa_start=0.2;
        iteration=6000;times=6;
        alpha_start=0.3;
        %%
        alpha_start=0.3;
        Nhat=N(1)/(2^(times/2));
        file=strcat('NCwhole_',num2str(Nhat),'times',num2str(Nhat),num2str(q),'.txt');
        NCserial=importdata(file);
        NCserial=NCserial.data;
        [~,Z_res]=composelike_newversion(kappa_start,Z,times,G_start,alpha_start);
        [neibcell,~]=composedecom(Z,times,q); %
        pairs=0;
        movemat=[0,1;0,-1;1,0;-1,0;1,1;1,-1;-1,-1;-1,1];
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
            kappa_new=normrnd(kappa_start,0.06);% good, maybe could be larger since accept ratio is still high
            if kappa_new>0 && kappa_new<0.4
                kappa_serial_start=kappa_start*alpha_start.^(0:(times-1));
                kappa_serial_new=kappa_new*alpha_start.^(0:(times-1));
                like_start=RCoDAlike(neibcell,kappa_serial_start,times,q);
                like_new=RCoDAlike(neibcell,kappa_serial_new,times,q);
%                 like_new=composelike(kappa_new,Z,times,G_start,alpha_start);
%                 like_start=composelike(kappa_start,Z,times,G_start,alpha_start);

                kappa_prob=like_new-like_start+pairs*(kappa_new-kappa_start)*alpha_start^times...
                        +ncintegnew(kappa_start*alpha_start^times,NCserial) - ...
                        ncintegnew(kappa_new*alpha_start^times,NCserial);
                if kappa_prob>log(rand(1))
                    kappa_start=kappa_new;
                end
            end
            %% update alpha
            alpha_new=normrnd(alpha_start,0.15);
            if alpha_new>0 && alpha_new<1
                kappa_serial_start=kappa_start*alpha_start.^(0:(times-1));
                kappa_serial_new=kappa_start*alpha_new.^(0:(times-1));
                like_alpha_start=RCoDAlike(neibcell,kappa_serial_start,times,q);
                like_alpha_new=RCoDAlike(neibcell,kappa_serial_new,times,q);
%                 like_alpha_new=composelike(kappa_start,Z,times,G_start,alpha_new);
%                 like_alpha_start=composelike(kappa_start,Z,times,G_start,alpha_start);
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
        kappa_est_rep(:,reps)=mean(kappa_mat(2000:end));
        kappa_var_rep(:,reps)=var(kappa_mat(2000:end));
        tem=kappa_mat(2000:end);
        kappa_interval(reps,:)=[quantile(tem,0.025),quantile(tem,0.975)];
    end
end

N=[32,32];q=3;
dsave(strcat('smalldecomncsecondcover',num2str(N(1)),num2str(q),'_M','.mat'));

exit

cover=zeros(1,3);
for id=1:4
    load(strcat('smalldecomncsecondcover323',num2str(id),'_C','.mat'));
    for i=1:3
        kappa_int=kappa_interval{i};
        kappa_point=kappa_true{i};
        for j=1:size(kappa_int,1)
            if kappa_point>kappa_int(j,1) && kappa_point<kappa_int(j,2)
                cover(i)=cover(i)+1;
            end
        end
    end
end
cover/200
latex(cover/200,'nomath','%.3f')
%%
for i=1:3
    kappa_int=kappa_est_rep{i};
    kappa_upper=kappa_est{i}+1.96*sqrt(kappa_var{i});
    kappa_low=kappa_est{i}-1.96*sqrt(kappa_var{i});
    kappa_perc(i)=sum((kappa_int<kappa_upper ).* (kappa_int>kappa_low));
end
kappa_perc

for i=1:3
    kappa_int=kappa_est_rep{i};
    kappa_upper=kappa_interval{i}(2);
    kappa_low=kappa_interval{i}(1);
    kappa_perc(i)=sum((kappa_int<kappa_upper ).* (kappa_int>kappa_low));
end
kappa_perc
latex(kappa_perc,'nomath','%.0f')

