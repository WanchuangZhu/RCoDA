clear
latticesize=[32 128 256];

%parpool('local',3)
rng(jobid) %%%%%%%%%%%% set different random seed for each generation.

kappa=0.1:0.1:0.3;
for jobid=1:200
    for sz=1:length(latticesize)
        m=latticesize(sz);n=m;
        chess=zeros(m,n);
        chess=conschess(chess);
        q=3;
        G_start=q;
        for i=1:length(kappa)
            kappa_true=kappa(i);
            Z=randi(q,[m,n]);
            for it=1:1000
                Z=potts_prop(Z,kappa_true,q,chess);
            end
            fadd=fullfile('secondorderdata',strcat(num2str(m),num2str(q)),...
                strcat('Znd',num2str(kappa_true),num2str(jobid),'.txt'));
            dlmwrite(fadd,Z);
        end
    end
end