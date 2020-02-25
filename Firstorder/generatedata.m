clear
latticesize=[32 128 256];

kappa=0.1:0.1:0.8;
for jobid=1:200 % 200 generations
    for sz=1:length(latticesize)
        m=latticesize(sz);n=m;
        q=2;
        G_start=q;
        for i=1:8
            kappa_true=kappa(i);
            Z=randi(q,[m,n]);
            for it=1:1000
                Z=potts_prop(Z,kappa_true,q);
            end
            fadd=fullfile('firstorderdata',strcat(num2str(m),num2str(q)),...
                strcat('Zst',num2str(kappa_true),num2str(jobid),'.txt'));
            % Give name for each simulated data, and put in different
            % folder
            dlmwrite(fadd,Z);
        end
    end
end
