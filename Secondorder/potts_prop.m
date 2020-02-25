function new_potts=potts_prop(potts,beta,q,chess)
%this function returns a new potts sample based on the previous potts. this
%is to say it is a one-step update for potts sample. here we use gibbs
%algorithm to updata.
%% input
%potts : old Potts model
%beta spatial paramter in Potts model
% q: number of component in Potts model
% chess: label the Potts model according to coding method in our paper
N=size(potts);
sigma=potts;
randTol = 1;
for choice = 1:4
    choose=(chess==choice);
    spins_new = prop_new_potts(sigma,q,choose);
    sigma_new = (1-choose).*sigma + choose.*spins_new;
    %sigma_new is derived from sigma by changing half of the components to random numbers
    % Calculate number of neighbored nodes equal to the spin of all nodes
    movemat=[0,1;0,-1;1,0;-1,0;1,1;1,-1;-1,-1;-1,1];
    % need to compare with 8 neighbouhoods     % we don't consider the side effect.
    neighbors=0;
    for i =1:8
        d=(sigma==circshift(sigma, movemat(i,:)));
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
        neighbors=neighbors+d;
    end
    
    neighbors_new=0;
    for i =1:8
        d=(sigma_new==circshift(sigma_new, movemat(i,:)));
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
        neighbors_new=neighbors_new+d;
    end
    p = exp(-beta*(neighbors-neighbors_new));
    transitions = (rand(N) < p ).* (rand(N) < randTol) .*...
        choose .* (spins_new - sigma);
    sigma = sigma + transitions;
end
new_potts=sigma;
end