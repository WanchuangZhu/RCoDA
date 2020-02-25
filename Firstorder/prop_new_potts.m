function Z=prop_new_potts(Z_old,G_start)
%% this is used to propose a new Potts model
%% input
% Z_old: old Potts model
%G_start: number of component in Potts model
N=size(Z_old);
chess=Z_old~=0;
Z=Z_old+randi(G_start-1,N);
chess2=G_start*(Z==G_start);
Z=chess.*(mod(Z,G_start))+chess2;
end