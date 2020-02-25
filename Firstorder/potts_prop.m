function new_potts=potts_prop(potts,beta,q)
%this function returns a new potts sample based on the previous potts. this
%is to say it is a one-step update for potts sample. here we use gibbs
%algorithm to updata.
%% input
% potts: old Potts model, wait to update
% beta: spatial parameter in Potts model
% q: number of component in Potts model
N=size(potts);
switch length(N)
    case 3
        % q=max(max(max(potts)));
        sigma=potts;
        randTol = 1;
        for choice = 0:1
            threedeye(:,:,1)=eye(2);
            threedeye(:,:,2)=circshift(eye(2),[0,1]);
            chess = (circshift(repmat(threedeye,N/2),[0 0 choice]));
            spins_new = prop_new_potts(sigma,q);
            sigma_new = (1-chess).*sigma + chess.*spins_new; 
            %sigma_new is derived from sigma by changing half of the components to random numbers
            
            % Calculate number of neighbored nodes equal to the spin of all nodes
            neighbors = (sigma==circshift(sigma, [ 0 1 0])) + ... % compare with their left neighbours
                (sigma==circshift(sigma, [ 0 -1 0])) + ...  %compare with their right neighbours
                (sigma==circshift(sigma, [ 1 0 0])) + ...  %compare with their backward neighbours
                (sigma==circshift(sigma, [ -1 0 0]))+... % compare with their foreward neighbours
                (sigma==circshift(sigma, [ 0 0 1]))+... %compare with their downward neighbours
                (sigma==circshift(sigma, [ 0 0 -1]));  %compare with their upward neighbours.the neighbors here is a matrix
            neighbors_new = (sigma_new==circshift(sigma_new, [ 0 1 0])) + ... % compare with their left neighbours
                (sigma_new==circshift(sigma_new, [ 0 -1 0])) + ...  %compare with their right neighbours
                (sigma_new==circshift(sigma_new, [ 1 0 0])) + ...  %compare with their backward neighbours
                (sigma_new==circshift(sigma_new, [ -1 0 0]))+... % compare with their foreward neighbours
                (sigma_new==circshift(sigma_new, [ 0 0 1]))+... %compare with their downward neighbours
                (sigma_new==circshift(sigma_new, [ 0 0 -1]));  %compare with their upward neighbours.the neighbors here is a matrix 
            % Calculate the change in energy of flipping a spin 
            %     DeltaE = 2 .* neighbors + 4 + 2 .* B .* sigma; 
            %     DeltaE = 2*sigma .* neighbors + 2 .* B .* sigma;
            %    DeltaE = 2 * sigma .* neighbors; %(critical value for beta = 0.4406868)
            % DeltaE = sigma .* neighbors;  %(critical value for beta = 0.8813736)
            % Calculate the transition probabilities 
            p = exp(-beta*(neighbors-neighbors_new));
            % Decide which transitions will occur 
            %transitions = (U < p ).*(rand(N) < randTol) * -2 + 1; 
            transitions = (rand(N) < p ).* (rand(N) < randTol) .*... 
                chess .* (spins_new - sigma);
            % Perform the transitions  
            sigma = sigma + transitions;
        end
        new_potts=sigma;
    case 2
        % q=max(max(potts));
        sigma=potts;
        randTol = 1;
        for choice = 0:1
            threedeye=eye(2);
            % threedeye(:,:,2)=circshift(eye(2),[0,1]);
            chess = (circshift(repmat(threedeye,N/2),[ 0 choice]));
            spins_new = prop_new_potts(sigma,q);
            sigma_new = (1-chess).*sigma + chess.*spins_new; 
            %sigma_new is derived from sigma by changing half of the components to random numbers
            
            % Calculate number of neighbored nodes equal to the spin of all nodes
            neighbors=0;
        d=(sigma==circshift(sigma, [0 1]));
        d(:,1)=0;
        neighbors=neighbors+d;
        d=(sigma==circshift(sigma, [0 -1]));
        d(:,size(sigma,2))=0;
        neighbors=neighbors+d;
        d=(sigma==circshift(sigma, [1 0]));
        d(1,:)=0;
        neighbors=neighbors+d;
        d=(sigma==circshift(sigma, [-1 0]));
        d(size(sigma,1),:)=0;
        neighbors=neighbors+d;
        % neighbors=neighbors/2;  
        
            neighbors_new=0;
        d=(sigma_new==circshift(sigma_new, [0 1]));
        d(:,1)=0;
        neighbors_new=neighbors_new+d;
        d=(sigma_new==circshift(sigma_new, [0 -1]));
        d(:,size(sigma_new,2))=0;
        neighbors_new=neighbors_new+d;
        d=(sigma_new==circshift(sigma_new, [1 0]));
        d(1,:)=0;
        neighbors_new=neighbors_new+d;
        d=(sigma_new==circshift(sigma_new, [-1 0]));
        d(size(sigma_new,1),:)=0;
        neighbors_new=neighbors_new+d;
        % neighbors_new=neighbors_new/2;
            %     DeltaE = 2 .* neighbors + 4 + 2 .* B .* sigma; 
            %     DeltaE = 2*sigma .* neighbors + 2 .* B .* sigma;
            %    DeltaE = 2 * sigma .* neighbors; %(critical value for beta = 0.4406868)
            % DeltaE = sigma .* neighbors;  %(critical value for beta = 0.8813736)
            % Calculate the transition probabilities 
            p = exp(-beta*(neighbors-neighbors_new));
            % Decide which transitions will occur 
            %transitions = (U < p ).*(rand(N) < randTol) * -2 + 1; 
            transitions = (rand(N) < p ).* (rand(N) < randTol) .*... 
                chess .* (spins_new - sigma);
            % Perform the transitions  
            sigma = sigma + transitions;
        end
        new_potts=sigma;
end