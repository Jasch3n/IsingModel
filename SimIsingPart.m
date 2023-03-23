function [configs, configGrids] = SimIsingPart(L, beta, J, N)
% [configs, configGrids] = SIMISINGPART(L, beta, J, N) uses MCMC to sample 
% from the distribution of a 2D Ising model in zero field on a LxL grid, 
% given beta and J. Compared to the SIMISING function, it is optimized by
% flipping multiple spins at once (partitioning the spin grid into
% non-intersecting sub-grids).
%
% The function returns [configs, configGrids] recorded from the
% simulation. Config is a list of grid encodings. ConfigGrids
% are the corresponding grids. It should be noted that for larger values of 
% L, the config numbers are too big for computers and NaN's will appear in
% the config array.
    c = 2*randi([0, 1],[L, L]) - 1; % initialize the grid randomly
    c(2,2) = -1;
    NSweeps = N; NSweepThreshold = 100; idx=1;
    
    configs = zeros(1, NSweeps - NSweepThreshold);
    % At each sweep, update the configuration of the whole grid using the 
    % Metropolis Hasting algorithm. The distribution of these
    % configurations eventually approach the distirbution used to construct
    % the MH algorithm.
    [X,Y] = meshgrid(1:L, 1:L);
    configGrids = zeros(L,L,NSweeps-NSweepThreshold);
    for itr = 1:NSweeps
        % Partition the grid into inter-crossing sub-grids
        red = (mod(X+Y, 2)==0) & ((X<L) & (Y<L)); red(L,L)=true; 
        blue = mod(X+Y, 2)==1 & ((X<L) & (Y<L));
        green = mod(X+Y, 2)==0 & ((X==L) | (Y==L)); green(L,L)=false;
        yellow = mod(X+Y, 2)==1 & ((X==L) | (Y==L));
        
        % Compute the spins to flip for each color
        redFlips = spinFlipPartition(c, beta, red);
        c(redFlips) = -1 * c(redFlips);

        blueFlips = spinFlipPartition(c, beta, blue);
        c(blueFlips) = -1 * c(blueFlips);

        greenFlips = spinFlipPartition(c, beta, green);
        c(greenFlips) = -1 * c(greenFlips);

        yellowFlips = spinFlipPartition(c, beta, yellow);
        c(yellowFlips) = -1 * c(yellowFlips);

        % Record the configuration sample.
        if itr > NSweepThreshold
            configs(idx) = encodeGrid(c);
            % This step maybe needs some justification...
            % But this is to prevent the model from getting stuck in
            % one particular ordered state when beta is large.
            % No matter the value of beta, these two configurations always
            % have the same probabilities!
            if beta>1 && (configs(idx) == 0 || configs(idx)==2^(L^2)-1)
                if rand() >= 0.5
                    c = -1 * c;
                end
            end
            configGrids(:,:,idx) = c;
            idx = idx + 1;
        end
    end
end