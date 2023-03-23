function [configs, configGrids] = SimIsing(L, beta, J, N)
% [pmf, histconfigs, configs] = SIMISING(L, beta, J) uses MCMC to sample 
% from the distribution of a 2D Ising model in zero field on a LxL grid, 
% given beta and J. 
%
% The function returns [pmf, histconfigs, configs] recorded from the
% simulation. Histconfigs are the edge values used to approximate the pmf,
% while configs gives the true observations during the simulation.
    c = 2*randi([0, 1],[L, L]) - 1; % initialize the grid randomly 
    NSweeps = N; NSweepThreshold = 100; idx=1;
    
    configs = zeros(1, NSweeps - NSweepThreshold);
    % At each sweep, update the configuration of the whole grid using the 
    % Metropolis Hasting algorithm. The distribution of these
    % configurations eventually approach the distirbution used to construct
    % the MH algorithm.
    configGrids = zeros(L,L,NSweeps-NSweepThreshold);
    for itr = 1:NSweeps
        % flips one spin at a time, 
        % uniformly choose the transition probability 
        for i=1:(L^2)
            x=randi(L); y=randi(L); dExy = deltaE(c, [x,y], J);%             if 1/(1 + 1/exp(-beta*dExy)) > rand()
%                 c(x,y) = -1 * c(x,y);
%             end
            if exp(-beta * dExy) > rand()
                c(x,y) = -1 * c(x,y);
            end
        end

        if itr > NSweepThreshold
            configs(idx) = encodeGrid(c);
            configGrids(:,:,idx) = c;
            idx = idx + 1;
        end
    end
end