L = 5; % size of the grid 
h = 0; % magnetic field
J = 1; % constant coupling constant

%% Test energy computation
clc;
testGrid = 2*randi([0, 1],[L, L]) - 1;
testGridNew = testGrid; 
testGridNew(2,2) = -testGridNew(2,2);

EnergyDiff1 = energy(testGridNew, J) - energy(testGrid, J);
EnergyDiff2 = deltaE(testGrid, [2,2], J);
sprintf("energy and deltaE functions are self-consistent: %d", EnergyDiff1==EnergyDiff2)

%% Metropolis-Hasting MCMC
L=3;
c = 2*randi([0, 1],[L, L]) - 1; % initialize the grid randomly 
NSweeps = 10e5; NSweepThreshold = 1000;
beta = 0.1; % beta = 1/kT

configs = zeros(1, NSweeps - NSweepThreshold);
for itr = 1:NSweeps
    % uniformly choose the transition probability 
    for x=1:L 
        for y=1:L
            if exp(-beta * deltaE(c, [x,y], J)) > rand()
                c(x,y) = -c(x,y);
            end
        end
    end
    if itr > NSweepThreshold
        configs(itr) = encodeGrid(c);
    end
end

figure(); hold on;
histogram(configs, 'normalization', 'pdf')

%% Theoretical Config
probs = zeros(1, 2^(L*L));
energies = zeros(1, 2^(L*L));
cs = 0.:(2^(L*L)-1);
for c = cs
    energies(c+1) = energy(decodeGrid(c, L), J); 
    probs(c+1) = exp(-beta * energy(decodeGrid(c, L), J));
end
probs = probs / sum(probs);
plot(cs, probs);

%% Measuring Ising
beta = [0.01 0.02];
Es = zeros(1,length(beta));
EMsqs = zeros(1, length(beta));
L=27; J=1; N=2000;

for i=1:length(beta)
    [~, exMsq, exNrgy, ~, ~] = MeasureIsing(L, beta(i), J, N);
    Es(i) = exNrgy;
    EMsqs(i) = exMsq;
end

%% Plotting pdf's
figure(); subplot(2,1,1);
histogram(Msq, 20, 'normalization','probability');
title("$M^2$",'interpreter','latex');
subplot(2,1,2);
histogram(E, 50, 'normalization','probability');
title("$E$",'interpreter','latex');

%% Measuring Ising Model, Fast
L = 27; J=1; N=10000;
betas = [0.1:0.1:1, 4];
% betas = [0.1 0.2];
snapshots = zeros(L,L,length(betas));
Es = zeros(length(betas), N-100); % the 100 here is for the sweep cutoff
Msqs = zeros(length(betas), N-100);
for i=1:length(betas)
    tic;
    [configs, configGrids] = SimIsingPart(L, betas(i), J, 10000);
    [E, Msq] = MeasureIsing(configGrids, J);
    toc;
    sprintf("Finished Simulation and Measurements for beta=%.3f", betas(i))
    snapshots(:,:,i) = configGrids(:,:,floor(2*N/3));
    Es(i,:) = E;
    Msqs(i,:) = Msq;
end

%% Analysis of Measurements
figure();
for i=1:length(betas)
    subplot(3,4,i); 
    title(sprintf("$\beta=%.3f$", betas(i)), 'interpreter', 'latex');
    histogram(Es(i,:), 30, 'normalization', 'pdf'); xlim([-2 1]);
    xlabel("$\beta$",'interpreter','latex');
end

figure();
for i=1:length(betas)
    subplot(3,4,i); 
    title(sprintf("$\beta=%.3f$", betas(i)), 'interpreter', 'latex');
    histogram(Msqs(i,:), 30, 'normalization', 'pdf'); xlim([-0.1 1.1]);
    xlabel("$\beta$",'interpreter','latex');
end

ExEs = zeros(1,length(betas));
ExMsqs = zeros(1, length(betas));
for i=1:length(betas)
    ExEs(i) = mean(Es(i,:));
    ExMsqs(i) = mean(Msqs(i,:));
end

figure(); subplot(2,1,1);
plot(betas, ExEs); 
title("$\langle E\rangle$",'interpreter','latex');

subplot(2,1,2);
plot(betas, ExMsqs); 
xlabel("$\beta$",'interpreter','latex'); 
title("$\langle M^2\rangle$",'interpreter','latex');

% Can also compute the heat capacity...

