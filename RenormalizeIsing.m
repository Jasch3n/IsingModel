betas = [0, 0.3, 0.4, 0.5, 0.6, 50]; J=1; N=3000; L=81;
snapshots = zeros(L,L,length(betas));
for i=1:length(betas)
    tic;
    [configs, configGrids] = SimIsingPart(L, betas(i), J, N);
    snapshots(:,:,i) = configGrids(:,:,floor(2*N/3));
    toc;
end

singleCoarseSnapshots = zeros(L/3,L/3,length(betas));
for i=1:length(betas)
    unscaled = snapshots(:,:,i);
    for x=1:3:L
        for y=1:3:L
            block = unscaled(x:(x+2), y:(y+2));
            singleCoarseSnapshots((x-1)/3 + 1,(y-1)/3 + 1,i) = sign(sum(sum(block)));
        end
    end
end

doubleCoarseSnapshots = zeros(L/9, L/9, length(betas));
for i=1:length(betas)
    unscaled = singleCoarseSnapshots(:,:,i);
    for x=1:3:(L/3)
        for y=1:3:(L/3)
            block = unscaled(x:(x+2), y:(y+2));
            doubleCoarseSnapshots((x-1)/3 + 1,(y-1)/3 + 1,i) = sign(sum(sum(block)));
        end
    end
end

%% plotting the snapshots
figure();
for i=1:length(betas)
    subplot(2,3,i);
    imagesc(snapshots(:,:,i)); clim([-1 1]);
    title(sprintf("$\\beta=%.3f$", betas(i)), 'interpreter', 'latex');
end

figure();
for i=1:length(betas)
    subplot(2,3,i);
    imagesc(singleCoarseSnapshots(:,:,i)); clim([-1 1]);
    title(sprintf("$\\beta=%.3f$", betas(i)), 'interpreter', 'latex');
end

figure();
for i=1:length(betas)
    subplot(2,3,i);
    imagesc(doubleCoarseSnapshots(:,:,i)); clim([-1 1]);
    title(sprintf("$\\beta=%.3f$", betas(i)), 'interpreter', 'latex');
end

%% Getting J' = R(J)
% betas = [0, 0.3, 0.4, 0.5, 0.6, 4]; J=1; N=6000; L=81;
betas = 0.1:0.1:1; J=1; N=3000; L=81;
nativeUnscaledMsq = zeros(size(betas));
coarseMsq = zeros(size(betas));
for i=1:length(betas)
    display("Simulating...")
    tic;
    [configs, configGrids] = SimIsingPart(L, betas(i), J, N);
    toc;
    [~,~,numGrids] = size(configGrids);
    display("Coarse Graining...")
    tic;
    coarseGrids = zeros(L/3, L/3, length(betas));
    for j=1:numGrids
        coarseGrids(:,:,j) = coarseGrain(configGrids(:,:,j));
    end
    toc;
    display("Measuring Coarse Grained Spins...")
    tic;
    [~,unscaledMsq] = MeasureIsing(configGrids, 1);
    [~,Msq] = MeasureIsing(coarseGrids, 1);
    nativeUnscaledMsq(i) = mean(unscaledMsq);
    coarseMsq(i) = mean(Msq);
    toc;
end
coarseMsq

%% Native
betas = [0.1:0.1:1]; J=1; N=5000; L=81;
nativeMsq = zeros(1, length(betas));
for i=1:length(betas)
    tic;
    [configs, configGrids] = SimIsingPart(L/3, betas(i), 1, N); % ???
    [E,Msq] = MeasureIsing(configGrids, 1);
    nativeMsq(i) = mean(Msq);
    toc;
end
nativeMsq

%% Getting R(J) vs J
interpBetas = linspace(0.1,1,1000);
nativeMsqInterp=interp1(betas, nativeMsq, interpBetas, 'makima');
RJ = zeros(1,length(betas));
for i=1:length(betas)
    diffs = abs(nativeMsqInterp - coarseMsq(i));
    [~, indices] = sort(diffs, "ascend");
    RJ(i) = interpBetas(indices(1));
end
interpRJ = interp1(betas, RJ - 1e-6*(1:length(betas)).*RJ, interpBetas,'makima');
figure(); hold on;
plot(betas, RJ);
plot(interpBetas, interpRJ);
plot(interpBetas, interpBetas);

%% Testing
[configs, configGrids] = SimIsingPart(81, 0.6, J, 10000);
[E,Msq]=MeasureIsing(configGrids,1);
mean(Msq)