function [E, Msq] = MeasureIsing(configGrids, J)
%     display("Measuring Ising")
%     % Sample from the distribution of the model
%     tic; [configs, configGrids] = SimIsing(L, beta, J, N); toc;
%     display("Measuring Energy...")
    % Measuring the Energy of the system <E>
    sizeGrids = size(configGrids); 
    L=sizeGrids(1); numSamples = sizeGrids(3);
%     tic; 
    E = zeros(1, numSamples);
    for i=1:length(E)
        E(i) = energyQ(configGrids(:,:,i), J) / (L*L);
    end 
%     toc;

%     display("Measuring Msq")
    % Measuring the magnetization <M^2>
%     tic; 
    Msq = zeros(1,numSamples);
    for i=1:length(Msq)
        M = sum(sum(configGrids(:,:,i))) / (L*L);
        Msq(i) = M*M;
    end
%     toc;
end