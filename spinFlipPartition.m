function spinsToFlip = spinFlipPartition(spins, beta, partition)
    deltaEs= energyPart(spins, 1);
    expDeltaEs = 1./(1 + 1./exp(-beta * deltaEs)); 
    spinsToFlip = (expDeltaEs.*partition) > rand(size(spins));
end