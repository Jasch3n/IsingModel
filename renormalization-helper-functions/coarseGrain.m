function scaledSpins = coarseGrain(spins)
% scaledSpins = COARSEGRAIN(spins) coarse grains an input square matrix
% with size multiple of 3 with 3x3 blocks using the majority rule.
    [Lx, Ly] = size(spins);
    assert(Lx == Ly, "Input to coarseGrain must be a square matrix.")
    assert(mod(Lx,3) == 0, "Input to coarseGrain must have size that is a multiple of 3!")
    L = Lx;
    scaledSpins = zeros(L/3, L/3); 
    for x=1:3:L
        for y=1:3:L
            block = spins(x:(x+2), y:(y+2));
            scaledSpins((x-1)/3 + 1,(y-1)/3 + 1) = sign(sum(sum(block)));
        end
    end
end