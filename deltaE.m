function dE = deltaE(spins, spinToFlip, J)
% DELTAE computes the difference in energy given which spin is being
% flipped in a configuration.
    dE = 0;
    [m, n] = size(spins);
    x = spinToFlip(1); y = spinToFlip(2);
    for a=[-1 1]
        if (x+a >= 1) && (x+a <= m)
            dE = dE + 2*J*spins(x,y)*spins(x+a, y);
        end
    end

    for b = [-1 1]
        if (y+b >= 1) && (y+b <= n)
            dE = dE + 2*J*spins(x,y)*spins(x, y+b);
        end
    end
end