function E = energy(spins, J)
% ENERGY(spins, J) computes the energy of a grid of Ising spins, 
% which are either +1 or -1.
    [m,n] = size(spins);
    E=0;
    for x=1:m
        for y=1:n
            for a=[-1 1]
                if (x+a >= 1) && (x+a <= m)
                    E = E - J*spins(x,y)*spins(x+a, y);
                end
            end
            for b = [-1 1]
                if (y+b >= 1) && (y+b <= n)
                    E = E - J*spins(x,y)*spins(x, y+b);
                end
            end

        end
    end
    E = E/2;
end