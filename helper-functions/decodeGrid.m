function grid = decodeGrid(code, L)
    grid = zeros(L, L);
    binRep = (reshape(reverse(dec2bin(code, L*L)), [L, L]) == '1');
    grid = 2*binRep - 1;
end