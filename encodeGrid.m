function code = encodeGrid(spins)
    binaryRep = (spins + 1) / 2;
    binaryRep = binaryRep(:);
    code = 0;
    for i=1:length(binaryRep)
        code = code + binaryRep(i) * 2^(i-1);
    end
end