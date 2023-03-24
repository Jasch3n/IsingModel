function JPrime = applyScaling(J, Js, RJs)
    THRESHOLD = 1e-3;
    indices = find(Js<=J+THRESHOLD & Js>=J-THRESHOLD);
    JIndex = indices(1);
    JPrime = RJs(JIndex);
end