function deltaEs = energyPart(c, J)
    % Quick calculation of energy changes...
    ssums = circshift(c,1,1) + circshift(c,-1,1) +...
        circshift(c,1,2) + circshift(c,-1,2);
    % Clean up extraneous sums on the boundary of the grid
    ssums(1,:) = ssums(1,:) - c(end,:); 
    ssums(end, :) = ssums(end, :) - c(1,:);
    ssums(:,1) = ssums(:,1) - c(:,end); 
    ssums(:, end) = ssums(:, end) - c(:,1);

    deltaEs = 2 * J * (c.*ssums);
end