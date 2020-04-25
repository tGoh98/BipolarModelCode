function desiredY = doFit(desiredX, currentY)
    % ADD DOCUMENTATION
    oldX = linspace(desiredX(1), desiredX(length(desiredX)), length(currentY));
    desiredY = interp1(oldX, currentY, desiredX);    
end