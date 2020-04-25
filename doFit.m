function desiredY = doFit(desiredX, currentY)
    % doFit Interpolates currentY to fit desiredX.
    %   Input:
    %       currentY - Array of values.
    %       desiredX - Array with number of desired elements.
    %   Output:
    %       desiredY - interpolated result of currentY
    %           ode used.
    %   Assumptions:
    %       size(desiredX) >= currentY.
    
    % Calculate vector with same number of elements as currentY
    oldX = linspace(desiredX(1), desiredX(length(desiredX)), length(currentY));
    % Interpolate
    desiredY = interp1(oldX, currentY, desiredX);
end