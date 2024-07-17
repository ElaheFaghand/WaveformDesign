function [roundedNumber] = number (number)
    % Number to be rounded
    % number = 0.25477;
    
    % Decimal places to round down to
    decimalPlaces = 2;
    
    % Call the roundDown function
    roundedNumber = roundDown(number, decimalPlaces);
    
    % Display the result
    % disp(roundedNumber);
    
    function roundedNumber = roundDown(number, decimalPlaces)
        % Calculate the factor to scale the number
        scaleFactor = 10^decimalPlaces;
        
        % Scale the number and round down using floor
        scaledNumber = floor(number * scaleFactor);
        
        % Scale back to the original number
        roundedNumber = scaledNumber / scaleFactor;
    end
end