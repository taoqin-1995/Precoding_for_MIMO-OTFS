function power_allocation = water_filling(channel_gains, total_power)
    % Water-filling algorithm to allocate power optimally across subchannels
    %
    % Inputs:
    % - channel_gains: vector of channel gains (squared singular values or SNR of each subchannel)
    % - total_power: total available power to distribute
    %
    % Output:
    % - power_allocation: vector of allocated power for each subchannel

    % Sort channel gains in descending order
    [sorted_gains, indices] = sort(channel_gains, 'descend');
    
    % Initialize parameters
    num_channels = length(channel_gains);
    power_allocation = zeros(num_channels, 1);
    
    % Compute inverse of channel gains for water level calculations
    inverse_gains = 1 ./ sorted_gains;
    
    % Initialize water level and remaining power
    water_level = 0;
    remaining_power = total_power;
    
    for i = 1:num_channels
        % Update water level based on current subchannel
        water_level = (remaining_power + sum(inverse_gains(1:i))) / i;
        
        % Check if current water level can be achieved
        if water_level > inverse_gains(i)
            % Allocate power to current subchannel set
            power_allocation(indices(1:i)) = water_level - inverse_gains(1:i);
            % Update remaining power
            remaining_power = total_power - sum(power_allocation);
        else
            break;
        end
    end
end
