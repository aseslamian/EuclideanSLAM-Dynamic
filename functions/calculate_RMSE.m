function RMSE = calculate_RMSE(x, y)

    if any(size(x) ~= size(y)) || size(x, 1) ~= 2
        error('Both vectors must have size (2, N)');
    end

% Calculate the squared differences
squared_diffs = (x(1,:) - y(1,:)).^2 + (x(2,:) - y(2,:)).^2;

% Calculate the mean
RMSE = sqrt(mean(squared_diffs));

end
