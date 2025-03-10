function theta = angleBetweenPoints(A, B, C)
    % Define vectors BA and BC
    BA = [A(1) - B(1), A(2) - B(2)];
    BC = [C(1) - B(1), C(2) - B(2)];

    % Compute the dot product of BA and BC
    dotProduct = dot(BA, BC);

    % Compute the magnitudes of BA and BC
    magBA = norm(BA);
    magBC = norm(BC);

    % Compute the cosine of the angle
    cosTheta = dotProduct / (magBA * magBC);

    % Calculate the angle in radians
    theta = acos(cosTheta);
end
