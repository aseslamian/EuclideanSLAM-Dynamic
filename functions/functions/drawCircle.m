function c = drawCircle(x_center, y_center, radius, varargin)
    %DRAWTRANSPARENTCIRCLE Draws a transparent circle on a MATLAB figure.
    %
    % Inputs:
    %   - x_center, y_center: Coordinates of the circle's center.
    %   - radius: Radius of the circle.
    %   - transparency: Value between 0 (completely transparent) and 1 (completely opaque).
    %   - varargin: Any additional properties for the circle.
    %
    % Outputs:
    %   - h: Handle to the circle (patch object).

    % Create circle coordinates
    theta = linspace(0, 2*pi, 100); % Create 100 points to form the circle
    x = x_center + radius * cos(theta);
    y = y_center + radius * sin(theta);

    % Use patch to draw a filled transparent circle
    c = patch(x, y, 'yellow',varargin{:});
end
