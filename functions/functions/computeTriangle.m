function [xs, ys] = computeTriangle(pt1, pt2, size)
    angle = atan2(pt2(2) - pt1(2), pt2(1) - pt1(1));
    
    xs = [pt2(1) + size * cos(angle)
          pt2(1) + size/2 * cos(angle - 2*pi/3)
          pt2(1) + size/2 * cos(angle + 2*pi/3)];
          
    ys = [pt2(2) + size * sin(angle)
          pt2(2) + size/2 * sin(angle - 2*pi/3)
          pt2(2) + size/2 * sin(angle + 2*pi/3)];
end