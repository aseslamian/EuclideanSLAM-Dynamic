function c = estDist(a, b, C)    
    % Use the Law of Cosines to compute the third side
    c = sqrt(a.^2 + b.^2 - 2*a.*b.*cos(C));
end
