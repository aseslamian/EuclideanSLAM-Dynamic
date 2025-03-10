function [z, index] = removedynamic(mismatch, z)

    if ~isempty(mismatch) && size(z,1) > 2
        mismatch = unique(mismatch);
        index = find(ismember(z(:,1), mismatch));
        z(index,:) = [];
    else
        index = [];
    end

end
