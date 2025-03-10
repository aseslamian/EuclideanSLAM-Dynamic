function  [mismatch] = detectMLM(mat1, mat2,sigma)

common = intersect(mat1(:,1:2),mat2(:,1:2),'rows');
c1 = zeros(size(common, 1), size(mat1, 2));
c2 = zeros(size(common, 1), size(mat2, 2));
tolerance = 1e-6;

counter1 = 0;
counter2 = 0;

for i = 1:size(common,1)
    row = common(i,:);
    
    for j = 1:size(mat1,1)
        if all(abs(row - mat1(j,1:2)) < tolerance)
            counter1 = counter1 + 1;
            c1(counter1,:) = mat1(j,:);
        end
    end
    
    for j = 1:size(mat2,1)
        if all(abs(row - mat2(j,1:2)) < tolerance)
            counter2 = counter2 + 1;
            c2(counter2,:) = mat2(j,:);
        end
    end
end

c1 = c1(1:counter1, :);
c2 = c2(1:counter2, :);

Diff =abs(c1 - c2);

Diff = round(Diff/sigma^2)*sigma^2;
[row,~] = find(Diff);
vec = mat1(row,1:2);
mismatch = mode(vec(:));

end