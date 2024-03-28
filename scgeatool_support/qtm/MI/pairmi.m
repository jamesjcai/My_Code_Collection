function M = pairmi(data)

n = size(data, 1);

M = zeros(n);

for j = 1:n
    for k = j:n
        M(j,k) = mutualinfo(data(j,:), data(k,:));
        M(k,j) = M(j,k);
    end
end

% for j = 1:n
%     for k = j:n
%         M(k,j) = M(j,k);
%     end
% end
