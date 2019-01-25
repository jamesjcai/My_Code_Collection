function p = q_to_curve(q)

[n,T] = size(q);

for i = 1:T
    qnorm(i) = norm(q(:,i));
end

for i = 1:n
    p(i,:) = [ cumtrapz( q(i,:).*qnorm )/(T) ] ;
end