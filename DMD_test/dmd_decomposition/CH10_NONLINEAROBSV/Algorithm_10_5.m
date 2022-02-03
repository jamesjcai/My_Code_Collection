q2=[q; (q.*abs(q).^2)];
y0 = Phi2\q2;

q_modes = zeros(r,length(t));  
for iter = 1:length(t)
q_modes(:,iter) =(y0.*exp(omega*(t(iter))));
end
q_dmd2 = Phi2*q_modes;  

q_dmd = q_dmd2(1:n,:); %% Koopman approximation
Phi = Phi2(1:n,:); %% Koopman eigenfunctions