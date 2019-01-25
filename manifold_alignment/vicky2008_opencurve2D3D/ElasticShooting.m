function q2n = ElasticShooting(q1,v)

d = sqrt(InnerProd_Q(v,v));
%  eps=.1;
if d < 0.00001
    q2n = q1;
else
%     q2{1} = cos(eps*d)*q1 + (sin(eps*d)/d)*v;
% % % %     q2{1} = ProjectC(q2{1});
%    v = Parallel_Transport_C(v,q1,q2{1});
%      d = sqrt(InnerProd_Q(v,v));
%    for j=2:10
%          q2{j} = cos(eps*d)*q2{j-1} + (sin(eps*d)/d)*v;
% % %        q2{j} = ProjectC(q2{j});
%        v = Parallel_Transport_C(v,q2{j-1},q2{j});
%         d = sqrt(InnerProd_Q(v,v));
%      end
%      q2n = q2{10};
%  end

q2n=q1+v;
end