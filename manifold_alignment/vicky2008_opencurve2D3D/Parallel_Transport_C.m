function w_new = Parallel_Transport_C(w,q1,q2)

lw = sqrt(InnerProd_Q(w,w));
if(lw < 0.0001)
    w_new = w;
else
    w_new = Project_Tangent(w,q2);
    w_new = w_new*lw/sqrt(InnerProd_Q(w_new,w_new));
end

return;
