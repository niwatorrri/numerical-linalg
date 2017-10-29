for k = 5:20
    Hk = hilb(k);
    fprintf('%d & %e & %e\\\\\n', k, ...
            norm(Hk, inf) * invnorm(Hk), cond(Hk, inf));
end