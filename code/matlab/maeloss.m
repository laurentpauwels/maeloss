function mae = maeloss(a_init,nu_it)
%a = [exp(a_init);1]; % compute weights that are between 0 and 1 and sum to 1
%a = a/sum(a);
%T = size(nu_it,1);
[T,p] = size(nu_it);
%a = zeros(p,1);
a(1:p-1,1) = a_init;
a(p,1) = 1 - sum(a_init);
nu_c = nu_it*a;
mae = sum(abs(nu_c))/T;
end