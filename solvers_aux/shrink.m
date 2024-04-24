function vec_res = shrink( num, vec )
% shrinking operator

% patch
num2 = abs(num);

% operator per se
vec_res = zeros(size(vec));
vec_res( vec> num2 ) = vec(vec> num2) - num2;
vec_res( vec<-num2 ) = vec(vec<-num2) + num2;

end