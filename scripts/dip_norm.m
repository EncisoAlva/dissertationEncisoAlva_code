function normJ = dip_norm( J )
% In the case of unconstrained dpoles, the magnitude of J needs to be
% extracted. This happens very often so it is made into a function.
nDips = size(J,1)/3;
normJ = zeros(nDips, size(J,2));
for ii = 1:nDips
  normJ(ii,:) = vecnorm( J(3*(ii-1)+[1,2,3],:), 2, 1 );
end

end