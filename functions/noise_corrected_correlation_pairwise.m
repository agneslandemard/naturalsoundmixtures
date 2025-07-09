function c = noise_corrected_correlation_pairwise(r1,r2,p1,p2,z_correct)
% if z-correct (default: false), will go in z-space to do mean of both
% correlations
if nargin < 5
    z_correct = true;
end
 
c1 = nancorr_pairwise(r1,p2);
c2 = nancorr_pairwise(r2,p1);

if z_correct
    n = tanh(0.5*(atanh(c1)+atanh(c2)));
else
    n = 0.5*(c1+c2);
end

d = sqrt(nancorr_pairwise(r1,r2).*nancorr_pairwise(p1,p2));
c = n./d;

end