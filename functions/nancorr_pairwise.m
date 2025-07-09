function C = nancorr_pairwise(A,B,corr_type)
if nargin < 3
    corr_type = 'pearson';
end
% computes correlation between matching columns of A and B, ignoring nans

% compute nan mask to apply to each matrix
nan_mask = isnan(A) | isnan(B);
A(nan_mask) = nan;
B(nan_mask) = nan;

if strcmp(corr_type,'spearman')
    [~, At] = sort(A, 1);
    [~, Bt] = sort(B, 1);
    % put nans back in
    for col = 1: size(A,2)
        At(isnan(A(:, col)), col) = nan;
        Bt(isnan(B(:, col)), col) = nan;
    end

else
    At = A;
    Bt = B;
end
An = bsxfun(@minus,At,nanmean(At,1)); %%% zero-mean
Bn = bsxfun(@minus,Bt,nanmean(Bt,1)); %%% zero-mean
An = bsxfun(@times,An,1./sqrt(nansum(An.^2,1))); %% L2-normalization
Bn = bsxfun(@times,Bn,1./sqrt(nansum(Bn.^2,1))); %% L2-normalization
C = nansum(An.*Bn,1);

end


