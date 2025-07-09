function new_struct = correct_struct(old_struct, mask, remove)
fn = fieldnames(old_struct);
new_struct = struct();
for f = 1:length(fn)
    if isnumeric(old_struct.(fn{f}))
        sz = size(old_struct.(fn{f}));
        mask_d = repmat(mask, [1 sz(2:end)]);
        if remove 
            new = reshape(old_struct.(fn{f})(mask_d), [sum(mask) sz(2:end)]);
        else
            new = old_struct.(fn{f});
            new(~ mask_d) = nan;
  
        end
        new_struct.(fn{f}) = new;
    end
end
end