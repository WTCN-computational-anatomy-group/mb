function to = CopyFields(from, to)
fn = fieldnames(from);
for i=1:numel(fn)
    to.(fn{i}) = from.(fn{i});
end
end
%==========================================================================