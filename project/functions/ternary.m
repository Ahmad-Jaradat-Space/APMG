function out = ternary(cond, val_true, val_false)
    if cond
        out = val_true;
    else
        out = val_false;
    end
end