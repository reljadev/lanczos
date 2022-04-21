function [w] = mvMult_times(v)
    w = mvMult(mvMult_transpose(v));
end

