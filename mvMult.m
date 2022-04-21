function [w] = mvMult(v)
    global Ua Ub Sa Sb Va Vb;
    w = Vb' * diag(v) * Va;
    w = Sb * w * Sa;
    w = sum((Ub * w).*Ua'.',2);
end

