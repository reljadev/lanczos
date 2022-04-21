function [w] = mvMult_transpose(x)
    global Ua Ub Sa Sb Va Vb;
    w = Ub' * diag(x) * Ua;
    w = Sb * w * Sa;
    w = sum((Vb * w).*Va'.',2);
end
