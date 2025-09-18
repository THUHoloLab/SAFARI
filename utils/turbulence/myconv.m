function C = myconv(A, B, delta)
% function C = myconv(A, B, delta)
    N = length(A);
    C = ift(ft(A, delta) .* ft(B, delta), 1/(N*delta));