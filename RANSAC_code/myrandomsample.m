function y = myrandomsample(n, k)
%MYRANDOMSAMPLE Random sample, with or without replacement.
%   Y = MYRANDOMSAMPLE(N,K) returns Y as a column vector of K values sampled
%   uniformly at random, without replacement, from the integers 1:N.

    r = (n-1)*rand(k, 1);
    y = 1+round(r);
%end