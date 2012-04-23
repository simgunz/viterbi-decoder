function [ Y, S, N] = buildMaps( g1, g2, startIndex )
%NEIGHBORS Summary of this function goes here
%   Detailed explanation goes here

dim = length(g1);
nu = dim-1;

A = [g1(1); g2(1)];
B = [g1(2:dim); g2(2:dim)];

C = zeros(dim-1,1);
C(1) = 1;
D = diag(ones(1,dim-2),-1);

Y = cell(2,2^nu);
S = zeros(2,2^nu);
N = cell(2^nu,1);

offset = ~startIndex;

for s=1:2^nu
    for u=1:2
        Y{u,s} = mod(A*(u-1) + B*toBin(s-1,nu),2);
        S(u,s) = toDec( mod(C*(u-1) + D*toBin(s-1,nu),2) ) + startIndex;
        N{S(u,s) + offset} = [N{S(u,s) + offset} [u - offset; s - offset]];
    end
end

end

function [ binS ] = toBin ( s , nu)
    binS = fliplr(de2bi(s,nu))';
end

function [ decS ] = toDec ( s )
    decS = bi2de(fliplr(s'));
end