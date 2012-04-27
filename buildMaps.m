function [ Y, S, N] = buildMaps( g1, g2 )
%NEIGHBORS Summary of this function goes here
%   Detailed explanation goes here

dim = length(g1);
nu = dim-1;             % Number of memory elements
ns = 2^nu;              % Number of possible states
A = [g1(1); g2(1)];
B = [g1(2:dim); g2(2:dim)];

C = zeros(nu,1);
C(1) = 1;
D = diag(ones(1,nu-1),-1);

Y = zeros(2,ns,2);          
S = zeros(2,ns);
N = zeros(ns,2,2);

predIndex = ones(1,ns); % Position in the N table for the next predecessor
for s=1:ns
    for u=1:2
        Y(u,s,:) = mod(A*(u-1) + B*toBin(s-1,nu),2);
        S(u,s) = toDec(mod(C*(u-1) + D*toBin(s-1,nu),2));
        N(S(u,s)+1,predIndex(S(u,s)+1),:) = [u-1, s-1];
        predIndex(S(u,s)+1) = predIndex(S(u,s)+1) + 1;
    end
end

end

function [ binS ] = toBin ( s , nu)
    binS = fliplr(de2bi(s,nu))';
end

function [ decS ] = toDec ( s )
    decS = bi2de(fliplr(s'));
end