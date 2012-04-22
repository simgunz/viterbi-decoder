function [ Y, S, N] = buildMaps( g1, g2 )
%NEIGHBORS Summary of this function goes here
%   Detailed explanation goes here

dim = length(g1);
nu = dim-1;

A = [g1(1); g2(1)];
B = [g1(2:dim); g2(2:dim)];

C = zeros(dim-1,1);
C(1) = 1;
D = diag(ones(1,dim-2),-1);

% Y = @(u,s) xor(A*u, B*s);
% S = @(u,s) xor(C*u, D*s);
Y = @(u,s) mod(A*u + B*s,2);
S = @(u,s) mod(C*u + D*s,2);

N=cell(2^nu,1);
for i=1:2^nu
    for j=1:2
        s = S(j-1,fliplr(de2bi(i-1,nu))');
        N{bi2de(fliplr(s'))+1} = [N{bi2de(fliplr(s'))+1} [j-1; i-1]];
    end
end
end


