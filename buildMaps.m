function [ Y, S, N] = buildMaps( g1, g2 )
%BUILDMAS Build the output map Y, the state update map S and the neighbors map N
% Y(u,s,k) Output map, store the vector [y1, y2]
% in position (u,s,:) associated to the input u and the state s
% S(u,s) State update map (matrix), store the state sNext
% in position (u,s) associated to the input u and the state s
% N(s,k,i) Neighbors map, store the k-th neighbor [u, sN] associated to
% the state s in position (s,k,:)

% To match matlab index convention the range of stored
% state is 1-2^nu and the range of stored input is 1-2
% To match C index convention the range of stored input is 0-1
% and the range of stored state is 0-(2^nu)-1
% The maps are produced in C format. To match matlab index system it is
% necessary to add one to every elements of S,N matrix

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

neighIndex = ones(1,ns); % Position in the N table for the next neighbor
for s=1:ns
    for u=1:2
        Y(u,s,:) = mod(A*(u-1) + B*toBin(s-1,nu),2);
        S(u,s) = toDec(mod(C*(u-1) + D*toBin(s-1,nu),2));
        N(S(u,s)+1,neighIndex(S(u,s)+1),:) = [u-1, s-1];
        neighIndex(S(u,s)+1) = neighIndex(S(u,s)+1) + 1;
    end
end

end

function [ binS ] = toBin ( s , nu)
    binS = fliplr(de2bi(s,nu))';
end

function [ decS ] = toDec ( s )
    decS = bi2de(fliplr(s'));
end