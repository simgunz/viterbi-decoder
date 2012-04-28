function [ u_output ] = convolutionalTxSystem( u_input, g1, g2, mex, gammaDB )
%CONVOLUTIONALTXSYSTEM Summary of this function goes here
%   Detailed explanation goes here

nu = length(g1) - 1;             % Number of memory elements
ns = 2^nu;

% Y(u,s,k) Output map, store the vector [y1, y2]
% in position (u,s,:) associated to the input u and the state s
% S(u,s) State update map (matrix), store the state sNext
% in position (u,s) associated to the input u and the state s
% N(s,k,i) Neighbors map, store the k-th neighbor [u, sN] associated to
% the state s in position (s,k,:) 

% For Matlab all the values of input and states are shifted by 1
% To match matlab index convention the range of stored
% state is 1-2^nu and the range of stored input is 1-2
% To match C index convention the range of stored input is 0-1
% and the range of stored state is 0-(2^nu)-1

[Y,S,N] = buildMaps(g1,g2);
    
u = [u_input zeros(1,nu)];         % Input extended with zero padding

if ~mex              % Match the correct index system
    S = S + 1;
    N = N + 1; 
    u = u + 1;                
end

%%%%%%% ENCODER %%%%%%%

% Output vector, each column is [y1; y2]

if mex
    y = vitencoder(u,Y,S);
else
    s = 1;                  % Initial state
    y = zeros(2,mu + nu);
    for k=1:(mu + nu)
        y(:,k) = Y(u(k),s,:);
        s = S(u(k),s);
    end
end


%%%%%%%%% P/S %%%%%%%%%

d = y(:)';


%%% BPAM MODULATOR %%%%

sTx = d;
sTx(sTx==0) = -1;    % Map 0 to -1


%%%%%%% CHANNEL %%%%%%%

r = awgn(sTx,gammaDB);                           % Received signal


%%%%%%%%% S/P %%%%%%%%%

r = reshape(r,2,length(r)/2);


%%% VITERBI DECODER %%%

if mex
    u_output = vitdecoder(r,Y,S,N);
else
    gamma = ones(1,ns)*(-inf);    % Metrics column vector
    gamma(1) = 0;
    gammanew = gamma;
    survivors = zeros(2,mu+nu,ns);       % Survivors array of vectors        
    M = zeros(2);
    for i=1:(mu + nu)                   % For every output [y1; y2]
        for j=1:ns                    % For every possible state                
            M(:,1) = squeeze(Y(N(j,1,1),N(j,1,2),:));
            M(:,2) = squeeze(Y(N(j,2,1),N(j,2,2),:));
            M(M==0) = -1;
            tempgamma = r(:,i)'*M + gamma(N(j,1,2):N(j,2,2));
            if tempgamma(1) > tempgamma(2) 
                bestK = 1;
            else
                bestK = 2;
            end
            gammanew(j) = tempgamma(bestK);
            survivors(:,i,j) = squeeze(N(j,bestK,:));                
        end
        gamma = gammanew - max(gammanew);            
    end
    u_output = zeros(1,mu+nu);
    ss = 1;    
    for k=fliplr(1:(nu+mu))        
        u_output(k) = survivors(1,k,ss);
        ss = survivors(2,k,ss);
    end
    u_output = u_output(1:mu) - 1;      % Eventually translate from matlab index to real value
end   
    
end

