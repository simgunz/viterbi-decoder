function [ u_output ] = convolutionalTxSystem( u_input, g1, g2, enableMex, gammaDB, finiteTable)
%CONVOLUTIONALTXSYSTEM Simulates a full transmission system with convolutional encoding

dim = length(de2bi(oct2dec(max(g1,g2))));
g1 = fliplr(de2bi(oct2dec(g1),dim));        % Convert the genarators from octal
g2 = fliplr(de2bi(oct2dec(g2),dim));        % to binary vectors

mu = length(u_input);            % Input length
nu = length(g1) - 1;             % Number of memory elements
ns = 2^nu;                       % Number of possible states

[Y,S,N] = buildMaps(g1,g2);

u = [u_input zeros(1,nu)];       % Input extended with zero padding

if ~enableMex              % Match the correct index system
    S = S + 1;
    N = N + 1;
    u = u + 1;
end


%%%%%%% ENCODER %%%%%%%

% y Output vector, each column is [y1; y2]

if enableMex
    y = vitencoder(u,Y,S);      % Decode with mex
else
    s = 1;                      % Initial state
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
sTx(sTx==0) = -1;    % Map 0 to -1 and produce the transmitted signal


%%%%%%% CHANNEL %%%%%%%

r = awgn(sTx,gammaDB);          % Received signal


%%%%%%%%% S/P %%%%%%%%%

r = reshape(r,2,length(r)/2);


%%% VITERBI DECODER %%%

if enableMex
    if ~finiteTable
        u_output = vitdecoder(r,Y,S,N);
    else
        u_output = vitdecoderFiniteTable(r,Y,S,N);
    end
else
    gamma = ones(1,ns)*(-inf);          % Metrics vector
    gamma(1) = 0;
    gammanew = gamma;
    survivors = zeros(2,mu+nu,ns);      % Survivors table
    M = zeros(2);
    for i=1:(mu + nu)                   % For every output [y1; y2]
        for j=1:ns                      % For every possible state
            M(:,1) = squeeze(Y(N(j,1,1),N(j,1,2),:));   % Get the output associated to
            M(:,2) = squeeze(Y(N(j,2,1),N(j,2,2),:));   % the neighbors of j
            M(M==0) = -1;                               % and modulate the output
            tempgamma = r(:,i)'*M + gamma([N(j,1,2),N(j,2,2)]);
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
    u_output = u_output(1:mu) - 1;      % Translate back from matlab index
                                        % system to real value
end

end

