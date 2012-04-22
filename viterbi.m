clc;
clear all;

%%%%%%% TUNABLE PARAMETERS

g1 = [1 0 1];
g2 = [1 1 1];

%u_input = [0 1 1 1 0 1 0 0 0];
u_input = round(rand(1,1000));    % Input sequence

Psymbol = 0.08;     % Psymbol = Q(1/std_deviation_w)
                    % Percentage of erroneous transmitted BPAM symbols (symbols received with
                    % opposite sign)

[Y,S,N] = buildMaps(g1,g2);     % Y(u,s) Output function, return the column vector [y1; y2]
                                % associated to the input u and the state s
                                % S(u,s) State update function, return the column vector
                                % [s1; ... s_nu] associated to the input u and the state s
                                % N Neighbors map, array of vectors, each associated to the state j.
                                % Every vector contains all the qbors of j in the form [u; s]
                                % The states go from 0 to (2^nu)-1

mu = length(u_input);         % Input length
nu = length(g1)-1;              % Number of memory elements
q = 2^nu;                       % Number of possible states

u = [u_input zeros(1,nu)];         % Input extended with zero padding
    
s = zeros(nu,1);        % State (column) vector
y = zeros(2,mu + nu);   % Output vector, each column is [y1; y2]


%%%%%%% ENCODER %%%%%%%

for k=1:(mu + nu)
    y(:,k) = Y(u(k),s);
    s = S(u(k),s);
end


%%%%%%%%% P/S %%%%%%%%%

d = y(:)';


%%% BPAM MODULATOR %%%%

sTx = d;
index = find(sTx-1);            % Find the elements equal to zero
sTx(index) = sTx(index) - 1;    % Set them to -1


%%%%%%% CHANNEL %%%%%%%

std_deviation_w = 1/(qfuncinv(Psymbol));            % Noise standard deviation
r = sTx + randn(1,length(sTx)) * std_deviation_w;   % Received signal
symbolErr = sum(sign(sTx) ~= sign(r));              % Number of erroneous tx symbols

%%%%%%%%% S/P %%%%%%%%%

r = reshape(r,2,length(r)/2);


%%% VITERBI DECODER %%%

gamma = ones(1,2^nu)*(-inf);    % Metrics column vector
gamma(1) = 0;
gammanew = gamma;               
survivors = cell(2^nu,1);       % Survivors array of vectors

for i=1:(mu + nu)                   % For every output [y1; y2]
    for j=1:2^nu                    % For every possible state, possible states go from 1 to 2^nu
                                    % in this way they can be used as index to access vectors
        maxgamma = -inf;        
        for k=1:length(N{j})            % For every neighbors of the current state j
            neighInput = N{j}(1,k);
            neighState = N{j}(2,k) + 1;         % Convert state from range 0 - (2^nu)-1 to range 1-2^nu
            neighStateBin = fliplr(de2bi(neighState - 1,nu))';  % Convert the state in binary form
            
            tempgamma = gamma(neighState) + sum(r(:,i)'*M(Y(neighInput,neighStateBin)));
            if tempgamma > maxgamma
                goodNeighState = neighState;
                newSurv = N{j}(1,k);
                maxgamma = tempgamma;
            end
        end
        survivorsNew{j} = [survivors{goodNeighState} newSurv];
        gammanew(j) = maxgamma;
    end
    gamma = gammanew - max(gammanew);
    survivors = survivorsNew;
end

u_output = survivors{1}(1:mu);

residualErr = sum(u_input ~= u_output);

disp(['Symbols errors = ',num2str(symbolErr)]);
disp(['Residual errors = ',num2str(residualErr)]);

