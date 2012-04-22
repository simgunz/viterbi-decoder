clc;
clear all;

%%%%%%% TUNABLE PARAMETERS

g1 = [1 0 1];
g2 = [1 1 1];

%u_input = [0 1 1 1 0 1 0 0 0];
u_input = round(rand(1,10000));    % Input sequence

Psymbol = 0.08;     % Psymbol = Q(1/std_deviation_w)
                    % Percentage of erroneous transmitted BPAM symbols (symbols received with
                    % opposite sign)

[Y,S,N] = buildMaps(g1,g2);     % Y(u,s) Output map (cell), store the column vector [y1; y2]
                                % in position (u,s) associated to the input u and the state s
                                % S(u,s) State update map (matrix), store the state sNext
                                % in position (u,s) associated to the input u and the state s
                                % N(s) Neighbors map (cell), store the array of neighbors [u; s]
                                % in position (s) associated to the state s
                                % To match matlab index convention the range of stored 
                                % state is 1-2^nu and the range of stored input is 1-2
                                % The real range of input is 0-1 and the real range of state
                                % is 0-(2^nu)-1
                                

mu = length(u_input);           % Input length
nu = length(g1) - 1;            % Number of memory elements
q = 2^nu;                       % Number of possible states

u = [u_input zeros(1,nu)];         % Input extended with zero padding
u = u + 1;                         % Translation to match matlab indexing
    
s = 1;                  % Initial state
y = zeros(2,mu + nu);   % Output vector, each column is [y1; y2]


%%%%%%% ENCODER %%%%%%%

for k=1:(mu + nu)
    y(:,k) = Y{u(k),s};
    s = S(u(k),s);
end


%%%%%%%%% P/S %%%%%%%%%

d = y(:)';


%%% BPAM MODULATOR %%%%

sTx = d;
sTx(sTx==0) = -1;    % Map 0 to -1


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
survivorsNew = cell(2^nu,1);
for i=1:(mu + nu)                   % For every output [y1; y2]
    for j=1:2^nu                    % For every possible state
                                    
        maxgamma = -inf;        
        for k=1:length(N{j})            % For every neighbors of the current state j
            neighInput = N{j}(1,k);
            neighState = N{j}(2,k);         
            M = Y{neighInput,neighState};
            M(M==0) = -1;            
            tempgamma = gamma(neighState) + sum(r(:,i)'*M);
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

u_output = survivors{1}(1:mu) - 1;      % Translate from matlab index to real value

residualErr = sum(u_input ~= u_output);

disp(['Symbols errors = ',num2str(symbolErr)]);
disp(['Residual errors = ',num2str(residualErr)]);

