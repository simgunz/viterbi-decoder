clc;
%clear all;

%%%%%%% TUNABLE PARAMETERS

g1 = [1 0 1];
g2 = [1 1 1];

mu = 1000000;                     % Input length
nu = length(g1) - 1;            % Number of memory elements

Psymbol = 0.08;     % Psymbol = Q(1/std_deviation_w)
                    % Percentage of erroneous transmitted BPAM symbols (symbols received with
                    % opposite sign)

mexEnabled = 1;     % Enable or disable c implementation of the encoder/decoder

iter=30;
% Y(u,s) Output map (cell), store the column vector [y1; y2]
% in position (u,s) associated to the input u and the state s
% S(u,s) State update map (matrix), store the state sNext
% in position (u,s) associated to the input u and the state s
% N(s) Neighbors map (cell), store the array of neighbors [u; s]
% in position (s) associated to the state s

% It is possible to create the maps to match the array index system
% of Matlab or C. For Matlab all the values of input and states are shifted by 1
% To match matlab index convention the range of stored
% state is 1-2^nu and the range of stored input is 1-2
% To match C index convention the range of stored input is 0-1
% and the range of stored state is 0-(2^nu)-1

if mexEnabled
    startIndex = 0;     % First array index is 0 (C)
else
    startIndex = 1;     % First array index is 1 (Matlab)
end


[Y,S,N] = buildMaps(g1,g2,startIndex);
tic
for i=1:iter

    %u_input = [0 1 1 1 0 1 0 0 0];
    u_input = round(rand(1,mu));       % Input sequence

    u = [u_input zeros(1,nu)];         % Input extended with zero padding
    u = u + startIndex;                % Match the correct index system

    %%%%%%% ENCODER %%%%%%%

    % Output vector, each column is [y1; y2]


    s = startIndex;                  % Initial state
    if mexEnabled
        y = vitencoder(u,Y,S);
    else
        y = zeros(2,mu + nu);
        for k=1:(mu + nu)
            y(:,k) = Y{u(k),s};
            s = S(u(k),s);
        end
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
    %% ELIMINA RIDIMENSIONAMENTO ARRAY
%     if mexEnabled
%         u_output = vitdecoder(r,S,Y);
%     else
%         gamma = ones(1,2^nu)*(-inf);    % Metrics column vector
%         gamma(1) = 0;
%         gammanew = gamma;
%         survivors = cell(2^nu,1);       % Survivors array of vectors
%         survivorsNew = cell(2^nu,1);
%         for i=1:(mu + nu)                   % For every output [y1; y2]
%             for j=1:2^nu                    % For every possible state
%
%                 maxgamma = -inf;
%                 for k=1:length(N{j})            % For every neighbors of the current state j
%                     neighInput = N{j}(1,k);
%                     neighState = N{j}(2,k);
%                     M = Y{neighInput,neighState};
%                     M(M==0) = -1;
%                     tempgamma = gamma(neighState) + sum(r(:,i)'*M);
%                     if tempgamma > maxgamma
%                         goodNeighState = neighState;
%                         newSurv = N{j}(1,k);
%                         maxgamma = tempgamma;
%                     end
%                 end
%                 survivorsNew{j} = [survivors{goodNeighState} newSurv];
%                 gammanew(j) = maxgamma;
%             end
%             gamma = gammanew - max(gammanew);
%             survivors = survivorsNew;
%         end
%
%         u_output = survivors{1}(1:mu) - 1;      % Translate from matlab index to real value
%     end
%
%     residualErr = sum(u_input ~= u_output);

end
toc
disp(['Symbols errors = ',num2str(symbolErr)]);
disp(['Residual errors = ',num2str(residualErr)]);

