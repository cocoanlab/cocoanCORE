function out = markov_chain(states_input)

% Do markov chain analysis using one time transition using state time
% series data.
%
% :Usage:
% ::
%
%    out = markov_chain(states_time_series)
%
% :Inputs:
%
%   **states_time_series:**
%      the number of states, interger, column vector
%      e.g., states = [1 2 3 1 2 2 2 2 3 3 3 3 1 1 1]';
%
%      If the input is a matrix with multiple columns, it is doing the
%      markov chain analysis on each column.
%
% :Outputs:
%
%   **out.trans_mat**   transition matrix (counts)
%
%   **out.trans_prob**  transition probability
%                       (rows are states at t-1, and columns are states at t)
%
%   **out.steady_state_prob**   steady state probability calculated from
%                               transition probability
%
%   **out.frequency_prob**   frequency-based calculation of the probablity
%                            for each state
%

n_col = size(states_input,2);

for j = 1:n_col
    
    states = states_input(:,j);
    
    n_state = numel(unique(states));
    trans_idx = [states(1:end-1) states(2:end)];
    out{j}.trans_mat = zeros(n_state,n_state);
    
    for i = 1:size(trans_idx,1)
        out{j}.trans_mat(trans_idx(i,1), trans_idx(i,2)) = out{j}.trans_mat(trans_idx(i,1), trans_idx(i,2))+1;
    end
    
    out{j}.trans_prob = out{j}.trans_mat./repmat(sum(out{j}.trans_mat,2),1,n_state);
    
    if any(sum(isnan(out{j}.trans_prob),2)==n_state)
        out{j}.trans_prob(sum(isnan(out{j}.trans_prob),2)==n_state,:) = repmat(1/n_state, 1, n_state);
    end
    
    % steady-state probabilities
    out{j}.steady_state_prob = out{j}.trans_prob^1000000;   % trans_prob's power
    
    for i = 1:size(out{j}.steady_state_prob,1)
        wh_eq(i,1) = sum(out{j}.steady_state_prob(1,:)-out{j}.steady_state_prob(i,:))< .0000001;
    end
    
    if any(~wh_eq)
        error('Steady state prob is not converged. Something is wrong.');
    else
        out{j}.steady_state_prob = out{j}.steady_state_prob(end,:);
    end
    
    out{j}.frequency_prob = mean(out{j}.trans_prob); % trans_prob's mean
    
end

end