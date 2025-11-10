function X0 = generateFeasibleRepairSchedule(numNodes, T, C)

% GENERATEFEASIBLEREPAIRSCHEDULE returns a feasible repair schedule.
%
%   X0 = generateFeasibleRepairSchedule(numNodes, T, C)
%
% For each time step, this helper function randomly selects an integer total
% r (between 0 and C) for repair crews and randomly distributes these r crews
% among the nodes, ensuring sum(X0(:,t)) <= C.
%
% Output:
%   X0 - A (numNodes x T) matrix with feasible allocations.

    X0 = zeros(numNodes, T);
    for t = 1:T
        totalCrews = randi([0, C]);  % Total crews for time step t.
        allocation = zeros(numNodes, 1);
        for k = 1:totalCrews
            idx = randi(numNodes);
            allocation(idx) = allocation(idx) + 1;
        end
        X0(:, t) = allocation;
    end

% GENERATEFEASIBLEREPAIRSCHEDULE generates a feasible repair schedule.
%
% Inputs:
%   numNodes - Number of nodes.
%   T        - Number of time steps.
%   C        - Global repair crew budget per time step.
%
% Output:
%   X0       - A (numNodes x T) matrix such that, for every time step t,
%              sum(X0(:,t)) <= C.
%
% This function works as follows:
%   For each time step, we randomly decide on a total number of crews r (an 
%   integer between 0 and C). We then distribute these r crews randomly over 
%   the nodes.

    % X0 = zeros(numNodes, T);
    % for t = 1:T
    %     % Randomly choose a total number of repair crews for this time step.
    %     totalCrews = randi([0, C]);  % Total crews allocated this time step.
    % 
    %     % If totalCrews is 0, then no node is allocated a crew.
    %     if totalCrews == 0
    %         X0(:,t) = zeros(numNodes,1);
    %     else
    %         % Generate a random partition of totalCrews among numNodes.
    %         % One simple method: assign each of the r units to a random node.
    %         allocation = zeros(numNodes, 1);
    %         for k = 1:totalCrews
    %             idx = randi(numNodes);
    %             allocation(idx) = allocation(idx) + 1;
    %         end
    %         X0(:,t) = allocation;
    %     end
    % end
end
