clear; clc; close all;

O = 1:5;
D = 6:10;        
J_out = 11:20;   
J_in = 21:25;    

p_int = 15; 
n = 25;     

s = []; t = []; base_weights = [];

for i = 1:5
    s = [s; O(i)]; t = [t; J_out(2*i-1)]; base_weights = [base_weights; 1.0];
end

for i = 1:5
    s = [s; J_out(2*i)]; t = [t; D(i)]; base_weights = [base_weights; 1.0];
end

for k = 1:10
    next_k = mod(k, 10) + 1;
    s = [s; J_out(k); J_out(next_k)]; base_weights = [base_weights; 3.0];
    t = [t; J_out(next_k); J_out(k)]; base_weights = [base_weights; 3.0];
end

for i = 1:5
    next_i = mod(i, 5) + 1;
    s = [s; J_in(i)]; t = [t; J_in(next_i)]; base_weights = [base_weights; 0.5];
end

for i = 1:5
    s = [s; J_out(2*i-1); J_in(i)]; t = [t; J_in(i); J_out(2*i-1)];
    base_weights = [base_weights; 1.0; 1.0];
    s = [s; J_out(2*i); J_in(i)];   t = [t; J_in(i); J_out(2*i)];
    base_weights = [base_weights; 1.0; 1.0];
end

m = length(s); 

A = [];
RoutePaths = cell(n, 1);
rank_current = 0;
route_idx = 1;
fprintf('Generating 25 LOGICAL linearly independent routes...\n');

base_weights = ones(m, 1);
for e = 1:m
    if s(e) > 20 && t(e) > 20
        base_weights(e) = 0.5; 
    end
end

for i = 1:length(O)
    for j = 1:length(D)
        orig = O(i); dest = D(j);
        weights = base_weights; 
        
        for attempt = 1:20
            G = digraph(s, t, weights);
            path = shortestpath(G, orig, dest, 'Method', 'positive');
            if isempty(path), break; end
            
            col = zeros(m, 1);
            for step = 1:(length(path)-1)
                u = path(step); v = path(step+1);
                idx = find(s == u & t == v, 1);
                col(idx) = 1;
            end
            
            temp_A = [A, col];
            if rank(temp_A) > rank_current
                rank_current = rank_current + 1;
                A = temp_A;
                RoutePaths{route_idx} = path;
                route_idx = route_idx + 1;
                break; 
            else
                for step = 1:(length(path)-1)
                    u = path(step); v = path(step+1);
                    idx = find(s == u & t == v, 1);
                    weights(idx) = weights(idx) + 0.5; 
                end
            end
        end
        
        if attempt == 20
             warning(['Linearly independent logical path not found for O_', num2str(i), ' -> D_', num2str(j)]);
        end
    end
end

fprintf('\n--- DEFINED DOMINANT ROUTES (n=25) ---\n');
used_nodes = [];

for k = 1:n
    path = RoutePaths{k};
    used_nodes = [used_nodes, path]; 
    path_str = '';
    
    for step = 1:length(path)
        node = path(step);
        
        if node >= 1 && node <= 5
            node_name = sprintf('O_%d', node);
        elseif node >= 6 && node <= 10
            node_name = sprintf('D_%d', node - 5);
        elseif node >= 11 && node <= 25
            node_name = sprintf('J_%d', node - 10); 
        end
        
        if step == 1
            path_str = node_name;
        else
            path_str = [path_str, ' -> ', node_name];
        end
    end
    
    fprintf('Route %2d: %s\n', k, path_str);
end

unique_used = unique(used_nodes);
internal_junctions = 11:25;
unused_junctions = setdiff(internal_junctions, unique_used);

fprintf('\n--- JUNCTION UTILIZATION REPORT ---\n');
if isempty(unused_junctions)
    fprintf('SUCCESS: 100%% of Junctions (J_1 to J_15) are actively used!\n');
else
    fprintf('WARNING: The following junctions are unused: ');
    for u = 1:length(unused_junctions)
        fprintf('J_%d ', unused_junctions(u) - 10);
    end
    fprintf('\n');
end
I_int = zeros(p_int, m);
J_all = [J_out, J_in];
for idx = 1:p_int
    vla = J_all(idx);
    for e = 1:m
        if s(e) == vla, I_int(idx, e) = -1; end 
        if t(e) == vla, I_int(idx, e) = 1; end  
    end
end

M = [A; I_int * A];
fprintf('\n--- STRUCTURAL OBSERVABILITY CHECK ---\n');
fprintf('Actual Rank of M  = %d (Expected: 25)\n', rank(M));


figure('Name', 'H-MASE Topology', 'Color', 'w', 'Position', [100, 100, 850, 850]);

nodeLabels = cell(1, 25);
for i = 1:25
    if i <= 5
        nodeLabels{i} = sprintf('O_%d', i);           
    elseif i <= 10
        nodeLabels{i} = sprintf('D_%d', i - 5);       
    else
        nodeLabels{i} = sprintf('J_%d', i - 10);      
    end
end


angles_in = linspace(pi/2, pi/2 + 2*pi, 6); angles_in(end) = []; 
angles_out = zeros(1, 10);
for i = 1:5
    angles_out(2*i-1) = angles_in(i) + 0.25; 
    angles_out(2*i) = angles_in(i) - 0.25;   
end

X = zeros(1, 25); Y = zeros(1, 25);


X(J_in) = 1.5 * cos(angles_in); Y(J_in) = 1.5 * sin(angles_in);


X(J_out) = 3.5 * cos(angles_out); Y(J_out) = 3.5 * sin(angles_out);

X(O) = 5.5 * cos(angles_out(1:2:9)); Y(O) = 5.5 * sin(angles_out(1:2:9));
X(D) = 5.5 * cos(angles_out(2:2:10)); Y(D) = 5.5 * sin(angles_out(2:2:10));

G_vis = digraph(s, t);
p_graph = plot(G_vis, 'XData', X, 'YData', Y, 'LineWidth', 1.2, ...
               'ArrowSize', 10, 'MarkerSize', 8, 'NodeLabel', nodeLabels);

highlight(p_graph, O, 'NodeColor', '#77AC30', 'Marker', 's'); 
highlight(p_graph, D, 'NodeColor', '#D95319', 'Marker', '^'); 
highlight(p_graph, J_out, 'NodeColor', '#4DBEEE', 'Marker', 'd'); 
highlight(p_graph, J_in, 'NodeColor', '#0072BD', 'Marker', 'o'); 


highlight(p_graph, J_in, J_in([2 3 4 5 1]), 'EdgeColor', '#EDB120', 'LineWidth', 3.5); 
for k=1:10
    highlight(p_graph, [J_out(k), J_out(mod(k,10)+1)], [J_out(mod(k,10)+1), J_out(k)], 'EdgeColor', '#7E2F8E', 'LineWidth', 2);
end

title('H-MASE Urban Hub Model');
axis equal; axis off; set(gca, 'FontSize', 12);
exportgraphics(gcf, 'Fig_NetworkGraph.eps', 'ContentType', 'vector');
edge_usage_count = sum(A, 2); 
used_edge_indices = find(edge_usage_count > 0);
unused_edge_indices = find(edge_usage_count == 0);

fprintf('\n--- EDGE UTILIZATION REPORT ---\n');
fprintf('Total Number of Edges (m): %d\n', m);
fprintf('Edges Used by Routes:      %d\n', length(used_edge_indices));
fprintf('Unused (Idle) Edges:       %d\n', length(unused_edge_indices));

utilization_ratio = (length(used_edge_indices) / m) * 100;
fprintf('Network Edge Coverage:     %.2f%%\n', utilization_ratio);

if ~isempty(unused_edge_indices)
    fprintf('\n--- LIST OF UNUSED EDGES ---\n');
    for i = 1:length(unused_edge_indices)
        idx = unused_edge_indices(i);
        u = s(idx); v = t(idx);
        
        labels = cell(1, 2);
        nodes = [u, v];
        for n_idx = 1:2
            curr = nodes(n_idx);
            if curr <= 5
                labels{n_idx} = sprintf('O_%d', curr);
            elseif curr <= 10
                labels{n_idx} = sprintf('D_%d', curr-5);
            else
                labels{n_idx} = sprintf('J_%d', curr-10);
            end
        end
        fprintf('Edge %2d: %s -> %s (Not in any dominant route)\n', idx, labels{1}, labels{2});
    end
else
    fprintf('\nSUCCESS: 100%% of physical links are covered by the dominant route set!\n');
end

[max_val, max_idx] = max(edge_usage_count);
u_max = s(max_idx); v_max = t(max_idx);
fprintf('\n--- TRAFFIC HOTSPOT ANALYSIS ---\n');
fprintf('Most heavily used link: Edge %d (Used by %d routes)\n', max_idx, max_val);



N = m + p_int; 

A_adj = zeros(N, N);
J_all = [J_out, J_in];

for i = 1:m
    for j = i+1:m
        if s(i) == s(j) || s(i) == t(j) || t(i) == s(j) || t(i) == t(j)
            A_adj(i, j) = 1;
            A_adj(j, i) = 1;
        end
    end
end

for k = 1:p_int
    vla_node = J_all(k); 
    vla_agent_id = m + k;
    for e = 1:m
        if s(e) == vla_node || t(e) == vla_node
            A_adj(vla_agent_id, e) = 1;
            A_adj(e, vla_agent_id) = 1;
        end
    end
end

for k1 = 1:p_int
    for k2 = k1+1:p_int
        node1 = J_all(k1);
        node2 = J_all(k2);
        for e = 1:m
            if (s(e) == node1 && t(e) == node2) || (s(e) == node2 && t(e) == node1)
                A_adj(m + k1, m + k2) = 1;
                A_adj(m + k2, m + k1) = 1;
                break;
            end
        end
    end
end

deg = sum(A_adj, 2);
W = zeros(N, N);
for i = 1:N
    for j = 1:N
        if i ~= j && A_adj(i,j) == 1
            W(i,j) = 1 / (max(deg(i), deg(j)) + 1);
        end
    end
    W(i,i) = 1 - sum(W(i,:)); 
end

T_max = 150;     
K = 80;
T_init = 40;


omega_max = 2.0;
delta_max = 1.0; 
Gamma = omega_max * 2.5;

M_dagger = cell(N, 1);
P = cell(N, 1);

for i = 1:N
    M_i = M(i, :);
    
    if norm(M_i) > 0
        M_dagger{i} = M_i' / (norm(M_i)^2);
    else
        M_dagger{i} = zeros(n, 1);
    end
    P{i} = eye(n) - M_dagger{i} * M_i;     
end

r_hat = ones(n, N) * 20; 

r_true = zeros(n, T_max);

r_true(:, 1) = 20 + 20 * rand(n, 1); 


tilde_omega = -omega_max + 2*omega_max * rand(m, 1); 
y_p_init = A * r_true(:, 1) + tilde_omega; 
b_init = [y_p_init; zeros(p_int, 1)];

r_hat = zeros(n, N);
for i = 1:N
    r_hat(:, i) = M_dagger{i} * b_init(i); 
end

t_fault = 80;
faulty_agent = 1;
f_magnitude = 80.0;

global_error = zeros(T_max, 1);
eta_log = zeros(T_max, N);
sigma_log = ones(T_max, N);

omega_max = 2.0; 
delta_max = 1.0; 

Gamma_steady = omega_max * 2.5;
Gamma_init = 200;              
decay_rate = 0.15;             
Gamma_t = zeros(T_max, 1);
for t = 1:T_max
    Gamma_t(t) = Gamma_steady + Gamma_init * exp(-decay_rate * t);
end

r_hat_avg_log = zeros(n, T_max);

for t = 1:T_max
    
    if t > 1
        Delta_r = -delta_max + 2*delta_max * rand(n, 1);
        r_true(:, t) = max(0, r_true(:, t-1) + Delta_r); 
    end
    
    tilde_omega = -omega_max + 2*omega_max * rand(m, 1); 
    
    y_p = A * r_true(:, t) + tilde_omega; 
    
    if t >= t_fault
        y_p(faulty_agent) = y_p(faulty_agent) + f_magnitude;
    end
    
    b = [y_p; zeros(p_int, 1)]; 
    
    for k = 1:K
        r_hat_next = zeros(n, N);
        
        for i = 1:N
            r_diff_i = zeros(n, 1);
            for j = 1:N
                if W(i,j) > 0
                    r_diff_i = r_diff_i + W(i,j) * r_hat(:, j);
                end
            end
            
            eta_i = norm(M(i,:) * r_diff_i - b(i));
            eta_log(t, i) = eta_i; 
            
            if eta_i > Gamma_t(t)
                sigma_i = 0;
            else
                sigma_i = 1;
            end
            sigma_log(t, i) = sigma_i;
            
           
            if sigma_i == 1
                r_hat_next(:, i) = P{i} * r_diff_i + M_dagger{i} * b(i);
            else
                r_hat_next(:, i) = r_diff_i; 
            end
        end
        r_hat = r_hat_next;
    end
    
    E_t_K = 0;
    for i = 1:N
        e_i = r_hat(:, i) - r_true(:, t);
        E_t_K = E_t_K + norm(e_i)^2;
    end
    global_error(t) = sqrt(E_t_K);
    r_hat_avg_log(:, t) = mean(r_hat, 2);
end


figure('Color', 'w', 'Position', [100, 100, 800, 600]);

subplot(3,1,1);
plot(1:T_max, global_error, 'k-', 'LineWidth', 1.5); hold on;
xline(t_fault, 'r--', 'LineWidth', 1.5);
ylabel('$\|e(t,K)\|$', 'Interpreter', 'latex', 'FontSize', 12);
grid on; set(gca, 'TickLabelInterpreter', 'latex');

subplot(3,1,2);
plot(1:T_max, eta_log(:, faulty_agent), 'b-', 'LineWidth', 1.5); hold on;
plot(1:T_max, Gamma_t, 'r--', 'LineWidth', 1.5); 
ylabel(['$\eta_{' num2str(faulty_agent) '}(t)$'], 'Interpreter', 'latex', 'FontSize', 12);
grid on; set(gca, 'TickLabelInterpreter', 'latex');
legend({'Residual', '$\Gamma(t)$ Threshold'}, 'Interpreter', 'latex', 'Location', 'NorthEast');

subplot(3,1,3);
stairs(1:T_max, sigma_log(:, faulty_agent), 'k-', 'LineWidth', 1.5);
ylim([-0.2 1.2]); yticks([0 1]);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(['$\sigma_{' num2str(faulty_agent) '}(t)$'], 'Interpreter', 'latex', 'FontSize', 12);
grid on; set(gca, 'TickLabelInterpreter', 'latex');


fprintf('\n=======================================================\n');
fprintf('   H-MASE THEORETICAL VALIDATION & COMPLIANCE REPORT\n');
fprintf('=======================================================\n');

pre_fault_errors = global_error(T_init:t_fault-1);
mean_pre_error = mean(pre_fault_errors);
std_pre_error = std(pre_fault_errors);

fprintf('\n[Theorem 1 - ISS Convergence]\n');
fprintf('Expected: Error bounded within an invariant set (> 0 due to noise).\n');
fprintf('Result:   Mean Error = %.2f, Std Dev = %.2f\n', mean_pre_error, std_pre_error);
if std_pre_error < mean_pre_error * 0.5 && mean_pre_error < 100
    fprintf('Status:   PASS (Error is bounded and stable)\n');
else
    fprintf('Status:   FAIL (Error is diverging or highly unstable)\n');
end

nominal_residuals = eta_log(T_init:t_fault-1, faulty_agent);
max_nominal_res = max(nominal_residuals);

fprintf('\n[Lemma 3 - Nominal Consistency]\n');
fprintf('Expected: Residual \\eta_i(t) < \\Gamma (%.2f) before fault.\n', Gamma);
fprintf('Result:   Maximum Nominal Residual = %.2f\n', max_nominal_res);
if max_nominal_res < Gamma
    fprintf('Status:   PASS (No false alarms generated)\n');
else
    fprintf('Status:   FAIL (False alarms triggered, \\Gamma is too low)\n');
end

fault_residuals = eta_log(t_fault:T_max, faulty_agent);
min_fault_res = min(fault_residuals);
isolation_status = sum(sigma_log(t_fault:T_max, faulty_agent)) == 0; 

fprintf('\n[Fault Isolation Logic - Eq. 20]\n');
fprintf('Expected: Residual \\eta_i(t) > \\Gamma and \\sigma_i = 0 after fault.\n');
fprintf('Result:   Minimum Fault Residual = %.2f\n', min_fault_res);
if min_fault_res > Gamma && isolation_status
    fprintf('Status:   PASS (Fault successfully detected and agent isolated)\n');
else
    fprintf('Status:   FAIL (Fault missed or isolation logic failed)\n');
end

post_fault_errors = global_error(t_fault+10:T_max);
mean_post_error = mean(post_fault_errors);
resilience_ratio = mean_post_error / mean_pre_error;

fprintf('\n[Theorem 2 - Resilience Check]\n');
fprintf('Expected: Global error remains bounded despite severe fault.\n');
fprintf('Result:   Post-fault Mean Error = %.2f (%.2fx of pre-fault)\n', mean_post_error, resilience_ratio);
if resilience_ratio < 2.0 
    fprintf('Status:   PASS (System maintained strong observability and resilience)\n');
else
    fprintf('Status:   FAIL (Fault corrupted the global estimation)\n');
end
fprintf('=======================================================\n');

steady_idx = 40:(t_fault-1);
rmse_val = 0;
mape_val = 0;
for t_idx = steady_idx
    rmse_val = rmse_val + (global_error(t_idx)^2 / (n * N));
    
    mean_true = mean(r_true(:, t_idx)) + 1e-4;
    mape_val = mape_val + ( (global_error(t_idx)/sqrt(n * N)) / mean_true );
end
rmse_steady = sqrt(rmse_val / length(steady_idx));
mape_steady = (mape_val / length(steady_idx)) * 100;

fprintf('\n[Estimation Performance Metrics]\n');
fprintf('Steady-State RMSE : %.2f vehicles\n', rmse_steady);
fprintf('Steady-State MAPE : %.2f%%\n', mape_steady);
fprintf('=======================================================\n');


figure('Name', 'Global Estimation Error', 'Color', 'w', 'Position', [150, 150, 800, 400]);

p1 = plot(1:T_max, global_error, 'b-', 'LineWidth', 1.5); hold on;

p2 = xline(t_fault, 'k--', 'LineWidth', 1.5);

xlabel('Time Step ($t$)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Global Estimation Error $\|e(t,K)\|$', 'Interpreter', 'latex', 'FontSize', 12);

grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
xlim([1, T_max]);
ylim([0, max(global_error) * 1.2]);


legend([p1, p2], ...
    {'H-MASE Global Error $\|e(t,K)\|$', ...
     'Fault Injection ($t_{fault}$)'}, ...
    'Interpreter', 'latex', 'Location', 'best', 'FontSize', 11);
exportgraphics(gcf, 'Fig_Global_Error.eps', 'ContentType', 'vector');

figure('Name', 'Route Densities Tracking', 'Color', 'w', 'Position', [100, 100, 900, 500]);
hold on;

colors = lines(n); 

for i = 1:n
    plot(1:T_max, r_true(i, :), '-.', 'Color', colors(i,:), 'LineWidth', 1.0);
    
    plot(1:T_max, r_hat_avg_log(i, :), '-', 'Color', colors(i,:), 'LineWidth', 1.5);
end

p_fault = xline(t_fault, 'k--', 'LineWidth', 2);

p_true = plot(NaN, NaN, 'k-.', 'LineWidth', 1.5);
p_est  = plot(NaN, NaN, 'k-',  'LineWidth', 1.5);

xlabel('Time Step ($t$)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Route Traffic Flows (vehicles/step)', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
xlim([1, T_max]);
ylim([0, max(max(r_true)) * 1.2]); 

legend([p_true, p_est, p_fault], ...
    {'True Route Densities ($r$)', ...
     'Average Estimated Densities ($\hat{r}$)', ...
     'Fault Injection ($t_{fault}$)'}, ...
    'Interpreter', 'latex', 'Location', 'NorthWest', 'FontSize', 11);
exportgraphics(gcf, 'Fig_Tracking.eps', 'ContentType', 'vector');
