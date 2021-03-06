clear
clc
close all

% Parameters
N=5;
f=.;
q=5;
iterations = 10000;
max_encoding_store = zeros(1,iterations);
W_mat_store = cell(1,iterations);
seq_mat_store = cell(1,iterations);

for i = 1:iterations
    W=zeros(N);
    % Generate sequence
    sequence = generate_seq(q, N, f);
    % Train weight matrix
    [W num_encoded] = train_matrix(W, sequence, q);
%     disp(['Able to encode ' num2str(num_encoded) ' patterns.'])
    
    max_encoding_store(i) = num_encoded;
    seq_mat_store{i} = sequence;
    W_mat_store{i} = W;
    
end

max(max_encoding_store)



function sequence = generate_seq(q, N, f)
    sequence = zeros(N,q);
    for i = 1:q
        pattern = rand(N,1)<f;
        sequence(:,i) = pattern;
    end
end

function pattern = generate_pattern(N, f)
    pattern = rand(N,1)<f;
end

function [out_mat num_encoded] = train_matrix(W, sequence, q)
    out_mat = W;
    for i = 2:q
        mask = sequence(:,i-1)*sequence(:,i)';
        out_mat = out_mat|mask;
        if seq_tester(out_mat, sequence(:,i-1), sequence(:,i))
            continue
        else
            break
        end
    end
    num_encoded = i-1;
end

function response = seq_tester(W, in_pattern, out_pattern)
    % Tests if W correctly encodes in pattern to out pattern.
    threshold = sum(in_pattern);
    if threshold == 0
        threshold = 1;
    end
    out_decoded = floor(in_pattern'*W/threshold);
    if out_decoded == out_pattern'
        response = true;
    else
        response = false;
    end
end






