% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2019/2020 Implemented by 
% Dr. Laurita dos Santos
% Dr. Debora C. Correa
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% based on: phylosophical transactions A 375 2017: Multiscale ordinal network
% analysis of human cardiac dynamics
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this code is part of paper intituled "Characterisation of neonatal cardiac dynamics using ordinal partition network"
% submitted for Medical & Biological Engineering & Computing - Springer


function [permutation_entropy, CP_entropy, gn_entropy] = ordinal_partition_network(file,m,tau)

N = length(file);

% ******************
%file = flip(file); %for backwards
% ******************

% building Z vectors
N_emb = (N-(m-1)*tau);
z = zeros(N_emb,m);
for ii=1:m
    z(:,ii) = file([1:N_emb] + tau*(ii-1));
end

% symbols - according to order in the vector - s_d = matrix with symbols
[lin, col] = size(z);
s_d = zeros(lin,col); % vector will receive symbols
for j=1:lin
    [~,sorted] = sort(z(j,:));
    for t=1:m
        pt(sorted(t)) = t;
    end
    s_d(j,:) = pt;
end

% verify all permuations - according to m (embedding dimension)
[words,~, id_words] = unique(s_d, 'rows'); 
numb = zeros(1,N_emb);
% comparison between symbols and permutation
for j=1:N_emb
    window = s_d(j,:);
    indexes = ismember(words,window,'rows');
    numb(j) = find(indexes==1);
end 

% build adjacent matrix - HERE is still allows self-loops
nWords = size(words,1);
matriz_adj = zeros(nWords, nWords);
% scanning the time series; incrementing the conditional probabilities in
% the transition matrix
for n=1:length(numb)-1
    matriz_adj(id_words(n), id_words(n + 1)) = matriz_adj(id_words(n), id_words(n + 1)) + 1;
end

% % permutation entropy:  it is estimated by counting the relative
% % ocurrence of each symbol in the symbolic dynamics S - equation 3.1
par = unique(numb);
n_par = length(par);
pi = zeros(1,n_par);
for l=1:n_par
    count = sum(numb==par(l));
    pi(l) = count/N_emb;
end

permutation_entropy = -1*sum(pi .*log2(pi));

% finding the stochastic matrix P with elements: equation 3.2
sumOfRows = sum(matriz_adj,2); non_zero_rows = (sumOfRows > 0);
matriz_P(non_zero_rows,:) = bsxfun(@rdivide, matriz_adj(non_zero_rows,:), sumOfRows(non_zero_rows));

% estimation of the stationary distribution from matriz_adj is - equation
% 3.4
matriz_pi = sum(matriz_adj,2)/sum(sum(matriz_adj));

% conditional permutation entropy (h^(CPE)) - equation 3.5
pij = matriz_P.*log2(matriz_P);
pij(isnan(pij))=0;
P = -1 .* pi .* sum(pij);
CP_entropy = sum(P);

% % % local and global node out-link entropy - REBUILD adjacency matrix
% % % finding the stochastic matrix P with elements WITHOUT self-loops -
% % % equation 3.6
matriz_adj2 = matriz_adj - diag(diag(matriz_adj));

% % % finding the NEW stochastic matrix P with elements: equation 3.2 - without
% self-loops
sumOfRows = sum(matriz_adj2,2); non_zero_rows = (sumOfRows > 0);
matriz_PT(non_zero_rows,:) = bsxfun(@rdivide, matriz_adj2(non_zero_rows,:), sumOfRows(non_zero_rows));

% % shannon entropy of row i gives the local node out-link entropy for node
% % i: equation 3.7 (shannon_pT) and 3.4 (matriz_pi2 - without self-loogps)
shannon = matriz_PT .* log2(matriz_PT);
shannon(isnan(shannon)) = 0;
shannon_pT = -1 * sum(shannon,2);
matriz_pi2 = sum(matriz_PT,2)/sum(sum(matriz_PT));

% % global node out-link entropy (GNE) - equation 3.8
gn_entropy = sum(shannon_pT .* matriz_pi2);



