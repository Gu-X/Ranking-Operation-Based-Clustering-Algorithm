%% Copyright (c) 2019, Xiaowei Gu

%% All rights reserved. Please read the "license.txt" for license terms.

%% This code is the Ranking Operation-based Clustering Algorithm described in:
%==========================================================================================================
%% X. Gu, P. Angelov, Z. Zhao, "A distance-type-insensitive clustering approach," 
%% Applied Soft Computing, vol. 77, pp. 622-634, 2019.
%==========================================================================================================
%% Please cite the paper above if this code helps.

%% For any queries about the code, please contact Dr. Xiaowei Gu
%% x.gu3@lancaster.ac.uk

%% Programmed by Xiaowei Gu
function [Centres,IDX]=ROC(data,a,M)
%%% Input
%% data    - data for clustering, each row represents a data sample
%% a       - control the sparsity of the sparse ranking matrices, the value range is (0,1).
%% M       - control the number of neighbours for local maxima identification, the value should be larger than 0.
%% The recommended values of a and M: a=0.9; M=10.
%%% Output
%% Centres - Centres of the clusters; each row is a centre.
%% IDX     - Cluster indices of data the samples; each row is the index of the data sample in the same row in the input data matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[L,W]=size(data);
S=zeros(L,L,W);
Threshold=round(a*L*(L-1)/2);
R=zeros(L,L,W);
for ii=1:1:W
    [~,seq1]=sort(pdist(data(:,ii),'minkowski',1),'descend');
    [~,seq2]=sort(seq1,'ascend');
    seq3=squareform(seq2);
    R(:,:,ii)=seq3;
    seq3(find(seq3<=Threshold))=0;
    S(:,:,ii)=seq3+eye(L);
end
R=sum(R,3);
S=sum(sum(S,3),1);
[~,seq4]=sort(R,'descend');
[~,seq5]=max([S;S(seq4(1:1:M,:))]);
centerID=find(seq5==1);
Centres=data(centerID,:);
[~,IDX]=max(R(:,centerID),[],2);
IDX(centerID)=1:1:length(centerID);
end
