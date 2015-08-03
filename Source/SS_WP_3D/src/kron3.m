function M = kron3(I,J,K)
  tmp = I.'*J;
  tmp = tmp(:)*K;
  M = reshape(tmp,[length(I),length(J),length(K)]);
  