### Viterbi implementation for hsmm
library(hsmm)

# from hsmm library
Viterbi <- function(tau, J, M, VT.d, VT.p.tpm, VT.pi.ini, VT.pdf, VT.hiddenStates){
  .C("Viterbi", tau, J, M, VT.d, VT.p.tpm, VT.pi.ini, VT.pdf, VT.hiddenStates, PACKAGE="hsmm")
}  


ViterbiImpl=function (pi, pdf, d,tpm)
{
  # Inicializamos parametros, vectores y matrices
  J = length(pi)
  tau = dim(pdf)[2]
  M = max(dim(d)) ## Chequear esto
  
  vt.pdf=as.vector(t(pdf))
  vt.tpm=array(as.vector(tpm))
  vt.d=as.vector(t(d))
  hs= rep(as.integer(0), times = tau)
  results = Viterbi(tau, J, M,vt.d, vt.tpm, pi, vt.pdf, hs)[[8]]+1
  
  return (results)
}
