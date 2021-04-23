

cross.entropy=function(gamma,true.states,quienes)
{
  out=0
  correctos=gamma[,quienes]
  for (i in 1:length(quienes))
  {
    out=out+log(correctos[true.states[quienes[i]],i])
  }
  return(-out)
}