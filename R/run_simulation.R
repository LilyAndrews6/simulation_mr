run_simulation<-function(beta_xy,nsnp,nid,rsq_gx,af,ux_eff,uy_eff){
  
  u <- rnorm(nid)
  g <- make_geno(nid=nid, nsnp=nsnp, af=af)
  effs <- choose_effects(nsnp=nsnp, totvar=rsq_gx)
  x <- make_phen(effs=c(effs, ux_eff), indep=cbind(g, u))
  y <- make_phen(effs=c(beta_xy, uy_eff), cbind(x, u))
  mr<-systemfit(y ~ x, method="2SLS", inst = ~ g)
  
  sim<-list(u=u,g=g,effs=effs,x=x,y=y,mr=mr)
  return(sim)
  
}
