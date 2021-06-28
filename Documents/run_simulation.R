run_simulation<-function(beta_xy,nsnp,nid,rsq_gx,af,ux_eff,uy_eff){
  
  u <- rnorm(nid) # Generate a confounder
  g <- make_geno(nid=nid, nsnp=nsnp, af=af) # Generate genotypes with allele frequencies of 0.5
  effs <- choose_effects(nsnp=nsnp, totvar=rsq_gx) # These SNPs instrument some exposure, and together explain 5% of the variance
  x <- make_phen(effs=c(effs, ux_eff), indep=cbind(g, u)) # Create X - influenced by snps and the confounder
  y <- make_phen(effs=c(beta_xy, uy_eff), cbind(x, u)) # Create Y - negatively influenced by X and positively influenced by the confounder
  mr<-systemfit(y ~ x, method="2SLS", inst = ~ g)
  
  sim<-list(u=u,g=g,effs=effs,x=x,y=y,mr=mr)
  return(sim)
  
}
