two_sample_mr_simulation<-function(nsnp,nid,rsq_gx,rsq_gy,maf){
  # Genotypes for 20 SNPs and 10000 individuals, where the SNPs have MAF = 0.5:
  g <- make_geno(nid, nsnp, maf)
  # These SNPs instrument some exposure, and together explain 5% of the variance
  effs <- choose_effects(nsnp, rsq_gx)
  # Create X
  x <- make_phen(effs, g)
  # Check that the SNPs explain 5% of the variance in x
  sum(cor(x, g)^2)
  # Create Y, which is negatively influenced by x, explaining 10% of the variance in Y
  # so if the variances of x and y are both 1 then the effect size is about -0.31
  beta_xy <- -sqrt(rsq_gy)
  y <- make_phen(beta_xy, x)
  #perform 2 sample MR on summary data
  dat <- get_effs(x, y, g)
  TwoSampleMR::mr(dat)
}