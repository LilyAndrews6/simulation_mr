summmary_level_data<-function(nsnp, af, rsq_gx, nid_x, b_xy, h2, nid_y){
#Code to create summary level data
#Generate effects that corresponds to the values
b_gx <- dplyr::tibble(af=af, snp=1:nsnp) %>%
  generate_gwas_params(h2=rsq_gx)
#use an LD correlation matrix to transform effects to reflect correlation structure 
bhat_gx <- generate_gwas_ss(b_gx, nid=nid_x)
#simulation of balanced pleiotropy
plei_index <- sample(1:nsnp, nsnp*0.5)
b_gy_plei <- dplyr::tibble(af=af[plei_index], snp=plei_index) %>%
  generate_gwas_params(h2)
#sum the genetic effects on y
b_gy <- b_gx
b_gy$beta <- b_gy$beta * b_xy
b_gy$beta[plei_index] <- b_gy$beta[plei_index] + b_gy_plei$beta
#generate summary data for g-y associations
bhat_gy <- generate_gwas_ss(b_gy, nid=nid_y)
}