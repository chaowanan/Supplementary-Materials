# SEMinR analysis for "When the Avatar Is First Worn" (informatics-4590796)
# Requires: install.packages("seminr") (version 2.6.0). Input: analysis_data.csv (I/N/F/P renamed pI/pN/pF/pP; CONG = participant and avatar gender match).
# Exp codes: A = not sure, B = never, C = used occasionally, D = use regularly. Dur codes: A = <1 h, B = 1-5 h, C = 6-10 h, D = >10 h per week.
suppressMessages(library(seminr)); set.seed(2026)
d <- read.csv("analysis_data.csv"); NB <- 5000; CORES <- 4
out <- list()
mm1 <- constructs(
  composite("I", single_item("pI")), composite("N", single_item("pN")),
  composite("F", single_item("pF")), composite("P", single_item("pP")),
  composite("PROP", multi_items("PROP",1:3), weights=mode_A),
  composite("AGEN", multi_items("AGEN",1:5), weights=mode_A),
  composite("CHNG", multi_items("CHNG",1:3), weights=mode_A),
  composite("PRIV", single_item("PRIV1")))
sm1 <- relationships(paths(from=c("I","N","F","P"), to=c("PROP","AGEN","CHNG","PRIV")))
stage1 <- function(dd) estimate_pls(dd, mm1, sm1, inner_weights=path_factorial)
m1 <- stage1(d); s1 <- summary(m1)
write.csv(s1$loadings, "r_s1_loadings.csv"); write.csv(s1$reliability, "r_s1_rel.csv")
write.csv(s1$validity$htmt, "r_htmt.csv"); write.csv(s1$validity$cross_loadings, "r_cross.csv")
write.csv(cor(m1$construct_scores), "r_lvcor.csv"); write.csv(s1$paths, "r_s1_paths.csv"); write.csv(s1$fSquare, "r_s1_f2.csv")
capture.output(print(s1$vif_antecedents), file="r_vif.txt")
b1 <- bootstrap_model(m1, nboot=NB, cores=CORES, seed=11); write.csv(summary(b1)$bootstrapped_paths, "r_s1_boot.csv")
# PLSpredict
set.seed(7); pp <- predict_pls(m1, technique=predict_DA, noFolds=10, reps=10)
res <- pp$items$PLS_out_of_sample_residuals; act <- pp$items$item_actuals[, colnames(res)]
q2i <- 1 - colSums(res^2)/colSums(sweep(act,2,colMeans(act))^2)
lmr <- pp$items$lm_out_of_sample_residuals
itm <- data.frame(item=colnames(res), Q2predict=q2i, RMSE_PLS=sqrt(colMeans(res^2)), RMSE_LM=sqrt(colMeans(lmr^2)))
write.csv(itm, "r_plspredict_items.csv", row.names=FALSE)
cs <- pp$composites; cons <- c("PROP","AGEN","CHNG","PRIV")
q2c <- sapply(cons, function(k) 1 - sum((cs$actuals_star[,k]-cs$composite_out_of_sample[,k])^2)/sum((cs$actuals_star[,k]-mean(cs$actuals_star[,k]))^2))
write.csv(data.frame(construct=cons, Q2predict=q2c), "r_plspredict_constructs.csv", row.names=FALSE)
# Stage 2
st2 <- function(sc, dims, covs=c(), seed=21) {
  dd <- as.data.frame(sc[, c("I","N","F","P",dims,covs)])
  lst <- list(composite("I",single_item("I")),composite("N",single_item("N")),composite("F",single_item("F")),composite("P",single_item("P")),
              composite("UP", dims, weights=mode_A))
  for (cv in covs) lst[[length(lst)+1]] <- composite(cv, single_item(cv))
  mm2 <- do.call(constructs, lst)
  m <- estimate_pls(dd, mm2, relationships(paths(from=c("I","N","F","P",covs), to="UP")), inner_weights=path_factorial)
  b <- bootstrap_model(m, nboot=NB, cores=CORES, seed=seed); s <- summary(m); sb <- summary(b)
  list(paths=s$paths, rel=s$reliability["UP",], load=s$loadings[dims,"UP"], w=s$weights[dims,"UP"], f2=s$fSquare[,"UP"],
       boot=sb$bootstrapped_paths, bload=sb$bootstrapped_loadings)
}
sc <- as.data.frame(m1$construct_scores)
dimsets <- list(M3=c("PROP","AGEN","CHNG"), M4=c("PROP","AGEN","CHNG","PRIV"), noCHNG=c("PROP","AGEN","PRIV"), PA=c("PROP","AGEN"))
for (nm in names(dimsets)) out[[nm]] <- st2(sc, dimsets[[nm]])
# UP-PRIV correlation (main model)
m3 <- estimate_pls(as.data.frame(sc), constructs(composite("I",single_item("I")),composite("N",single_item("N")),composite("F",single_item("F")),composite("P",single_item("P")),composite("UP",c("PROP","AGEN","CHNG"),weights=mode_A)),
      relationships(paths(from=c("I","N","F","P"),to="UP")), inner_weights=path_factorial)
out$rUPPRIV <- cor(m3$construct_scores[,"UP"], sc$PRIV)
# covariates n=113
d113 <- d[d$Sex %in% c("F","M"),]; rownames(d113) <- NULL
m113 <- stage1(d113); sc113 <- as.data.frame(m113$construct_scores)
sc113$CONG <- d113$CONG
# Exp: A = not sure, B = never, C = used occasionally (reference), D = use regularly
sc113$EXP_A <- as.numeric(d113$Exp=="A"); sc113$EXP_B <- as.numeric(d113$Exp=="B"); sc113$EXP_D <- as.numeric(d113$Exp=="D")
out$C0 <- st2(sc113, dimsets$M3, seed=31); out$C1 <- st2(sc113, dimsets$M3, c("CONG","EXP_A","EXP_B","EXP_D"), seed=41)
saveRDS(out, "r_out.rds")
for (nm in c("M3","M4","noCHNG","PA","C0","C1")) { cat("\n==",nm,"\n"); print(round(out[[nm]]$boot,3)); print(round(out[[nm]]$paths,3)); print(round(out[[nm]]$rel,3)); print(round(out[[nm]]$f2,3))}
cat("\nM3 loadings/weights\n"); print(round(out$M3$bload,3)); print(round(out$M3$w,3)); cat("r UP-PRIV", out$rUPPRIV, "\n")
