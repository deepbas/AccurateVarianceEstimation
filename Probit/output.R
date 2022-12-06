##############################
## Analysing the output
##############################
load("var_out")

nreps 	<- length(sims_for_lupus)

lug1 <- matrix(0, nrow = 1, ncol = 3)
lug2 <- matrix(0, nrow = 1, ncol = 3)
lug3 <- matrix(0, nrow = 1, ncol = 3)

main_out <- sims_for_lupus

lug1 <- colMeans(Reduce(rbind, lapply(main_out, function(x) x[[1]][[1]][1,])))
lug2 <- colMeans(Reduce(rbind, lapply(main_out, function(x) x[[1]][[1]][2,])))
lug3 <- colMeans(Reduce(rbind, lapply(main_out, function(x) x[[1]][[1]][3,])))


##############################
## Analyzing the running output
##############################

#####
# For batches
#####

load("var1_running")

nreps 	<- length(sims_for_lupus)
no.n <- length(sims_for_lupus[[1]])

lug1 <- matrix(0, nrow = no.n, ncol = 2)
lug2 <- matrix(0, nrow = no.n, ncol = 2)
lug3 <- matrix(0, nrow = no.n, ncol = 2)
poli <- numeric(length = no.n)

# Sigma <- true_Sigmas
# phi <- phis


for(i in 1:no.n)
{
	main_out <- sims_for_lupus

	lug1[i, ] <- colMeans(Reduce(rbind, lapply(main_out, function(x) x[[i]][[1]][[1]][1,])))
	lug2[i, ] <- colMeans(Reduce(rbind, lapply(main_out, function(x) x[[i]][[1]][[1]][2,])))
	lug3[i, ] <- colMeans(Reduce(rbind, lapply(main_out, function(x) x[[i]][[1]][[1]][3,])))
	poli[i]		<- Reduce(mean, lapply(main_out, function(x) x[[i]][[1]][[2]]))
}

pdf("plots/var1_run_r1.pdf", height = 6, width = 6)
plot(nseq, lug1[,1], type = "n", 
	ylim = range(lug1),
	ylab = "Batch Size", xlab = "Sample size", main = "r = 1")
for(i in 1:2) lines(nseq, lug1[,i], col = i)
legend("topleft", bty = "n", legend = c("Higher-order", "First-order"), lty = 1 , col = 1:2)
dev.off()

pdf("plots/var1_run_r2.pdf", height = 6, width = 6)
plot(nseq, lug2[,1], type = "n", 
	ylim = range(lug2, poli),
	ylab = "Batch Size", xlab = "Sample size", main = "r = 2")
for(i in 1:2) lines(nseq, lug2[,i], col = i)
lines(nseq, poli, col = i, lty = 2)
legend("topleft", bty = "n", legend = c("Higher-order", "First-order", "Politis"), lty = c(1,1,1,2) , col = c(1:3, 3))
dev.off()

pdf("plots/var1_run_r3.pdf", height = 6, width = 6)
plot(nseq, lug2[,1], type = "n", 
	ylim = range(lug3),
	ylab = "Batch Size", xlab = "Sample size", main = "r = 3")
for(i in 1:2) lines(nseq, lug3[,i], col = i)
legend("topleft", bty = "n", legend = c("Higher-order", "First-order"), lty = 1 , col = 1:3)
dev.off()



#####
# Running Determinant
#####


lug1 <- matrix(0, nrow = no.n, ncol = 3)
lug2 <- matrix(0, nrow = no.n, ncol = 3)
lug3 <- matrix(0, nrow = no.n, ncol = 3)
poli <- numeric(length = no.n)

for(i in 1:no.n)
{
	main_out <- sims_for_lupus

	lug1[i, 1] <- mean(Reduce(c, lapply(main_out, function(x) det(x[[i]][[2]][[1]]) )))
	lug1[i, 2] <- mean(Reduce(c, lapply(main_out, function(x) det(x[[i]][[3]][[1]]) )))
	lug1[i, 3] <- mean(Reduce(c, lapply(main_out, function(x) det(x[[i]][[4]][[1]]) )))

	lug2[i, 1] <- mean(Reduce(c, lapply(main_out, function(x) det(x[[i]][[2]][[2]]) )))
	lug2[i, 2] <- mean(Reduce(c, lapply(main_out, function(x) det(x[[i]][[3]][[2]]) )))
	lug2[i, 3] <- mean(Reduce(c, lapply(main_out, function(x) det(x[[i]][[4]][[2]]) )))
	poli[i]		 <- mean(Reduce(c, lapply(main_out, function(x) det(x[[i]]$politBM) )))

	lug3[i, 1] <- mean(Reduce(c, lapply(main_out, function(x) det(x[[i]][[2]][[3]]) )))
	lug3[i, 2] <- mean(Reduce(c, lapply(main_out, function(x) det(x[[i]][[3]][[3]]) )))
	lug3[i, 3] <- mean(Reduce(c, lapply(main_out, function(x) det(x[[i]][[4]][[3]]) )))
}


lims <- range(c(lug1,  lug3, det(Sigma)))
pdf(paste0("plots/det_r", 1,".pdf"), height = 5, width = 5)
plot(nseq, lug1[, 1], type = "l", ylim = lims,
	ylab = "Determinant", xlab = "Sample Size")
lines(nseq, lug1[, 2], col = "red")
lines(nseq, lug1[, 3], col = "green")
abline(h = det(Sigma), col = "purple")
legend("bottomright", bty = "n", legend = c("Exact", "Higher-order", "First-order", "Truth"),
 lty = 1, col = c(1:3, "purple"))
dev.off()

pdf(paste0("plots/det_r", 2,".pdf"), height = 5, width = 5)
plot(nseq, lug2[, 1], type = "l", ylim = lims,
	ylab = "Determinant", xlab = "Sample Size")
lines(nseq, lug2[, 2], col = "red")
lines(nseq, lug2[, 3], col = "green")
lines(nseq, poli, col = "green", lty = 2)
legend("bottomright", bty = "n", legend = c("Exact", "Higher-order", "First-order", "Politis", "Truth"), 
	lty = c(1,1,1,2, 1) , col = c(1:3, 3, "purple"))
abline(h = det(Sigma), col = "purple")
dev.off()

pdf(paste0("plots/det_r", 3,".pdf"), height = 5, width = 5)
plot(nseq, lug3[, 1], type = "l", ylim = lims,
	ylab = "Determinant", xlab = "Sample Size")
lines(nseq, lug3[, 2], col = "red")
lines(nseq, lug3[, 3], col = "green")
legend("bottomright", bty = "n", legend = c("Exact", "Higher-order", "First-order", "Truth"),
 lty = 1, col = c(1:3, "purple"))
abline(h = det(Sigma), col = "purple")
dev.off()



#####
# Relative frobenius norm
#####


lug1 <- matrix(0, nrow = no.n, ncol = 3)
lug2 <- matrix(0, nrow = no.n, ncol = 3)
lug3 <- matrix(0, nrow = no.n, ncol = 3)
poli <- numeric(length = no.n)

Sigma <- as.matrix(Sigma)
for(i in 1:no.n)
{
	main_out <- sims_for_n

	lug1[i, 1] <- mean(Reduce(c, lapply(main_out, function(x) norm(x[[i]][[2]][[1]] - Sigma, "F")/norm(Sigma, "F") )))
	lug1[i, 2] <- mean(Reduce(c, lapply(main_out, function(x) norm(x[[i]][[3]][[1]] - Sigma, "F")/norm(Sigma, "F") )))
	lug1[i, 3] <- mean(Reduce(c, lapply(main_out, function(x) norm(x[[i]][[4]][[1]] - Sigma, "F")/norm(Sigma, "F") )))

	lug2[i, 1] <- mean(Reduce(c, lapply(main_out, function(x) norm(x[[i]][[2]][[2]] - Sigma, "F")/norm(Sigma, "F") )))
	lug2[i, 2] <- mean(Reduce(c, lapply(main_out, function(x) norm(x[[i]][[3]][[2]] - Sigma, "F")/norm(Sigma, "F") )))
	lug2[i, 3] <- mean(Reduce(c, lapply(main_out, function(x) norm(x[[i]][[4]][[2]] - Sigma, "F")/norm(Sigma, "F") )))
	poli[i]		 <- mean(Reduce(c, lapply(main_out, function(x) norm(x[[i]]$politBM - Sigma, "F")/norm(Sigma, "F") )))

	lug3[i, 1] <- mean(Reduce(c, lapply(main_out, function(x) norm(x[[i]][[2]][[3]] - Sigma, "F")/norm(Sigma, "F") )))
	lug3[i, 2] <- mean(Reduce(c, lapply(main_out, function(x) norm(x[[i]][[3]][[3]] - Sigma, "F")/norm(Sigma, "F") )))
	lug3[i, 3] <- mean(Reduce(c, lapply(main_out, function(x) norm(x[[i]][[4]][[3]] - Sigma, "F")/norm(Sigma, "F") )))
}


lims <- range(c(0, lug1,  lug3))
pdf(paste0("plots/norm_r", 1,".pdf"), height = 5, width = 5)
plot(nseq, lug1[, 1], type = "l", ylim = lims,
	ylab = "Relative Frobenius Norm", xlab = "Sample Size")
lines(nseq, lug1[, 2], col = "red")
lines(nseq, lug1[, 3], col = "green")
legend("top", bty = "n", legend = c("Exact", "Higher-order", "First-order"),
 lty = 1, col = 1:3)
dev.off()

pdf(paste0("plots/norm_r", 2,".pdf"), height = 5, width = 5)
plot(nseq, lug2[, 1], type = "l", ylim = lims,
	ylab = "Relative Frobenius Norm", xlab = "Sample Size")
lines(nseq, lug2[, 2], col = "red")
lines(nseq, lug2[, 3], col = "green")
lines(nseq, poli, col = "green", lty = 2)
legend("top", bty = "n", legend = c("Exact", "Higher-order", "First-order", "Politis"), lty = c(1,1,1,2) , col = c(1:3, 3))
dev.off()

pdf(paste0("plots/norm_r", 3,".pdf"), height = 5, width = 5)
plot(nseq, lug3[, 1], type = "l", ylim = lims,
	ylab = "Relative Frobenius Norm", xlab = "Sample Size")
lines(nseq, lug3[, 2], col = "red")
lines(nseq, lug3[, 3], col = "green")
legend("top", bty = "n", legend = c("Exact", "Higher-order", "First-order"),
 lty = 1, col = 1:3)
dev.off()


