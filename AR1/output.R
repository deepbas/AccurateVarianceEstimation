##############################
## Analysing the output
##############################
load("ar1_batch")

nreps 	<- length(sims_for_rho[[1]])

lug1 <- matrix(0, nrow = length(rho), ncol = 3)
lug2 <- matrix(0, nrow = length(rho), ncol = 3)
lug3 <- matrix(0, nrow = length(rho), ncol = 3)

for(i in 1:length(rho))
{
	Sigma <- true_Sigmas[[i]]
	phi <- phis[[i]]
	rho_this <- rho[i]

	main_out <- sims_for_rho[[i]]

	lug1[i, ] <- Reduce(c, lapply(main_out, function(x) x[[1]][1,]))
	lug2[i, ] <- Reduce(c, lapply(main_out, function(x) x[[1]][2,]))
	lug3[i, ] <- Reduce(c, lapply(main_out, function(x) x[[1]][3,]))
}

pdf("plots/ar1_batchesr1.pdf", height = 5, width = 5)
plot(rho, lug1[,1], type = "n", 
	ylim = range(lug1),
	ylab = "Batch Size", xlab = expression(rho))
for(i in 1:3) lines(rho, lug1[,i], col = i)
legend("topleft", legend = c("Exact", "Higher-order", "First-order"), 
	lty = 1 , col = 1:3, bty = "n")
dev.off()

pdf("plots/ar1_batchesr2.pdf", height = 5, width = 5)
plot(rho, lug2[,1], type = "n", 
	ylim = range(lug2),
	ylab = "Batch Size", xlab = expression(rho))
for(i in 1:3) lines(rho, lug2[,i], col = i)
legend("topleft", legend = c("Exact", "Higher-order", "First-order"),
 lty = 1 , col = 1:3, bty = "n")
dev.off()

pdf("plots/ar1_batchesr3.pdf", height = 5, width = 5)
plot(rho, lug2[,1], type = "n", 
	ylim = range(lug3),
	ylab = "Batch Size", xlab = expression(rho))
for(i in 1:3) lines(rho, lug3[,i], col = i)
legend("topleft", legend = c("Exact", "Higher-order", "First-order"),
 lty = 1 , col = 1:3, bty = "n")
dev.off()




rm(list = ls())


##############################
## Analysing the running output
##############################
load("ar1_running")

nreps 	<- length(sims_for_n)
no.n <- length(sims_for_n[[1]])

lug1 <- matrix(0, nrow = no.n, ncol = 3)
lug2 <- matrix(0, nrow = no.n, ncol = 3)
lug3 <- matrix(0, nrow = no.n, ncol = 3)

rho <- .99

Sigma <- true_Sigmas
phi <- phis
rho_this <- rho


for(i in 1:no.n)
{
	main_out <- sims_for_n

	lug1[i, ] <- colMeans(do.call(rbind, lapply(main_out, function(x) x[[i]][[1]][1,])))
	lug2[i, ] <- colMeans(do.call(rbind, lapply(main_out, function(x) x[[i]][[1]][2,])))
	lug3[i, ] <- colMeans(do.call(rbind, lapply(main_out, function(x) x[[i]][[1]][3,])))
}

pdf("plots/ar1_run_r1.pdf", height = 5, width = 5)
plot(nseq, lug1[,1], type = "n", 
	ylim = range(lug1),
	ylab = "Batch Size", xlab = "Sample size", main = "r = 1")
for(i in 1:3) lines(nseq, lug1[,i], col = i)
legend("topleft", legend = c("Exact", "Higher-order", "First-order"), 
	lty = 1 , col = 1:3, bty = "n")
dev.off()

pdf("plots/ar1_run_r2.pdf", height = 5, width = 5)
plot(nseq, lug2[,1], type = "n", 
	ylim = range(lug2),
	ylab = "Batch Size", xlab = "Sample size", main = "r = 2")
for(i in 1:3) lines(nseq, lug2[,i], col = i)
legend("topleft", legend = c("Exact", "Higher-order", "First-order"), 
	lty = 1 , col = 1:3, bty = "n")
dev.off()

pdf("plots/ar1_run_r3.pdf", height = 5, width = 5)
plot(nseq, lug2[,1], type = "n", 
	ylim = range(lug3),
	ylab = "Batch Size", xlab = "Sample size", main = "r = 3")
for(i in 1:3) lines(nseq, lug3[,i], col = i)
legend("topleft", legend = c("Exact", "Higher-order", "First-order"), 
	lty = 1 , col = 1:3, bty = "n")
dev.off()
