
##############################
## Analyzing the running output
##############################

#####
# For batches
#####

load("stock_running")

nreps 	<- length(sims_for_stocks)
no.n <- length(sims_for_stocks[[1]])

lug1 <- matrix(0, nrow = no.n, ncol = 2)
lug2 <- matrix(0, nrow = no.n, ncol = 2)
lug3 <- matrix(0, nrow = no.n, ncol = 2)

for(i in 1:no.n)
{
	main_out <- sims_for_stocks

	lug1[i, ] <- colMeans(Reduce(rbind, lapply(main_out, function(x) x[[i]][[1]][1,])))
	lug2[i, ] <- colMeans(Reduce(rbind, lapply(main_out, function(x) x[[i]][[1]][2,])))
	lug3[i, ] <- colMeans(Reduce(rbind, lapply(main_out, function(x) x[[i]][[1]][3,])))
}

pdf("plots/sto_batch_r1.pdf", height = 6, width = 6)
plot(nseq, lug1[,1], type = "n", 
	ylim = range(lug1),
	ylab = "Batch Size", xlab = "Sample size")
for(i in 1:2) lines(nseq, lug1[,i], col = i+1)
legend("topleft", bty = "n", legend = c("Higher-order", "First-order"), lty = 1 , col = 2:3)
dev.off()

pdf("plots/sto_batch_r2.pdf", height = 6, width = 6)
plot(nseq, lug2[,1], type = "n", 
	ylim = range(lug2),
	ylab = "Batch Size", xlab = "Sample size")
for(i in 1:2) lines(nseq, lug2[,i], col = i+1)
legend("topleft", bty = "n", legend = c("Higher-order", "First-order"), lty = 1 , col = 2:3)
dev.off()

pdf("plots/sto_batch_r3.pdf", height = 6, width = 6)
plot(nseq, lug2[,1], type = "n", 
	ylim = range(lug3),
	ylab = "Batch Size", xlab = "Sample size", main = "r = 3")
for(i in 1:2) lines(nseq, lug3[,i], col = i+1)
legend("topleft", bty = "n", legend = c("Higher-order", "First-order"), lty = 1 , col = 2:3)
dev.off()



#####
# Running Determinant
#####


lug1 <- matrix(0, nrow = no.n, ncol = 3)
lug2 <- matrix(0, nrow = no.n, ncol = 3)
lug3 <- matrix(0, nrow = no.n, ncol = 3)

for(i in 1:no.n)
{
	main_out <- sims_for_stocks

	lug1[i, 1] <- mean(Reduce(c, lapply(main_out, function(x) det(x[[i]][[2]][[1]]) )))
	lug1[i, 2] <- mean(Reduce(c, lapply(main_out, function(x) det(x[[i]][[3]][[1]]) )))

	lug2[i, 1] <- mean(Reduce(c, lapply(main_out, function(x) det(x[[i]][[2]][[2]]) )))
	lug2[i, 2] <- mean(Reduce(c, lapply(main_out, function(x) det(x[[i]][[3]][[2]]) )))

	lug3[i, 1] <- mean(Reduce(c, lapply(main_out, function(x) det(x[[i]][[2]][[3]]) )))
	lug3[i, 2] <- mean(Reduce(c, lapply(main_out, function(x) det(x[[i]][[3]][[3]]) )))
}


lims <- range(c(lug1, lug2, lug3))
pdf(paste0("plots/sto_det_r", 1,".pdf"), height = 5, width = 5)
plot(nseq, lug1[, 1], type = "l", ylim = lims,
	ylab = "Determinant", xlab = "Sample Size", col = 2)
lines(nseq, lug1[, 2], col = 3)
legend("bottomright", bty = "n", legend = c("Higher-order", "First-order"),
 lty = 1, col = 2:3)
dev.off()

pdf(paste0("plots/sto_det_r", 2,".pdf"), height = 5, width = 5)
plot(nseq, lug2[, 1], type = "l", ylim = lims,
	ylab = "Determinant", xlab = "Sample Size",  col = 2)
lines(nseq, lug2[, 2], col = 3)
legend("topright", bty = "n", legend = c("Higher-order", "First-order"), 
	lty = 1 , col = 2:3)
dev.off()

pdf(paste0("plots/sto_det_r", 3,".pdf"), height = 5, width = 5)
plot(nseq, lug3[, 1], type = "l", ylim = lims,
	ylab = "Determinant", xlab = "Sample Size", col = 2)
lines(nseq, lug3[, 2], col = 3)
legend("bottomright", bty = "n", legend = c("Higher-order", "First-order"),
 lty = 1, col = 2:3)
dev.off()


