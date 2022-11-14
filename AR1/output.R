##############################
## Analysing the output
##############################
load("ar1_out")

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


	lug1[i, ] <- colMeans(Reduce(rbind, lapply(main_out, function(x) x[[1]][[1]][1,])))
	lug2[i, ] <- colMeans(Reduce(rbind, lapply(main_out, function(x) x[[1]][[1]][2,])))
	lug3[i, ] <- colMeans(Reduce(rbind, lapply(main_out, function(x) x[[1]][[1]][3,])))
}

pdf("plots/ar1_batches.pdf", height = 6, width = 6)
plot(rho, lug1[,1], type = "n", 
	ylim = range(cbind(lug1, lug2, lug3)),
	ylab = "Batch Size", xlab = expression(rho))
for(i in 1:3)
{
  lines(rho, lug1[,i], col = i)
  lines(rho, lug2[,i], col = i, lty = 2)
  lines(rho, lug3[,i], col = i, lty = 3)
}
legend("topleft", legend = c("Exact", "Higher-order", "First-order", 
                             "r = 1", "r = 2", "r = 3"), lty = c(1,1,1, 1, 2,3) , col = rep(1:3,2), ncol=2)
dev.off()

lug1 <- matrix(0, nrow = length(rho), ncol = 3)
lug2 <- matrix(0, nrow = length(rho), ncol = 3)
lug3 <- matrix(0, nrow = length(rho), ncol = 3)

for(i in 1:length(rho))
{
	Sigma <- true_Sigmas[[i]]
	phi <- phis[[i]]
	rho_this <- rho[i]

	main_out <- sims_for_rho[[i]]


	exact <- t(sapply(main_out, function(x) unlist(x[[2]])))
	second <- t(sapply(main_out, function(x) unlist(x[[3]])))
	first <- t(sapply(main_out, function(x) unlist(x[[4]])))
}

lug1 <- (cbind(exact[,1], second[,1], first[,1]))# - true_Sigmas[i][[1]][1,1])^2
lug2 <- (cbind(exact[,2], second[,2], first[,2]))# - true_Sigmas[i][[1]][1,1])^2
lug3 <- (cbind(exact[,3], second[,3], first[,3]))# - true_Sigmas[i][[1]][1,1])^2

colnames(lug1) <- c("exact", "second-order", "first-order")
colnames(lug2) <- c("exact", "second-order", "first-order")
colnames(lug3) <- c("exact", "second-order", "first-order")

par(mfrow = c(3,1))
boxplot(lug1, main = "r = 1")
abline(h = true_Sigmas[i][[1]][1,1])

boxplot(lug2, main = "r = 2")
abline(h = true_Sigmas[i][[1]][1,1])


boxplot(lug3, main = "r = 3")
abline(h = true_Sigmas[i][[1]][1,1])

