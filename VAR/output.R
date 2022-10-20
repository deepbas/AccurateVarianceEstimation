##############################
## Analysing the output
##############################
load("var_out")

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

plot(rho, lug1[,1], type = "n", ylim = range(cbind(lug1, lug2, lug3)))
for(i in 1:3)
{
	lines(rho, lug1[,i], col = i)
	lines(rho, lug2[,i], col = i, lty = 2)
	lines(rho, lug3[,i], col = i, lty = 3)
}
legend("topleft", legend = c("Exact", "Higher-order", "First-order"), lty = 1, col = 1:3)