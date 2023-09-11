
##############################
## Analyzing the running output
##############################

#####
# For batches
#####
library(ggplot2)
library(tidyr)
library(gridExtra)

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


data <- data.frame(nseq = nseq, Higher_order = lug1[,1], First_order = lug1[,2])
data_melted <- gather(data, key = "Type", value = "Batch_Size", -nseq)

g1 <- ggplot(data_melted, aes(x = nseq, y = Batch_Size, color = Type, linetype = Type)) +
  geom_line(size = 1) +
  geom_point(aes(shape = Type)) + # Optional: Adds point markers
  scale_color_manual(values = c("Higher_order" = "blue", "First_order" = "red")) +
  scale_linetype_manual(values = c("Higher_order" = "dashed", "First_order" = "solid")) +
  scale_shape_manual(values = c("Higher_order" = 16, "First_order" = 17)) + # Optional: Sets point markers
  scale_x_continuous(trans = 'log10') + # Log-scale for x-axis if needed
  scale_y_continuous(trans = 'log10') + # Log-scale for y-axis if needed
  labs(
    title = "Batch Size vs Sample Size",
    x = "Sample Size",
    y = "Batch Size",
    caption = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "top left",
    legend.title = element_blank(),
    plot.title = element_text(size = 8) # Set title size
  )
  



data2 <- data.frame(nseq = nseq, Higher_order = lug2[,1], First_order = lug2[,2])
data2_melted <- gather(data2, key = "Type", value = "Batch_Size", -nseq)

g2 <- ggplot(data2_melted, aes(x = nseq, y = Batch_Size, color = Type, linetype = Type)) +
  geom_line(size = 1) +
  geom_point(aes(shape = Type)) + # Optional: Adds point markers
  scale_color_manual(values = c("Higher_order" = "blue", "First_order" = "red")) +
  scale_linetype_manual(values = c("Higher_order" = "dashed", "First_order" = "solid")) +
  scale_shape_manual(values = c("Higher_order" = 16, "First_order" = 17)) + # Optional: Sets point markers
  scale_x_continuous(trans = 'log10') + # Log-scale for x-axis if needed
  scale_y_continuous(trans = 'log10') + # Log-scale for y-axis if needed
  labs(
    title = "Batch Size vs Sample Size",
    x = "Sample Size",
    y = "Batch Size",
    caption = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "top left",
    legend.title = element_blank(),
    plot.title = element_text(size = 8) # Set title size
  )




# Assuming nseq and lug3 are already loaded into your environment
data3 <- data.frame(nseq = nseq, Higher_order = lug3[,1], First_order = lug3[,2])
data3_melted <- gather(data3, key = "Type", value = "Batch_Size", -nseq)

g3 <- ggplot(data3_melted, aes(x = nseq, y = Batch_Size, color = Type, linetype = Type)) +
  geom_line(size = 1) +
  geom_point(aes(shape = Type)) +
  scale_color_manual(values = c("Higher_order" = "blue", "First_order" = "red")) +
  scale_linetype_manual(values = c("Higher_order" = "dashed", "First_order" = "solid")) +
  scale_shape_manual(values = c("Higher_order" = 16, "First_order" = 17)) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  labs(
    title = "Batch Size vs Sample Size (r = 3)",
    x = "Sample Size",
    y = "Batch Size",
    caption = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "top left",
    legend.title = element_blank(),
    plot.title = element_text(size = 8) # Set title size
  )

# Save the plot to a PDF file

ggsave("plots/sto_batch_r1_new.pdf", plot = g1, height = 6, width = 6)
ggsave("plots/sto_batch_r2_new.pdf", plot = g2, height = 6, width = 6)
ggsave("plots/sto_batch_r3_new.pdf", plot = g3, height = 6, width = 6)



gfinal <- grid.arrange(g1, g2, g3, ncol = 3,
             top = "Your Title Here",
             bottom = "Your Footer Here")
ggsave("plots/sto_batch_r_all.pdf", plot = gfinal, height = 3, width = 8)



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
data1 <- data.frame(nseq = nseq, Higher_order = lug1[,1], First_order = lug1[,2])
data1_melted <- gather(data1, key = "Type", value = "Determinant", -nseq)

gg1 <- ggplot(data1_melted, aes(x = nseq, y = Determinant, color = Type, linetype = Type)) +
  geom_line(size = 1) +
  geom_point(aes(shape = Type)) +
  scale_color_manual(values = c("Higher_order" = "blue", "First_order" = "red")) +
  scale_linetype_manual(values = c("Higher_order" = "dashed", "First_order" = "solid")) +
  scale_shape_manual(values = c("Higher_order" = 16, "First_order" = 17)) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  labs(
    x = "Sample Size",
    y = "Determinant",
    title = paste0("Determinant vs Sample Size (r = ", 1, ")")
  ) +
  theme_minimal() +
  theme(
    legend.position = "top left",
    legend.title = element_blank(),
    plot.title = element_text(size = 8) # Set title size
  )



data2 <- data.frame(nseq = nseq, Higher_order = lug2[,1], First_order = lug2[,2])
data2_melted <- gather(data2, key = "Type", value = "Determinant", -nseq)

gg2 <- ggplot(data2_melted, aes(x = nseq, y = Determinant, color = Type, linetype = Type)) +
  geom_line(size = 1) +
  geom_point(aes(shape = Type)) +
  scale_color_manual(values = c("Higher_order" = "blue", "First_order" = "red")) +
  scale_linetype_manual(values = c("Higher_order" = "dashed", "First_order" = "solid")) +
  scale_shape_manual(values = c("Higher_order" = 16, "First_order" = 17)) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  labs(
    x = "Sample Size",
    y = "Determinant",
    title = paste0("Determinant vs Sample Size (r = ", 2, ")")
  ) +
  theme_minimal() +
  theme(
    legend.position = "top left",
    legend.title = element_blank(),
    plot.title = element_text(size = 8) # Set title size
  )




data3 <- data.frame(nseq = nseq, Higher_order = lug3[,1], First_order = lug3[,2])
data3_melted <- gather(data3, key = "Type", value = "Determinant", -nseq)

gg3 <- ggplot(data3_melted, aes(x = nseq, y = Determinant, color = Type, linetype = Type)) +
  geom_line(size = 1) +
  geom_point(aes(shape = Type)) +
  scale_color_manual(values = c("Higher_order" = "blue", "First_order" = "red")) +
  scale_linetype_manual(values = c("Higher_order" = "dashed", "First_order" = "solid")) +
  scale_shape_manual(values = c("Higher_order" = 16, "First_order" = 17)) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  labs(
    x = "Sample Size",
    y = "Determinant",
    title = paste0("Determinant vs Sample Size (r = ", 3, ")")
  ) +
  theme_minimal() +
  theme(
    legend.position = "top left",
    legend.title = element_blank(),
    plot.title = element_text(size = 8) # Set title size
  )

ggsave(paste0("plots/sto_det_r", 1, ".pdf"), plot = gg1, height = 5, width = 5)
ggsave(paste0("plots/sto_det_r", 2, ".pdf"), plot = gg2, height = 5, width = 5)
ggsave(paste0("plots/sto_det_r", 3, ".pdf"), plot = gg3, height = 5, width = 5)


ggfinal <- grid.arrange(gg1, gg2, gg3, ncol = 3,
                       top = "Your Title Here",
                       bottom = "Your Footer Here")
ggsave("plots/sto_det_r_all.pdf", plot = ggfinal, height = 3, width = 8)
