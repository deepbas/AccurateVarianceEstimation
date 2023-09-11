##############################
## Analyzing the output
##############################
load("ar1_batch")
library(ggplot2)
library(tidyr)

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

data1 <- data.frame(rho = rho, Exact = lug1[,1], Higher_order = lug1[,2], First_order = lug1[,3])
data1_melted <- gather(data1, key = "Method", value = "Batch_Size", -rho)


data2 <- data.frame(rho = rho, Exact = lug2[,1], Higher_order = lug2[,2], First_order = lug2[,3])
data2_melted <- gather(data2, key = "Method", value = "Batch_Size", -rho)


data3 <- data.frame(rho = rho, Exact = lug3[,1], Higher_order = lug3[,2], First_order = lug3[,3])
data3_melted <- gather(data3, key = "Method", value = "Batch_Size", -rho)


# 
# ggsave("plots/ar1_batchesr1_new.pdf", plot = gg1, height = 5, width = 5)
# ggsave("plots/ar1_batchesr2_new.pdf", plot = gg2, height = 5, width = 5)
# ggsave("plots/ar1_batchesr3_new.pdf", plot = gg3, height = 5, width = 5)


# Common plot theme
common_theme <- theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    plot.margin = margin(15, 15, 15, 20)
  )



# First plot
gg1 <- ggplot(data1_melted, aes(x = rho, y = Batch_Size, color = Method)) +
  geom_line(size = 1) +
  geom_point(aes(shape = Method)) +  # Add shapes
  labs(
    x = expression(rho),
    y = "b",
    subtitle = "r=1"
  ) + common_theme +
  ylim(c(0,750)) +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )

# The second plot
gg2 <- ggplot(data2_melted, aes(x = rho, y = Batch_Size, color = Method)) +
  geom_line(size = 1) +
  geom_point(aes(shape = Method)) +  # Add shapes
  labs(
    x = expression(rho),
    y = "b",
    subtitle = "r=2"
  ) + common_theme +
  ylim(c(0,750)) +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )

# The third plot
gg3 <- ggplot(data3_melted, aes(x = rho, y = Batch_Size, color = Method)) +
  geom_line(size = 1) +
  geom_point(aes(shape = Method)) +  # Add shapes
  labs(
    x = expression(rho),
    y = "b",
    subtitle = "r=3"
  ) + common_theme +
  ylim(c(0,750)) +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )


library(ggplot2)
library(gridExtra)
library(cowplot)
library(ggpubr)

# Modified function to return the combined grob
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = "bottom"))
  legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  
  main_plot <- do.call(arrangeGrob, c(plots, ncol = length(plots)))
  
  combined <- arrangeGrob(main_plot, legend, 
                          heights = unit.c(unit(1, "npc") - lheight, lheight),
                          ncol = 1)
  return(combined)
}

# Generate the combined grob
combined_plot <- grid_arrange_shared_legend(gg1, gg2, gg3)

# Convert to a ggplot object using cowplot's ggdraw and as.grob
combined_gg <- ggdraw() + draw_grob(combined_plot)

# Use ggsave to save the ggplot object
ggsave("plots/combined_batchesAR1.pdf", plot = combined_gg, width = 11, height = 5)
