
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

# Create melted data frames similar to the previous example
data1 <- data.frame(Sample_size = nseq, Exact = lug1[,1], Higher_order = lug1[,2], First_order = lug1[,3])
data1_melted <- gather(data1, key = "Method", value = "Batch_Size", -Sample_size)

data2 <- data.frame(Sample_size = nseq, Exact = lug2[,1], Higher_order = lug2[,2], First_order = lug2[,3])
data2_melted <- gather(data2, key = "Method", value = "Batch_Size", -Sample_size)

data3 <- data.frame(Sample_size = nseq, Exact = lug3[,1], Higher_order = lug3[,2], First_order = lug3[,3])
data3_melted <- gather(data3, key = "Method", value = "Batch_Size", -Sample_size)

# Shared plot theme
common_theme <- theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    plot.margin = margin(15, 15, 15, 20)
  )

# Plot for r=1
# Individual ggplot object for r=1
gg1 <- ggplot(data1_melted, aes(x = Sample_size, y = Batch_Size, color = Method)) +
  geom_line(size = 1) +
  geom_point(aes(shape = Method)) +
  labs(x = "n", y = "b", subtitle = "r=1") +
  common_theme + ylim(c(0,800)) +
  theme(axis.text.y = element_text(angle = 90, size = 8)) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(angle = 90, size = 8),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )

# Plot for r=2
gg2 <- ggplot(data2_melted, aes(x = Sample_size, y = Batch_Size, color = Method)) +
  geom_line(size = 1) +
  geom_point(aes(shape = Method)) +
  labs(
    x = "n",
    y = "b",
    subtitle = "r=2"
  ) + common_theme +  ylim(c(0,800)) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(angle = 90, size = 8),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )

# Plot for r=3
gg3 <- ggplot(data3_melted, aes(x = Sample_size, y = Batch_Size, color = Method)) +
  geom_line(size = 1) +
  geom_point(aes(shape = Method)) +
  labs(
    x = "n",
    y = "b",
    subtitle = "r=3"
  ) + common_theme +  ylim(c(0,800)) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(angle = 90, size = 8),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )

# Save each individual plot
ggsave("plots/ar1_run_r1_ggplot.pdf", plot = gg1, height = 5, width = 5)
ggsave("plots/ar1_run_r2_ggplot.pdf", plot = gg2, height = 5, width = 5)
ggsave("plots/ar1_run_r3_ggplot.pdf", plot = gg3, height = 5, width = 5)




# Function to handle shared legend across multiple plots
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
combined_plot <- grid_arrange_shared_legend(gg1,gg2,gg3)

# Convert to a ggplot object using cowplot's ggdraw and as.grob
combined_gg <- ggdraw() + draw_grob(combined_plot)

# Use ggsave to save the ggplot object
ggsave("plots/bVsn_rho99.pdf", plot = combined_gg, width = 11, height = 5)


