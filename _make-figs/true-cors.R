library(ggcorrplot)
library(ggplot2)
library(ggpubr)
library(cowplot)


t.cor = readRDS("../_data/_ground-truths//correlations.RDS")

plotCorr = function(corr, title){
  ggcorrplot(corr, type = "upper",
             outline.col = "black",
             colors = c("firebrick1", "gray95", "deepskyblue2"),
             lab = T, lab_size = 2.5) +
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text = element_blank(), 
          axis.title.x = element_blank(),  
          axis.title.y = element_blank(), 
          axis.ticks = element_blank(),
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank()) + scale_x_reverse()
}

plot_list = lapply(1:length(t.cor), function(i) plotCorr(t.cor[[i]], paste(i)))

legend_data = data.frame(cor_value = seq(0, 1, length.out = 100))
legend_data$y = 1

legend_plot = ggplot(legend_data, aes(x = cor_value, y = y, fill = cor_value)) +
  geom_tile() +
  scale_fill_gradient2(low = "gray95", 
                       high = "deepskyblue4", 
                       mid = "deepskyblue1", 
                       midpoint = 0.5) +
  theme_void() +
  theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(title = "Correlation", 
                                title.position = "top", 
                                barwidth = 10, 
                                barheight = 1))

legend_only = get_legend(legend_plot)
labels = c(letters[1:14], "legend")  

plot_list_with_legend = c(plot_list, list(legend_only))
labs = c(expression(A[1]), expression(A[2]), expression(A[3]), 
         expression(B[1]), expression(B[2]), expression(B[3]), 
         expression(C[1]), expression(C[2]), expression(C[3]), 
         expression(D), expression(E), "")

plot_to_save = do.call(plot_grid, 
                       list(plot_list_with_legend, 
                            ncol = 3, 
                            nrow = 5, 
                            labels = labs))


plot_to_save = plot_to_save + 
  theme(plot.background = element_rect(fill = "white", colour = "white"))
ggsave("../_manuscript/_figs/true_cors.png", 
       plot_to_save, 
       width = 6.6, 
       height = 8, 
       dpi = 300)
