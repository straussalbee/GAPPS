#Set custom theme for all figures.
theme_dara <- theme_bw() +
  theme(
    title = element_text(size=7,colour="black",face="plain"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#ECECEC"),
    panel.border = element_rect(colour="gray80"),
    axis.text = element_text(size=6,colour="black"),
    axis.ticks = element_blank(),
    axis.ticks.length = unit(0, "mm"),
    axis.title = element_text(size=7,colour="black"),
    axis.title.x = element_text(vjust=0.5),
    axis.title.y = element_text(vjust=0.5),
    panel.margin = unit(0,"mm"), #note this can only have one parameter to be compatible with facet_grid
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size=7,colour="black",face="plain"),
    plot.margin = unit(c(1,1,1,1),"mm"),
    legend.text=element_text(size=6,colour="black",face="plain"),
    legend.title=element_text(size=7,colour="black",face="plain"),
    legend.background = element_blank(),
    legend.key.size = unit(1, "cm"),
    legend.justification=c(1,0)
  )

theme_set(theme_dara)
