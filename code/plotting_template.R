# Thi Mui Pham, t.m.pham@hsph.harvard.edu
# ============================================================================ #
# GGplot templates for figures in manuscript
# ============================================================================ #
library(ggthemes)
theme_template <- function(base_size=14, base_family="helvetica"){
  (theme_foundation(base_size=base_size, base_family=base_family) + 
    theme_bw() + 
    theme(axis.line = element_line(colour="black"), 
          axis.title = element_text(size=18, face="bold"),
          axis.text = element_text(size=14),
          strip.text = element_text(size=22, face=c("bold.italic")),
          legend.title = element_text(size=20, face="bold"),
          legend.text = element_text(size=18),
          title = element_text(size=24, face="bold"),
          plot.background = element_rect(fill='transparent', color=NA), 
          legend.position = "none", 
          legend.background = element_rect(fill='transparent', color=NA), #transparent legend bg
          legend.box.background = element_rect(fill='transparent'), 
          panel.grid.minor = element_line(color="gray94", linetype = "solid"), 
          panel.grid.major = element_line(color="gray94", linetype = "solid"))
  )
}

theme_template_time_NI <- function(base_size=14, base_family="helvetica"){
  (theme_template()+
     theme(strip.text = element_text(size=22, face=c("bold")),
           axis.text.x = element_text(size=14, angle=0, hjust = 0.5))
   )
}

theme_template_time <- function(base_size=14, base_family="helvetica"){
  (theme_template()+
     theme(axis.text.x = element_text(size=14, angle=0, hjust = 0.5))
  )
}


theme_rmd_template <- function(base_size=14, base_family="helvetica"){
  (theme_foundation(base_size=base_size, base_family=base_family) + 
     theme_bw() + 
     theme(axis.line = element_line(colour="black"), 
           axis.title = element_text(size=24),
           axis.text = element_text(size=20),
           strip.text = element_text(size=18),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA), 
           panel.grid.minor = element_line(color="gray94", linetype = "solid"), 
           panel.grid.major = element_line(color="gray94", linetype = "solid"))
  )
}

theme_template_ioslides <- function(base_size=14, base_family="helvetica"){
  (theme_foundation(base_size=base_size, base_family=base_family) + 
     theme_bw() + 
     theme(axis.line = element_line(colour="black"), 
           axis.title = element_text(size=14),
           axis.text = element_text(size=10),
           strip.text = element_text(size=12),
           legend.title = element_text(size=12),
           legend.text = element_text(size=10),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA), 
           panel.grid.minor = element_line(color="gray94", linetype = "solid"), 
           panel.grid.major = element_line(color="gray94", linetype = "solid"))
  )
}
