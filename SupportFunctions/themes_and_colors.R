library(scales)
set_theme <- theme_bw()+theme(axis.title = element_text(size = 14))
color_palette_discrete <-  color_palette_discrete <- ggsci::pal_nejm("default")(8)
color_palette_discrete[1:2]  <- c("#E69F00", "#009E73")
color_palette_blue_white_red <- scale_colour_gradient2()
color_palette_blue_white <- scale_colour_gradient(high = muted("blue"), low = "white")

male <- "#20854EFF"
female <- "#E18727FF"

blue_color <- muted("blue")