library(plotly)
library(htmlwidgets)

load("cov/data/plotly_illustration_figures.RData")

saveWidget(figure31a, "interactive_html_figures/figure31a.html", selfcontained = T)

