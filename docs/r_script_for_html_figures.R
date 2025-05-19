library(plotly)
library(htmlwidgets)
#library(htmltools)

load("cov/data/plotly_illustration_figures.RData")

saveWidget(figure31a, "docs/figure31.html", selfcontained = T) # for the html figura 3.1a and 3.1b are identical
saveWidget(figure32a, "docs/figure32a.html", selfcontained = T)
saveWidget(figure32b, "docs/figure32b.html", selfcontained = T)
saveWidget(figure37a, "docs/figure37a.html", selfcontained = T)
saveWidget(figure37b, "docs/figure37b.html", selfcontained = T)

load("cov/data/figure_3_9.RData")

saveWidget(figure39a, "docs/figure39a.html", selfcontained = T)
saveWidget(figure39b, "docs/figure39b.html", selfcontained = T)

load("cov/data/figure_3_14.RData")

saveWidget(figure_3_14a, "docs/figure314a.html", selfcontained = T)
saveWidget(figure_3_14b, "docs/figure314b.html", selfcontained = T)
