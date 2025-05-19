library(plotly)
library(htmlwidgets)
library(htmltools)

load("cov/data/plotly_illustration_figures.RData")

saveWidget(tagList(h1("Figure 3.1"), figure31a), "docs/figure31.html", selfcontained = T) # for the html figura 3.1a and 3.1b are identical
saveWidget(tagList(h1("Figure 3.2 (a)"), figure32a), "docs/figure32a.html", selfcontained = T)
saveWidget(tagList(h1("Figure 3.2 (b)"), figure32b), "docs/figure32b.html", selfcontained = T)
saveWidget(tagList(h1("Figure 3.7 (a)"), figure37a), "docs/figure37a.html", selfcontained = T)
saveWidget(tagList(h1("Figure 3.7 (b)"), figure37b), "docs/figure37b.html", selfcontained = T)

