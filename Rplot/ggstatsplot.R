library(ggstatsplot)
library(help="ggstatsplot")
library(ggplot2)
head(movies_long)
## # A tibble: 6 x 8
##   title                        year length budget rating  votes mpaa  genre
##   <chr>                       <int>  <int>  <dbl>  <dbl>  <int> <fct> <fct>
## 1 Shawshank Redemption, The    1994    142     25    9.1 149494 R     Drama
## 2 Lord of the Rings: The Ret~  2003    251     94    9   103631 PG-13 Acti~
## 3 Lord of the Rings: The Fel~  2001    208     93    8.8 157608 PG-13 Acti~
## 4 Lord of the Rings: The Two~  2002    223     94    8.8 114797 PG-13 Acti~
## 5 Pulp Fiction                 1994    168      8    8.8 132745 R     Drama
## 6 Schindler's List             1993    195     25    8.8  97667 R     Drama
ggbetweenstats(
  data = movies_long,
  x = mpaa, # > 2 groups
  y = rating,
  type = "p", # default
  messages = FALSE
)
ggbetweenstats(
  data = movies_long,
  x = mpaa,
  y = rating
)
ggbetweenstats(
  data = movies_long,
  x = mpaa,
  y = rating,
  type = "np",
  mean.ci = TRUE,
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  p.adjust.method = "fdr",
  messages = FALSE
)
ggbetweenstats(
  data = movies_long,
  x = mpaa,
  y = rating,
  type = "r",
  conf.level = 0.99,
  pairwise.comparisons = TRUE,
  pairwise.annotation = "p", 
  outlier.tagging = TRUE,
  outlier.label = title,
  outlier.coef = 2,
  ggtheme = hrbrthemes::theme_ipsum_tw(),
  palette = "Darjeeling2",
  package = "wesanderson",
  messages = FALSE
)
ggwithinstats(
  data = WRS2::WineTasting,
  x = Wine, # > 2 groups
  y = Taste,
  pairwise.comparisons = TRUE,
  pairwise.annotation = "p",
  ggtheme = hrbrthemes::theme_ipsum_tw(),
  ggstatsplot.layer = FALSE,
  messages = FALSE
)
ggscatterstats(
  data = movies_long,
  x = budget,
  y = rating,
  type = "p", # default #<<<
  conf.level = 0.99,
  marginal=F,
  messages = TRUE
)
