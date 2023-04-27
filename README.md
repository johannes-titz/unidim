Analysis for Titz: Are psychologists constructing their scales in the
right way…
================
Johannes Titz
04/27/2023

- [Part I: Theory](#part-i-theory)
  - [Prepare](#prepare)
  - [Generate Data](#generate-data)
  - [Plot](#plot)
  - [FA model](#fa-model)
- [Part II: Empirical Example](#part-ii-empirical-example)
  - [Violations](#violations)
  - [unidim analysis](#unidim-analysis)
  - [FA](#fa)
  - [plot](#plot-1)

This analysis accompanies the paper: Titz, J. Are psychologists
constructing their scales in the right way? Insights from studying
perfect unidimensionality.

# Part I: Theory

## Prepare

I use librarian to manage dependencies, so please install it if you do
not have it yet (uncomment the first line).

``` r
# install.packages("librarian")
library(librarian)
shelf(partitions, pbapply, tidyverse, psych, johannes-titz/zysno, xtable, Matrix, plot.matrix)
# use multiple cores
cores <- round(parallel::detectCores() * 0.8)
# create vectors from phi matrix
vec_from_tbl <- function(a, b, c, d) {
  mtrx <- matrix(c(a, b, c, d), ncol = 2)
  v1 <- rep(c(0, 1), colSums(mtrx))
  v2 <- rep(c(0, 1), rowSums(mtrx))
  cbind(v1, v2)
}
```

## Generate Data

``` r
n <- 200
# all combos, even the ones where phi is NA
grid <- as.matrix(t(compositions(n = n, m = 4)))
colnames(grid) <- paste0("Var", 1:4)
# rowwise does not work with the following code, I do not know why
# multiple cores seem to give no advantage, so is not used (argument cl)
# we cache this manually as it takes some time
grid <- xfun::cache_rds({
  grid <- as_tibble(grid) %>%
    rowwise() %>%
    mutate(errors = Var3,
           expected_errors = (Var1 + Var3) * (Var3 + Var4) / (Var1+ Var2+ Var3+ Var4),
           h = 1 - errors / expected_errors)
  phi <- pbapply(grid[, 1:4], 1, function(x) phi(x, digits = 4))
  cbind(grid, phi)
}, file = "sim1", rerun = F)
grid <- as_tibble(grid)
```

``` r
grid_unidim <-  grid %>%
  filter(Var3 == 0) # only unidimensional
grid_row <- grid %>%
  filter(phi >= 0, Var2 >= Var3) # error cell is b
```

## Plot

``` r
p <- ggplot(grid_row, aes(phi, h)) + 
  #geom_point(alpha = 0.05) +
  geom_hex(bins = 100, show.legend = FALSE) + 
  theme_classic() +
  scale_x_continuous("Phi", breaks = seq(0, 1, 0.1)) + 
  scale_y_continuous("H", breaks = seq(0, 1, 0.1) ) + 
  coord_fixed()
p
```

![](README_files/figure-gfm/phihs-1.png)<!-- -->

``` r
pdf("plots/phiH.pdf", width = 5, height = 5)
p
dev.off()
```

    ## png 
    ##   2

## FA model

It is not hard to create unidimensional model that results in a perfect
bifactor solution with factor analysis (for Pearson Correlation):

``` r
m1 <- vec_from_tbl(10, 80, 0, 10)
m2 <- vec_from_tbl(11, 78, 0, 11)
m3 <- vec_from_tbl(15, 70, 0, 15)

d <- cbind(m1, m2, m3[, -2])
# check that it is really unidimensional
zysnotize(d)
```

    ## $error_matrix
    ##      [,1] [,2] [,3] [,4] [,5]
    ## [1,]    0    0    0    0    0
    ## [2,]    0    0    0    0    0
    ## [3,]    0    0    0    0    0
    ## [4,]    0    0    0    0    0
    ## [5,]    0    0    0    0    0
    ## 
    ## $expected_error_matrix
    ##        [,1]   [,2]     [,3]     [,4]     [,5]
    ## [1,]  81.00  81.00  88.1100  88.1100 114.7500
    ## [2,]  81.00  81.00  88.1100  88.1100 114.7500
    ## [3,]  88.11  88.11  95.8441  95.8441 124.8225
    ## [4,]  88.11  88.11  95.8441  95.8441 124.8225
    ## [5,] 114.75 114.75 124.8225 124.8225 162.5625
    ## 
    ## $scalability_matrix
    ##      [,1] [,2] [,3] [,4] [,5]
    ## [1,]    1    1    1    1    1
    ## [2,]    1    1    1    1    1
    ## [3,]    1    1    1    1    1
    ## [4,]    1    1    1    1    1
    ## [5,]    1    1    1    1    1
    ## 
    ## $scalability
    ## [1] 1
    ## 
    ## $sum_errors
    ## [1] 0
    ## 
    ## $sum_expected_errors
    ## [1] 2016.858

``` r
ncol <- ncol(d)
colnames(d) <- paste("Item", seq(ncol(d)))
f1 <- factanal(d, 1)
f2 <- factanal(d, 2)

f1$loadings
```

    ## 
    ## Loadings:
    ##        Factor1
    ## Item 1 0.950  
    ## Item 2 0.119  
    ## Item 3 0.997  
    ## Item 4 0.125  
    ## Item 5 0.839  
    ## 
    ##                Factor1
    ## SS loadings      2.632
    ## Proportion Var   0.526

``` r
f2$loadings
```

    ## 
    ## Loadings:
    ##        Factor1 Factor2
    ## Item 1 0.950          
    ## Item 2         0.946  
    ## Item 3 0.997          
    ## Item 4         0.993  
    ## Item 5 0.837          
    ## 
    ##                Factor1 Factor2
    ## SS loadings      2.617   1.886
    ## Proportion Var   0.523   0.377
    ## Cumulative Var   0.523   0.900

``` r
print(xtable::xtable(cbind(f1$loadings, f2$loadings), label = "tab:fa1b", 
                     caption = "Cross tables of five unidimensional items"),
      file = "tables/fa1b.tex", booktabs = TRUE)

combs <- Map(function(x, y) table(d[, x], d[, y]),
             rep(1:ncol, each = ncol), rep(1:ncol, ncol))
lower <- seq(1, 21, 5)
upper <- seq(5, 25, 5)
res <- Map(function(x, y) Reduce(cbind, combs[x:y]),
    lower, upper)
res2 <- Reduce(rbind, res)
colnames(res2) <- paste0("$", rep(1:5, each = 2), "_", c(0, 1), "$")
rownames(res2) <- paste0("$", rep(1:5, each = 2), "_", c(0, 1), "$")

res2
```

    ##       $1_0$ $1_1$ $2_0$ $2_1$ $3_0$ $3_1$ $4_0$ $4_1$ $5_0$ $5_1$
    ## $1_0$    90     0    10    80    89     1    11    79    85     5
    ## $1_1$     0    10     0    10     0    10     0    10     0    10
    ## $2_0$    10     0    10     0    10     0    10     0    10     0
    ## $2_1$    80    10     0    90    79    11     1    89    75    15
    ## $3_0$    89     0    10    79    89     0    11    78    85     4
    ## $3_1$     1    10     0    11     0    11     0    11     0    11
    ## $4_0$    11     0    10     1    11     0    11     0    11     0
    ## $4_1$    79    10     0    89    78    11     0    89    74    15
    ## $5_0$    85     0    10    75    85     0    11    74    85     0
    ## $5_1$     5    10     0    15     4    11     0    15     0    15

``` r
print(xtable::xtable(res2, label = "tab:fa1",
                     caption = "Example of factor analysis for five unidimensional items"),
      booktabs = TRUE,
      file = "tables/fa1.tex",
      sanitize.text.function = function(x) {x})
```

# Part II: Empirical Example

Load data

``` r
d <- read.csv2("data.csv")
```

Find order, calculate socre, show df

``` r
o <- order(colMeans(d[,-1]), decreasing = TRUE)
d <- d[, c(1, o+1)]
d$score <- rowSums(d[, -1])
df <- as.data.frame(d)
df[order(df$score), ]
```

    ##      id X1 X4 X3 X2 X5 score
    ## 13   13  0  0  0  0  0     0
    ## 22   22  0  0  0  0  0     0
    ## 46   46  0  0  0  0  0     0
    ## 51   51  0  0  0  0  0     0
    ## 54   54  0  0  0  0  0     0
    ## 57   57  0  0  0  0  0     0
    ## 64   64  0  0  0  0  0     0
    ## 66   66  0  0  0  0  0     0
    ## 73   73  0  0  0  0  0     0
    ## 80   80  0  0  0  0  0     0
    ## 85   85  0  0  0  0  0     0
    ## 94   94  0  0  0  0  0     0
    ## 97   97  0  0  0  0  0     0
    ## 99   99  0  0  0  0  0     0
    ## 101 101  0  0  0  0  0     0
    ## 104 104  0  0  0  0  0     0
    ## 116 116  0  0  0  0  0     0
    ## 121 121  0  0  0  0  0     0
    ## 126 126  0  0  0  0  0     0
    ## 133 133  0  0  0  0  0     0
    ## 17   17  1  0  0  0  0     1
    ## 40   40  1  0  0  0  0     1
    ## 41   41  0  0  0  0  1     1
    ## 70   70  1  0  0  0  0     1
    ## 108 108  0  0  1  0  0     1
    ## 2     2  1  1  0  0  0     2
    ## 5     5  1  1  0  0  0     2
    ## 7     7  1  1  0  0  0     2
    ## 8     8  0  0  1  1  0     2
    ## 9     9  1  1  0  0  0     2
    ## 18   18  1  1  0  0  0     2
    ## 26   26  1  1  0  0  0     2
    ## 29   29  1  1  0  0  0     2
    ## 37   37  1  1  0  0  0     2
    ## 56   56  1  1  0  0  0     2
    ## 65   65  1  1  0  0  0     2
    ## 74   74  1  1  0  0  0     2
    ## 76   76  0  0  1  1  0     2
    ## 81   81  1  1  0  0  0     2
    ## 86   86  1  1  0  0  0     2
    ## 89   89  1  0  1  0  0     2
    ## 93   93  1  1  0  0  0     2
    ## 100 100  1  1  0  0  0     2
    ## 107 107  1  1  0  0  0     2
    ## 113 113  1  1  0  0  0     2
    ## 120 120  1  0  0  1  0     2
    ## 123 123  1  1  0  0  0     2
    ## 15   15  1  1  1  0  0     3
    ## 34   34  1  1  1  0  0     3
    ## 36   36  1  1  0  1  0     3
    ## 45   45  1  1  1  0  0     3
    ## 47   47  1  1  1  0  0     3
    ## 61   61  1  1  1  0  0     3
    ## 83   83  1  1  1  0  0     3
    ## 109 109  1  1  1  0  0     3
    ## 112 112  1  1  1  0  0     3
    ## 117 117  1  1  0  0  1     3
    ## 128 128  1  1  1  0  0     3
    ## 1     1  1  1  1  1  0     4
    ## 4     4  1  1  1  1  0     4
    ## 10   10  1  1  1  1  0     4
    ## 12   12  1  1  1  1  0     4
    ## 14   14  1  1  1  1  0     4
    ## 20   20  1  1  1  1  0     4
    ## 23   23  1  1  1  1  0     4
    ## 24   24  1  1  1  1  0     4
    ## 25   25  1  1  1  1  0     4
    ## 30   30  1  1  1  1  0     4
    ## 31   31  1  1  1  1  0     4
    ## 32   32  1  1  1  1  0     4
    ## 33   33  1  1  1  1  0     4
    ## 35   35  1  1  1  1  0     4
    ## 42   42  1  1  1  1  0     4
    ## 43   43  1  1  1  1  0     4
    ## 49   49  1  1  1  1  0     4
    ## 53   53  1  1  1  1  0     4
    ## 55   55  1  1  1  1  0     4
    ## 62   62  1  1  1  1  0     4
    ## 63   63  1  1  1  1  0     4
    ## 68   68  1  1  1  1  0     4
    ## 79   79  1  1  1  1  0     4
    ## 84   84  1  1  1  1  0     4
    ## 87   87  1  1  1  1  0     4
    ## 95   95  1  1  1  1  0     4
    ## 98   98  1  1  1  1  0     4
    ## 102 102  1  1  1  1  0     4
    ## 103 103  1  1  1  1  0     4
    ## 106 106  1  1  1  1  0     4
    ## 110 110  1  1  1  1  0     4
    ## 118 118  1  1  1  1  0     4
    ## 124 124  1  1  1  1  0     4
    ## 127 127  1  1  1  1  0     4
    ## 130 130  1  1  1  1  0     4
    ## 131 131  1  1  1  1  0     4
    ## 3     3  1  1  1  1  1     5
    ## 6     6  1  1  1  1  1     5
    ## 11   11  1  1  1  1  1     5
    ## 16   16  1  1  1  1  1     5
    ## 19   19  1  1  1  1  1     5
    ## 21   21  1  1  1  1  1     5
    ## 27   27  1  1  1  1  1     5
    ## 28   28  1  1  1  1  1     5
    ## 38   38  1  1  1  1  1     5
    ## 39   39  1  1  1  1  1     5
    ## 44   44  1  1  1  1  1     5
    ## 48   48  1  1  1  1  1     5
    ## 50   50  1  1  1  1  1     5
    ## 52   52  1  1  1  1  1     5
    ## 58   58  1  1  1  1  1     5
    ## 59   59  1  1  1  1  1     5
    ## 60   60  1  1  1  1  1     5
    ## 67   67  1  1  1  1  1     5
    ## 69   69  1  1  1  1  1     5
    ## 71   71  1  1  1  1  1     5
    ## 72   72  1  1  1  1  1     5
    ## 75   75  1  1  1  1  1     5
    ## 77   77  1  1  1  1  1     5
    ## 78   78  1  1  1  1  1     5
    ## 82   82  1  1  1  1  1     5
    ## 88   88  1  1  1  1  1     5
    ## 90   90  1  1  1  1  1     5
    ## 91   91  1  1  1  1  1     5
    ## 92   92  1  1  1  1  1     5
    ## 96   96  1  1  1  1  1     5
    ## 105 105  1  1  1  1  1     5
    ## 111 111  1  1  1  1  1     5
    ## 114 114  1  1  1  1  1     5
    ## 115 115  1  1  1  1  1     5
    ## 119 119  1  1  1  1  1     5
    ## 122 122  1  1  1  1  1     5
    ## 125 125  1  1  1  1  1     5
    ## 129 129  1  1  1  1  1     5
    ## 132 132  1  1  1  1  1     5

### Violations

``` r
is_monotonic <- function(x) {
  all(diff(as.numeric(x)) <= 0)
}

df$monotonic <- apply(df[, 2:6], 1, is_monotonic)
dg <- df[df$monotonic == FALSE, ]

dh <- dg[order(dg$score), ]
dh
```

    ##      id X1 X4 X3 X2 X5 score monotonic
    ## 41   41  0  0  0  0  1     1     FALSE
    ## 108 108  0  0  1  0  0     1     FALSE
    ## 8     8  0  0  1  1  0     2     FALSE
    ## 76   76  0  0  1  1  0     2     FALSE
    ## 89   89  1  0  1  0  0     2     FALSE
    ## 120 120  1  0  0  1  0     2     FALSE
    ## 36   36  1  1  0  1  0     3     FALSE
    ## 117 117  1  1  0  0  1     3     FALSE

``` r
print(
      xtable::xtable(dh[, 2:7], caption = "Participant answers that violate unidimensionality", 
                     label = "tab:viol", digits = 0),
      booktabs = TRUE,
      file = "tables/viol.tex"
)
```

### unidim analysis

``` r
lv <- loevenize(as.matrix(d[, 2:6]))
lv
```

    ## $error_matrix
    ##    col
    ## row 1 2 3 4 5
    ##   1 0 0 3 2 1
    ##   2 0 0 4 3 1
    ##   3 3 4 0 2 2
    ##   4 2 3 2 0 2
    ##   5 1 1 2 2 0
    ## 
    ## $expected_error_matrix
    ##    col
    ## row         1        2        3        4         5
    ##   1  0.000000 18.76692 15.87970 14.25564  7.398496
    ##   2 18.766917  0.00000 19.18797 17.22556  8.939850
    ##   3 15.879699 19.18797  0.00000 26.72932 13.872180
    ##   4 14.255639 17.22556 26.72932  0.00000 16.646617
    ##   5  7.398496  8.93985 13.87218 16.64662  0.000000
    ## 
    ## $h_matrix
    ##    col
    ## row         1         2         3         4         5
    ##   1 0.0000000 1.0000000 0.8110795 0.8597046 0.8648374
    ##   2 1.0000000 0.0000000 0.7915361 0.8258402 0.8881413
    ##   3 0.8110795 0.7915361 0.0000000 0.9251758 0.8558266
    ##   4 0.8597046 0.8258402 0.9251758 0.0000000 0.8798555
    ##   5 0.8648374 0.8881413 0.8558266 0.8798555 0.0000000
    ## 
    ## $h
    ## [1] 0.8741365
    ## 
    ## $sum_errors
    ## [1] 20
    ## 
    ## $sum_expected_errors
    ## [1] 158.9023

### FA

``` r
fa1 <- factanal(d[,c(-1, -7)], 1)
fa2 <- factanal(d[,c(-1, -7)], 2)

tab <- cbind(unclass(loadings(fa1)), unclass(loadings(fa2)))
print(
      xtable::xtable(tab, caption = "Factor analysis",
                     label = "tab:fa2", digits = 3),
      booktabs = TRUE,
      file = "tables/fa2.tex"
)
tab
```

    ##      Factor1   Factor1   Factor2
    ## X1 0.9176135 0.3067621 0.8365763
    ## X4 0.9596311 0.3224612 0.9439382
    ## X3 0.6267584 0.7436943 0.3651334
    ## X2 0.5772825 0.9632071 0.2297556
    ## X5 0.3429633 0.4628628 0.1732022

### plot

``` r
tab <- round(cor(d[, -c(1, 7)]), 2)

a <- tril(lv$h_matrix, -1) # strict lower triangular matrix (omit diagonals)
b <- triu(tab, 1) # strict upper triangular matrix
c <- a + b
diag(c) <- NA

pdf("plots/itempairs.pdf", width = 7, height = 7)
plot(as.matrix(c), fmt.cell='%.2f', main = "", key = NULL,
xlab = "Item", ylab = "Item", col = sapply(seq(1, 0, -0.1), gray, alpha = 0.5))
dev.off()
```

    ## png 
    ##   2

``` r
plot(as.matrix(c), fmt.cell='%.2f', main = "", key = NULL,
xlab = "Item", ylab = "Item", col = sapply(seq(1, 0, -0.1), gray, alpha = 0.5))
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->
