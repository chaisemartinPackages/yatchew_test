# yatchew_test

Yatchew (1997), de Chaisemartin and D'Haultfoeuille (2024) linearity test.

## Setup

### Stata (SSC)
```r
ssc install yatchew_test, replace
```

### R (CRAN)
```r
install.packages("YatchewTest")

```

## Syntax
### Stata 
```stata
yatchew_test varlist(min = 2 numeric) [if] [, het_robust path_plot order(#)]
```

### R (Method dispatch for generic data.frame objects)
```r
yatchew_test(data, Y, D, het_robust = FALSE, path_plot = FALSE, order = 1)
```

## Overview

This program implements the linearity test proposed by Yatchew (1997) and its heteroskedasticity-robust version proposed by de Chaisemartin and D'Haultfoeuille (2024).  In this overview, we sketch the intuition behind the two tests, as to motivate the use of the package and its options. Please refer to Yatchew (1997) and Section 3 of de Chaisemartin and D'Haultfoeuille (2024) for further details.

Yatchew (1997) proposes a useful extension of the test with multiple independent variables. The program implements this extension when the *D* argument has length $> 1$. It should be noted that the power and consistency of the test in the multivariate case are not backed by proven theoretical results. We implemented this extension to allow for testing and exploratory research. Future theoretical exploration of the multivariate test will depend on the demand and usage of the package.

### Univariate Yatchew Test
Let $Y$ and $D$ be two random variables. Let $m(D) = E[Y|D]$. The null hypothesis of the test is that $m(D) = \alpha_0 + \alpha_1 D$ for two real numbers $\alpha_0$ and $\alpha_1$. This means that, under the null, $m(.)$ is linear in $D$. The outcome variable can be decomposed as $Y = m(D) + \varepsilon$, with $E[\varepsilon|D] = 0$ and $\Delta Y = \Delta \varepsilon$ for $\Delta D \to 0$.  In a dataset with $N$ i.i.d. realisations of $(Y, D)$, one can test this hypothesis as follows: 
1. sort the dataset by $D$;
2. denote the corresponding observations by $(Y_{(i)}, D_{(i)})$, with $i \in \lbrace 1, ..., N\rbrace$;
3. approximate $\hat{\sigma}^2_{\text{diff}}$, i.e. the variance of the first differenced residuals $\varepsilon_{(i)} - \varepsilon_{(i-1)}$, by the variance of $Y_{(i)} - Y_{(i-1)}$;
4. compute $\hat{\sigma}^2_{\text{lin}}$, i.e. the variance of the residuals from an OLS regression of $Y$ on $D$. 

Heuristically, the validity of step (3) derives from the fact that $Y_{(i)} - Y_{(i-1)}$ = $m(D_{(i)}) - m(D_{(i-1)})$ + $\varepsilon_{(i)} - \varepsilon_{(i-1)}$ and the first difference term is close to zero for $D_{(i)} \approx D_{(i-1)}$. Sorting at step (1) ensures that consecutive $D_{(i)}s$  are as close as possible, and when the sample size goes to infinity the distance between consecutive observations goes to zero. Then, Yatchew (1997) shows that under homoskedasticity and regularity conditions $$T := \sqrt{G}\left(\dfrac{\hat{\sigma}^2_{\text{lin}}}{\hat{\sigma}^2_{\text{diff}}}-1\right) \stackrel{d}{\longrightarrow} \mathcal{N}\left(0,1\right).$$
Then, one can reject the linearity of $m(.)$ with significance level $\alpha$ if $T > \Phi(1-\alpha)$. 

If the homoskedasticity assumption fails, this test leads to overrejection. De Chaisemartin, D'Haultfoeuille & Gurgand (2024) propose a heteroskedasticity-robust version of the test statistic above. This version of the Yatchew (1997) test can be implemented by running the command with the option het_robust (Stata) or het_robust = TRUE (R).

## Multivariate Yatchew Test
Let $\textbf{D}$ be a vector of $K$ random variables. Let $g(\textbf{D}) = E[Y|\textbf{D}]$. Denote with $||.,.||$ the Euclidean distance between two vectors. The null hypothesis of the multivariate test is $g(\textbf{D}) = \alpha_0 + A'\textbf{D}$, with $A = (\alpha_1,..., \alpha_K)$, for $K+1$ real numbers $\alpha_0$, $\alpha_1$, ..., $\alpha_K$. This means that, under the null, $g(.)$ is linear in $\textbf{D}$. Following the same logic as the univariate case, in a dataset with $N$ i.i.d. realisations of $(Y, \textbf{D})$ we can approximate the first difference $\Delta \varepsilon$ by $\Delta Y$ valuing $g(.)$ between consecutive observations. The program runs a nearest neighbor algorithm to find the sequence of observations such that the Euclidean distance between consecutive positions is minimized. The algorithm has been programmed in C++ and it has been integrated in Stata via OS-specific plugins and in R thanks to the Rcpp library. 

The program follows a very simple nearest neighbor approach:
1. collect all the Euclidean distances between all the possible unique pairs of rows in $\textbf{D}$ in the matrix $M$, where $M_{n,m} = ||\textbf{D}_n,\textbf{D}_m||$ with $n,m \in \lbrace 1, ..., N\rbrace$;
2. setup the queue to $Q = \lbrace 1, ..., N\rbrace$, the (empty) path vector $I = \lbrace\rbrace$ and the starting index $i = 1$;
3. remove $i$ from $Q$ and find the column index $j$ of M such that $M_{i,j} = \min_{c \in Q} M_{i,c}$;
4. append $j$ to $I$ and start again from step 3 with $i = j$ until $Q$ is empty.

To improve efficiency, the program collects only the $N(N-1)/2$ Euclidean distances corresponding to the lower triangle of matrix $M$ and chooses $j$ such that $M_{i,j} = \min_{c \in Q} 1\lbrace c < i\rbrace M_{i,c} + 1\lbrace c > i\rbrace M_{c,i}$. The output of the algorithm, i.e. the vector $I$, is a sequence of row numbers such that the distance between the corresponding rows $\textbf{D}_{i}s$ is minimized. The program also uses two refinements suggested in Appendix A of Yatchew (1997):
* The entries in $\textbf{D}$ are normalized in $[0,1]$;
* The algorithm is applied to sub-cubes, i.e. partitions of the $[0,1]^K$ space, and the full path is obtained by joining the extrema of the subpaths.

By convention, the program computes $(2\lceil \log_{10} N \rceil)^K$ subcubes, where each univariate partition is defined by grouping observations in $2\lceil \log_{10} N \rceil$ quantile bins. If $K = 2$, the user can visualize in a ggplot graph the exact path across the normalized $\textbf{D}_{i}s$ by running the command with the option path_plot (Stata) or path_plot = TRUE (R).

Once the dataset is sorted by $I$, the program resumes from step (2) of the univariate case.

## Options

**data**: (data.frame, R only) A dataframe.

**Y**: (char, R only) Dependent variable. In Stata, this argument is the first varname in varlist.

**D**: (char, R only) Independent variable(s). In Stata, this argument is the collection of varnames after the first in varlist.

**het_robust**: (logical) If FALSE, the test is performed under the assumption of homoskedasticity (Yatchew, 1997). If TRUE, the test is performed using the heteroskedasticity-robust test statistic proposed by de Chaisemartin and D'Haultfoeuille (2024).

**path_plot**: (logical) if TRUE and argument D has length 2, the assigned object will include a plot of the sequence of $(D_{1i}, D_{2i})$ that minimizes the euclidean distance between each pair of consecutive observations (see Overview for further details).

**order**: (nonnegative integer $k$) If this option is specified, the program tests whether the conditional expectation of $Y$ given $D$ is a linear function of a $k$-degree polynomial in $D$. The command tests the hypothesis that the conditional mean of $Y$ given $D$ is constant whenever `order = 0` is specified.

## References 

de Chaisemartin, C., D'Haultfoeuille, X. ([2024](https://ssrn.com/abstract=4284811)). Two-way Fixed Effects and Difference-in-Difference Estimators in Heterogeneous Adoption Designs.

Yatchew, A. ([1997](https://doi.org/10.1016/S0165-1765(97)00218-8)). An elementary estimator of the partial linear model.
