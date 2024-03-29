#' Main function 
#' @md
#' @description Heteroskedasticity Robust Test of Linearity (Yatchew, 1997; de Chaisemartin and D'Haultfoeuille, 2024)
#' @param data A data object.
#' @param ... Undocumented.
#' @returns Method dispatch depending on the data object class. 
#' @export
yatchew_test <- function(data, ...) {
        UseMethod("yatchew_test")
    }

#' @md
#' @title General yatchew_test method for unclassed dataframes
#' @name yatchew_test.data.frame
#' @description General yatchew_test method for unclassed dataframes
#' @param data (data.frame) A dataframe.
#' @param Y (char) Dependent variable.
#' @param D (char) Independent variable.
#' @param het_robust (logical) If FALSE, the test is performed under the assumption of homoskedasticity (Yatchew, 1997). If TRUE, the test is performed using the heteroskedasticity-robust test statistic proposed by de Chaisemartin and D'Haultfoeuille (2024).
#' @param path_plot (logical) if TRUE, the assigned object will include a plot of the sequence of \eqn{(D_{1i}, D_{2i})}s that minimizes the euclidean distance between each pair of consecutive observations (see Overview for further details).
#' @param ... Undocumented.
#' @section Overview:
#' This program implements the linearity test proposed by Yatchew (1997) and its heteroskedasticity-robust version proposed by de Chaisemartin & D'Haultfoeuille (2024). 
#' In this overview, we sketch the intuition behind the two tests, as to motivate the use of the package and its options. 
#' Please refer to Yatchew (1997) and Section 3 of de Chaisemartin & D'Haultfoeuille (2024) for further details.
#' 
#' Yatchew (1997) proposes a useful extension of the test with multiple independent variables. 
#' The program implements this extension when the \code{D} argument has length \eqn{> 1}.
#' It should be noted that the power and consistency of the test in the multivariate case are not backed by proven theoretical results.
#' We implemented this extension to allow for testing and exploratory research.
#' Future theoretical exploration of the multivariate test will depend on the demand and usage of the package.
#' 
#' ## Univariate Yatchew Test
#' Let \eqn{Y} and \eqn{D} be two random variables with continuous and bounded support. 
#' The core function of the \code{yatchew_test()} program is to check that the \eqn{Y} argument is linear in the \code{D} argument. 
#' The null hypothesis of the test is \eqn{E[Y|D] = m(D)}, where \eqn{m(.)} is a continuous linear function. 
#' The outcome variable can be decomposed as \eqn{Y = m(D) + \varepsilon}, with \eqn{E[\varepsilon|D] = 0}. 
#' If the null holds, \eqn{\Delta Y = \Delta \varepsilon} for \eqn{\Delta D \to 0}. 
#' In a dataset with \eqn{N} sample realisations of \eqn{Y} and \eqn{D}, one can test this hypothesis as follows: 
#' 1. sort the dataset by \eqn{D};
#' 2. denote the corresponding observations by \eqn{(Y_{(i)}, D_{(i)})}, with \eqn{g \in \{1, ..., N\}};
#' 3. compute \eqn{\hat{\sigma}^2_{\text{diff}}}, i.e. the variance of the residuals \eqn{\varepsilon_{(i)} - \varepsilon_{(i-1)}} = \eqn{Y_{(i)} - Y_{(i-1)}};
#' 4. compute \eqn{\hat{\sigma}^2_{\text{lin}}}, i.e. the variance of the residuals from an OLS regression of \eqn{Y} on \eqn{D}. 
#' 
#' The feasibility of step (3) derives from the fact that \eqn{Y_{(i)} - Y_{(i-1)}} = \eqn{m(D_{(i)}) - m(D_{(i-1)})} + \eqn{\varepsilon_{(i)} - \varepsilon_{(i-1)}} and, under the null, the first difference term is close to zero for \eqn{D_{(i)} \approx D_{(i-1)}}. 
#' Sorting at step (1) ensures that consecutive \eqn{D_{(i)}}s are as close as possible.
#' Lastly, Yatchew (1997) shows that under homoskedasticity and regularity conditions
#' 
#' \deqn{T := \sqrt{G}\left(\dfrac{\hat{\sigma}^2_{\text{lin}}}{\hat{\sigma}^2_{\text{diff}}}-1\right) \stackrel{d}{\longrightarrow} \mathcal{N}\left(0,1\right)}
#' 
#' To this end, one can reject the linearity of \eqn{m(.)} with significance level \eqn{\alpha} if \eqn{T > \Phi(1-\alpha)}. 
#' 
#' If the homoskedasticity assumption fails, this test leads to overrejection. 
#' De Chaisemartin & D'Haultfoeuille (2024) propose a heteroskedasticity-robust version of the test statistic above. 
#' This version of the Yatchew (1997) test can be implemented by running the command with the option \code{het_robust = TRUE}.
#' 
#' ## Multivariate Yatchew Test
#' Let \eqn{\textbf{D}} is a \eqn{N \times M} random matrix and denote with \eqn{||.,.||} the Euclidean distance between two vectors.
#' The new null hypothesis of the test is \eqn{E[Y|\textbf{D}] = g(\textbf{D})}, where \eqn{g:\mathbb{R}^M \to \mathbb{R}} is a continuous linear function. 
#' Following the same logic as the univariate case, we want to reduce \eqn{\Delta Y} to \eqn{\Delta \varepsilon} by valuing \eqn{g(.)} between consecutive observations.
#' However, multivariate sorting does not ensure that \eqn{||\textbf{D}_{(i)},\textbf{D}_{(i-1)}|| \to 0}. 
#' To this end, the program runs a nearest neighbor algorithm to find the sequence of observations such that the Euclidean distance between consecutive positions is minimized.
#' The algorithm has been programmed in C++ and it has been integrated in R thanks to the \code{Rcpp} library. The program follows a very simple nearest neighbor approach:
#' 1. collect all the Euclidean distances between all the possible unique pairs of rows in \eqn{\textbf{D}} in the matrix \eqn{M}, where \eqn{M_{n,m} = ||\textbf{D}_n,\textbf{D}_m||} with \eqn{n,m \in \{1, ..., N\}} row indices of \eqn{\textbf{D}};
#' 2. setup the queue to \eqn{Q = \{1, ..., N\}}, the (empty) path vector \eqn{I = \{\}} and the starting index \eqn{i = 1};
#' 3. remove \eqn{i} from \eqn{Q} and find the column index \eqn{j} of M such that \eqn{M_{i,j} = \min_{c \in Q} M_{i,c}};
#' 4. append \eqn{j} to \eqn{I} and start again from step 3 with \eqn{i = j} until \eqn{Q} is empty.
#' 
#' To improve efficiency, the program collects only the \eqn{N(N-1)/2} Euclidean distances corresponding to the lower triangle of matrix \eqn{M} and chooses \eqn{j} as \eqn{\min_{c \in Q} 1\{c < i\} M_{i,c} + 1\{c > i\} M_{c,i}}.
#' The output of the algorithm is a sequence of row numbers such that the distance between the corresponding rows \eqn{\textbf{D}_i}s is minimized.
#' The program also uses two refinements suggested in Appendix A of Yatchew (1997):
#' * The entries in \eqn{\textbf{D}} are normalized in \eqn{[0,1]};
#' * The algorithm above is applied to sub-cubes, i.e. partitions of the \eqn{[0,1]^M} space, and the full path is obtained by joining the extrema of the subpaths.
#' 
#' By convention, the program computes \eqn{(2\lceil \log_{10} N \rceil)^M} subcubes, where each univariate partition is defined by grouping observations in \eqn{2\lceil \log_{10} N \rceil} quantile bins. If \eqn{M = 2}, the user can visualize in a ggplot graph the exact path across the normalized \eqn{\textbf{D}_i}s by running the command with the option \code{path_plot = TRUE}.
#' 
#' The vector \eqn{I} at the end of the algorithm is finally used to sort the dataset and the program resumes from step (2) of the univariate case.
#' @references 
#' de Chaisemartin, C., d'Haultfoeuille, X., Gurgand, M. (2024). Two-way Fixed Effects and Difference-in-Difference Estimators in Heterogeneous Adoption Designs.
#' 
#' Yatchew, A. (1997). An elementary estimator of the partial linear model.
#' 
#' @examples
#' df <- as.data.frame(matrix(NA, nrow = 1E3, ncol = 0))
#' df$x <- rnorm(1E3)
#' df$b <- runif(1E3)
#' df$y <- 2 + df$b * df$x
#' yatchew_test(data = df, Y = "y", D = "x")
#' @returns A list with test results.
#' @importFrom stats as.formula lm pnorm var
#' @export
yatchew_test.data.frame <- function(data, Y, D, het_robust = FALSE, path_plot = FALSE, ...){    
    args <- list()
    for (v in names(formals(yatchew_test))) {
        if (!(v %in% c("data", "..."))) {
            args[[v]] <- get(v)
        }
    }
    res <- list(args = args)
    grout <- NULL

    if (isTRUE(path_plot) & length(D) != 2) {
        stop("path_plot can be requested only with 2 treatment variables.")
    }

    data <- subset(data, !is.na(data[[Y]]))
    for (v in D) {
        data <- subset(data, !is.na(data[[v]]))
    }
    if (length(D) == 1) {
       data <- data[order(data[[D]], data[[Y]]), ]
       reg_formula <- paste0(Y,"~",D)
    } else {
        data <- nearest_neighbor_sort(data, D)
        if (length(D) == 2 & isTRUE(path_plot)) {
            grout <- path_plot(data, D)
        }
        reg_formula <- paste0(Y,"~",D[1])
        for (j in 1:length(D)) {
            reg_formula <- paste0(reg_formula,"+",D[j])
        }
    }

    # Variance of residuals from linear regression
    coef_lin <- lm(as.formula(reg_formula), data = data)$coefficients
    data$e_lin <- data[[Y]] - (coef_lin[1] + as.matrix(data[D], ncol = length(D)) %*% coef_lin[2:length(coef_lin)])
    var_lin <- var(data$e_lin, na.rm = TRUE)

    # Variance of residuals from nonparametric model
    data$e_diff <- data[[Y]] - c(NA, data[[Y]][1:(nrow(data)-1)])
    var_diff <- 0.5 * mean(data$e_diff ^ 2, na.rm = TRUE)

    if (isFALSE(het_robust)) {
        # Hypothesis test under Homoskedasticity
        T_test <- sqrt(nrow(data)) * ((var_lin/var_diff)-1)
    } else {
        # Hypothesis test under Heteroskedasticity
        num <- sqrt(nrow(data)) * (var_lin - var_diff)
        denom <- mean((data$e_lin[2:nrow(data)] * data$e_lin[1:(nrow(data)-1)])^2, na.rm = TRUE)
        T_test <- num/sqrt(denom)
    }
    results <- matrix(NA, nrow = 1, ncol = 5)
    results[1,1] <- var_lin
    results[1,2] <- var_diff
    results[1,3] <- T_test
    results[1,4] <- 1-pnorm(T_test)
    results[1,5] <- nrow(data)
    colnames(results) <- c("\U03C3\U00B2_lin", "\U03C3\U00B2_diff", "T", "p-value", "N")    
    rownames(results) <- ""
    res <- append(res, list(results))
    names(res)[length(res)] <- "results"
    if (!is.null(grout)) {
        res <- append(res, list(grout))
        names(res)[length(res)] <- "plot"        
    }
    class(res) <- "yatchew_test"
    return(res)
}