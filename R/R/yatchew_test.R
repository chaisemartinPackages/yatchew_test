#' Main function 
#' @md
#' @description Heteroskedasticity Robust Test of Linearity (Yatchew, 1997; de Chaisemartin and D'Haultfoeuille, 2024)
#' @param data (data.frame) A dataframe.
#' @param Y (char) Dependent variable.
#' @param D (char) Independent variable.
#' @param het_robust (logical) If FALSE, the test is performed under the assumption of homoskedasticity, i.e. \eqn{E[V(\varepsilon| X^2)] = \sigma^2}, as in Yatchew (1997). If TRUE, the test is performed using the heteroskedasticity-robust test statistic proposed by de Chaisemartin and D'Haultfoeuille (2024).
#' @param ... Undocumented.
#' @returns A list with test results. 
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
#' @export
yatchew_test <- function(data, Y, D, het_robust, ...) {
        UseMethod("yatchew_test")
    }

#' @title General yatchew_test method for unclassed dataframes
#' @name yatchew_test.data.frame
#' @description General yatchew_test method for unclassed dataframes
#' @param data (data.frame) A dataframe.
#' @param Y (char) Dependent variable.
#' @param D (char) Independent variable.
#' @param het_robust (logical) If FALSE, the test is performed under the assumption of homoskedasticity, i.e. \eqn{E[V(\varepsilon| X^2)] = \sigma^2}, as in Yatchow (1997). If TRUE, the test is performed using the heteroskedasticity-robust test statistic proposed by de Chaisemartin and D'Haultfoeuille (2024).
#' @param path_plot (logical)
#' @param ... Undocumented.
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