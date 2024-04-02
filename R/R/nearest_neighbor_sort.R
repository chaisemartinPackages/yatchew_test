#' Internal function to generate a sorting id based on euclidean distance nearest neighbors.
#' @param df A dataframe
#' @param D Independent variables vector
#' @importFrom stats quantile
#' @returns A dataframe with a new id variable
#' @noRd
nearest_neighbor_sort <- function(df, D) {
    n_cutoff <- ceiling(log(nrow(df), 10)) * 2
    cutoffs <- 1:n_cutoff / n_cutoff
    for (v in 1:length(D)) {
        max_var <- max(df[[D[v]]], na.rm = TRUE)
        min_var <- min(df[[D[v]]], na.rm = TRUE)
        df[[paste0(D[v],"_norm")]] <- (df[[D[v]]] - min_var) / (max_var - min_var)

        df[[paste0("cube_dim_",v)]] <- 1
        for (j in 1:length(cutoffs)) {
             df[[paste0("cube_dim_",v)]] <- ifelse(df[[paste0(D[[v]],"_norm")]] > quantile(df[[paste0(D[[v]],"_norm")]], cutoffs[j], type = 2),  df[[paste0("cube_dim_",v)]] + 1 ,  df[[paste0("cube_dim_",v)]])
        }
    }

    df$cube_var <- ""
    # Adjustments #
    for (v in 1:length(D)) {
        if (v > 1) {
            df[[paste0("cube_dim_",v)]] <- ifelse(df[[paste0("cube_dim_",v-1)]] %% 2 == 0, n_cutoff - df[[paste0("cube_dim_",v)]] + 1 , df[[paste0("cube_dim_",v)]])
        }
        df[[paste0("cube_dim_",v,"_str")]] <- sapply(df[[paste0("cube_dim_",v)]], function(x) intToUtf8(x + 64))
        df$cube_var <- paste0(df$cube_var, df[[paste0("cube_dim_",v,"_str")]])
    }

    for (v in 1:length(D)) {
        df[[paste0("cube_dim_",v)]] <- NULL
        df[[paste0("cube_dim_",v,"_str")]] <- NULL
    }
    df$rown <- 1:nrow(df)
    df <- df[order(df$cube_var, df$rown), ]

    df_sort <- NULL
    for (l in levels(factor(df$cube_var))) {
        df_temp <- subset(df, df$cube_var == l)
        mat <- matrix(NA, nrow = nrow(df_temp), ncol = length(D))
        for (v in 1:length(D)) {
            mat[,v] <- df_temp[[paste0(D[v],"_norm")]]
        }
        df_temp <- df_temp[msort(mat), ]
        df_sort <- rbind(df_sort, df_temp)
    }

    rownames(df_sort) <- 1:nrow(df_sort)
    return(df_sort)
}
