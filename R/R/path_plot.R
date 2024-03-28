#' Internal function to render a ggplot with the continuous path
#' @param df A dataframe
#' @param D A vector with two vars
#' @import ggplot2
#' @returns A ggplot object
path_plot <- function(df, D) {
    x_base <- df[[paste0(D[1],"_norm")]][c(1,nrow(df))]
    y_base <- df[[paste0(D[2],"_norm")]][c(1,nrow(df))]
    map <- ggplot(data = df, aes(x = .data[[paste0(D[1],"_norm")]], y = .data[[paste0(D[2],"_norm")]])) + 
    geom_path() +
    geom_point(aes(x = x_base[1], y = y_base[1]), size = 2, color = "green") + 
    geom_point(aes(x = x_base[2], y = y_base[2]), size = 2, color = "red") +
    labs(caption = "First node in green, last node in red.")
    return(map)
}