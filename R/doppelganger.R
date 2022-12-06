utils::globalVariables(".")
#' Find Doppelgangers
#' @name doppelganger
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @description Analyze correlations to select a subset of non-redundant variables.
#' @param data Data table source.
#' @param variables Variables included in \code{data} to consider.
#' @param priority Ranking method to prioritize variables. When "Centrality" is
#' selected, variables with higher mean absolute correlation will have higher
#' priority (see details). On the opposite, when "Peripherality" is selected,
#' the algorithm tries to keep the most peripheral variables, dropping variables
#' with elevated levels of correlation. "Raw order" prioritize variables using
#' their position within the vector \code{variables}.
#' @param threshold Correlation cut-off (absolute value) to identify doppelgangers.
#' Variables with absolute correlation values equal or greater to this value will be
#' considered aliases of each other.
#' @details The variable priority is established by calculating a centrality index
#' for each variable and sorting them accordingly. When \code{priority} is "Centrality",
#' the function sort variables from the most central to the most peripheral.
#' When \code{priority} is "Peripherality", the function sort variables from the most
#' peripheral to the most central. The centrality index is calculated as the
#' weighted mean of the absolute values of its correlations. Weights are obtained
#' by counting the data used to calculate each correlation. This weighting strategy
#' is applied in order to prioritize variables with less missing data or where
#' missing data match missing data of other variables.
#' @examples
#' \dontrun{
#' data(measures)
#' dg <- doppelganger(measures, threshold=0.7)
#' var_grp <- list(
#'     Selected = apply(sapply(dg$keep, "==", dg$variables), 2, which),
#'      Dropped = apply(sapply(dg$drop, "==", dg$variables), 2, which)
#' )
#' qgraph::qgraph(
#'     dg$cor_matrix, layout="spring", #vsize=2,
#'     minimum=0.3, edge.labels=TRUE,
#'     groups=var_grp, colors=c("orange","gray90")
#' )
#' }
#' @export
doppelganger <- function(data, variables=NULL, priority=c("centrality","peripherality","raw_order"), threshold=0.9) {
    # Find and select variables
    if(is.null(variables)) {
        variables <- colnames(data)
    }
    n_vars <- length(variables)
    data <- data[,variables]
    # Calculate the correlation matrix
    corrs <- cor(data, use="pairwise.complete.obs")
    diag(corrs) <- NA
    # Empty matrix to collect sizes of correlations
    sizes <- matrix(nrow=n_vars, ncol=n_vars, dimnames=list(variables, variables))
    # Reorder the variables according their centrality
    priority <- match.arg(priority)
    if(priority=="raw_order") {
        degree <- rep.int(NA, n_vars)
        names(degree) <- variables
    } else {
        # Calculate the size of each correlation
        check <- !is.na(data)
        couples <- combn(seq_len(ncol(check)), 2)
        weights <- sapply(
            seq_len(ncol(couples)),
            function(i) {
                sum(check[,couples[,i][1]] & check[,couples[,i][2]])
            }
        )
        sizes[lower.tri(sizes)] <- weights -> sizes[upper.tri(sizes)]
        # Calculate the representativeness of each variable
        degree <- colSums(abs(corrs)*sizes, na.rm=TRUE)/colSums(sizes, na.rm=TRUE)
        # Sort degree according to the direction of the algorithm
        direction <- priority=="centrality"
        degree <- sort(degree, decreasing=direction)
        # Reorder the vector of variables
        variables <- names(degree)
        # Reorder the matrices
        corrs <- corrs[variables, variables]
        sizes <- sizes[variables, variables]
    }
    # Prepare the output list
    out <- list(
        variables=variables,
        priority=priority,
        threshold=threshold,
        size_matrix=sizes,
        cor_matrix=corrs,
        cor_table=NULL,
        degree=degree,
        keep=variables,
        drop=character(0)
    )
    # Convert the correlation matrix in a tabular format
    corrs[lower.tri(corrs)] <- NA
    out$cor_table <- corrs <- corrs %>%
        data.frame(check.names=FALSE) %>%
        tibble::rownames_to_column("v1") %>%
        tidyr::pivot_longer(names_to="v2", values_to="cor", -.data$v1) %>%
        dplyr::filter(!is.na(.data$cor)) %>%
        dplyr::mutate(is_alias=abs(.data$cor)>=threshold)
    # Isolate the couples of aliases from the correlation table
    alias <-  corrs %>%
        dplyr::filter(!is.na(.data$cor) & .data$is_alias) %>%
        dplyr::select(.data$v1, .data$v2)
    # Iterate the alias table removing doppelgangers
    if(nrow(alias)>0) {
        corrs <- corrs %>%
            dplyr::filter(.data$is_alias) %>%
            dplyr::select(-.data$is_alias)
        drop <- NULL # variables to drop
        while(nrow(alias)>0 & length(variables)>0) {
            # Search aliases of the first listed variable
            var1_alias <- corrs %>%
                dplyr::filter(
                    (.data$v1==variables[1] | .data$v2==variables[1])
                ) %>%
                unlist(use.names=FALSE) %>%
                unique() %>%
                {variables[variables %in% .]}
            # If there is any alias for the first variable:
            if(length(var1_alias)>0) {
                # Isolate doppelgangers
                var_delete <- var1_alias[var1_alias != variables[1]]
                # Add doppelgangers to the drop list
                drop <- c(drop, var_delete)
                # Remove var1 and its aliases from the running list of variables
                variables <- variables[!(variables %in% var1_alias)]
                # Clear doppelgangers from the running alias table
                alias <- alias %>%
                    dplyr::filter(
                        !((.data$v1 %in% var_delete) | (.data$v2 %in% var_delete))
                    )
            } else {
                # No doppelgangers are found
                variables <- variables[-1]
            }
        }
        to_drop <- out$variables %in% drop
        out$keep <- out$variables[!to_drop]
        out$drop <- out$variables[to_drop]
    }
    class(out) <- c("doppelganger", "list")
    return(out)
}

#' Extract Correlations
#' @name inspect_corrs
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @description Given one or more variables, extracts the correlations with other
#' variables in the source data set.
#' @param dg_object An object of class \code{doppelganger}.
#' @param variable The variable of interest.
#' @export
inspect_corrs <- function(dg_object, variable) {
    dg_object$cor_table %>%
        dplyr::filter(.data$v1 %in% variable | .data$v2 %in% variable) %>%
        dplyr::mutate(
            v1 = ifelse(.data$v1==variable, "", .data$v1),
            v2 = ifelse(.data$v2==variable, "", .data$v2),
            variable = paste0(.data$v1, .data$v2)
        ) %>%
        dplyr::select(.data$variable, .data$is_alias, .data$cor) %>%
        dplyr::arrange(dplyr::desc(abs(.data$cor)))
}
