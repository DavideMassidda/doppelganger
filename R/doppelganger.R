#' Find Doppelgangers
#' @name doppelganger
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @description Analyze correlations to select a subset of non-redundant variables.
#' @param data Data table source.
#' @param variables Variables included in \code{data} to consider.
#' @param priority Ranking method to prioritize variables. When "Centrality" is
#' selected, variables with higher mean absolute correlation will have higher
#' priority. "Raw order" prioritize variables using their position within the
#' vector \code{variables}.
#' @param threshold Correlation cut-off (absolute value) to identify doppelgangers.
#' Variables with absolute correlation values equal or greater to this value will be
#' considered aliases of each other.
#' @details The variable priority is established by calculating a centrality index
#' for each variable, and sorting then accordingly.
#' @examples
#' \dontrun{
#' data(measures)
#' dg <- doppelganger(measures, threshold=0.7)
#' qgraph::qgraph(cor(measures), minimum=0.65, edge.labels=TRUE, layout="spring")
#' }
#' @export
doppelganger <- function(data, variables=NULL, priority=c("centrality","raw_order"), threshold=0.9) {
    # Find variables (if not provided)
    if(is.null(variables)) {
        variables <- colnames(data)
    }
    # Calculate the correlation matrix
    corrs <- suppressWarnings(
        data %>%
            dplyr::select(tidyselect::all_of(variables)) %>%
            cor(use="pairwise.complete.obs") %>%
            `diag<-`(NA)
    )
    # Reorder the variables according their importance
    n_vars <- length(variables)
    importance <- rep.int(NA, n_vars)
    names(importance) <- variables
    priority <- match.arg(priority)
    if(priority!="order") {
        weights <- abs(corrs)
        #if(priority=="separation") {
        #    to_correct <- !is.na(weights) & weights<threshold
        #    if(any(to_correct)) {
        #        weights[to_correct] <- 1-weights[to_correct]
        #    }
        #}
        importance <- sort(colSums(weights, na.rm=TRUE)/n_vars, decreasing=TRUE)
        variables <- names(importance)
    }
    # Reorder the correlation matrix according to the variable importance
    corrs <- corrs[variables, variables]
    # Convert the correlation matrix in a tabular format
    corrs[lower.tri(corrs)] <- NA
    corrs <- corrs %>%
        data.frame(check.names=FALSE) %>%
        tibble::rownames_to_column("v1") %>%
        tidyr::pivot_longer(names_to="v2", values_to="cor", -.data$v1) %>%
        dplyr::filter(!is.na(.data$cor)) %>%
        dplyr::mutate(is_alias=abs(.data$cor)>=threshold)
    # Isolate the couples of aliases from the correlation table
    alias <-  corrs %>%
        dplyr::filter(!is.na(.data$cor) & .data$is_alias) %>%
        dplyr::select(.data$v1, .data$v2)
    # Prepare the output list
    out <- list(
        variables=variables,
        priority=priority,
        importance=importance,
        threshold=threshold,
        corrs=corrs,
        keep=variables,
        drop=character(0)
    )
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
#' @description Given a variable, extracts the correlations with other variables in a data set.
#' @param dg_object An object of class \code{doppelganger}.
#' @param variable The variable of interest.
#' @export
inspect_corrs <- function(dg_object, variable) {
    dg_object$corrs %>%
        dplyr::filter(.data$v1==variable | .data$v2==variable) %>%
        dplyr::mutate(
            v1 = ifelse(.data$v1==variable, "", .data$v1),
            v2 = ifelse(.data$v2==variable, "", .data$v2),
            variable = paste0(.data$v1, .data$v2)
        ) %>%
        dplyr::select(.data$variable, .data$is_alias, .data$cor) %>%
        dplyr::arrange(dplyr::desc(abs(.data$cor)))
}

#' Simulated data
#' @name measures
#' Simulated data matrix to show hoe the package \code{doppelganger} works.
#' @docType data
#' @keywords data
#' @export
