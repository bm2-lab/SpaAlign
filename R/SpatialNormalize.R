SpatialNormalize <-
function (expression_profile, ST_method = c("osmFISH", "merFISH", "STARmap", "Others")) {
    ST_method <- ST_method[1]
    if (ST_method %in% c("osmFISH","merFISH", "STARmap")) {
        normalize_scale <- function(x) {
            return(x/sum(x))
        }
        expression_profile[is.na(expression_profile)] <- 0
        expression_profile <- expression_profile[apply(abs(expression_profile), 
            1, sum) > 0, ]
        expression_profile <- apply(expression_profile, 2, normalize_scale)
        expression_profile[is.na(expression_profile)] <- 0
        expression_profile[is.nan(expression_profile)] <- 0
        expression_profile <- log(expression_profile + 1)
    }
    else {
        normalize_scale_factor <- function(x) {
            return(10000 * (x/sum(x)))
        }
        expression_profile[is.na(expression_profile)] <- 0
        expression_profile <- expression_profile[apply(abs(expression_profile), 
            1, sum) > 0, ]
        expression_profile <- apply(expression_profile, 2, normalize_scale_factor)
        expression_profile[is.na(expression_profile)] <- 0
        expression_profile[is.nan(expression_profile)] <- 0
        expression_profile <- log(expression_profile + 1)
        return(expression_profile)
    }
}
