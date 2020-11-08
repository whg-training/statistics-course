linear_regression.ll <- function(
    outcome,
    params = list(
        beta,
        sigma2
    ),
    X
) {
    predictor = X %*% params$beta
    gaussian.ll( outcome, params = list( mean = predictor, sigma2 = params$sigma2 ))
}
