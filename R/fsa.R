#' @import dplyr
#' @import tidyr
NULL
#> NULL


#' Your title goes here: example: sensitivity point analysis
#'
#' Your description goes here: example: Creates a plot of the crayon colors in \code{\link{brocolors}}
#'
#' @param X an R array consisting all observed covariates
#' @param Y an R array
#' @return A number.
#' @examples
#' fsa_point()

fsa_point = function(X, Z, Y, epsilon1, epsilon0, truncpscore=c(0,1)){
  # pscore model
  pscore_logit = glm(Z ~ as.matrix(X), family='binomial')
  pscore = predict(pscore_logit, type='response')
  pscore = pmax(truncpscore[1], pmin(truncpscore[2], pscore))

  # outcome model
  mu1_lm = lm(Y ~ as.matrix(X), weights = Z)
  mu1 = predict(mu1_lm)
  mu0_lm = lm(Y ~ as.matrix(X), weights = 1-Z)
  mu0 = predict(mu0_lm)

  # 1. pred estimator
  pred_y1 = Z * Y + (1 - Z) * mu1 / epsilon1
  pred_y0 = Z * mu0 * epsilon0 + (1 - Z) * Y
  pred = mean(pred_y1 - pred_y0)

  # 2. ht estimator
  ht_y1 = (pscore * epsilon1 + 1 - pscore) * Z * Y / (pscore * epsilon1)
  ht_y0 = (pscore * epsilon0 + 1 - pscore) * (1 - Z) * Y / (1 - pscore)
  ht = mean(ht_y1 - ht_y0)

  # 3. hajek estimator
  haj_y1 = ht_y1 / mean((Z / pscore))
  haj_y0 = ht_y0 / mean(((1 - Z) / (1 - pscore)))
  haj = mean(haj_y1 - haj_y0)

  # 4. dr estimator
  dr_y1 = ht_y1 - (Z - pscore) * mu1 / (pscore * epsilon1)
  dr_y0 = ht_y0 - (pscore - Z) * mu0 * epsilon0 / (1 - pscore)
  dr = mean(dr_y1 - dr_y0)

  c(pred, ht, haj, dr)
}


#' Your title goes here: example: sensitivity boot
#'
#' Your description goes here: example: Creates a plot of the crayon colors in \code{\link{brocolors}}
#'
#' @param X an R array
#' @param Y an R array
#' @return A number.
#' @examples
#' fsa_boot()

fsa_boot = function(X, Z, Y, epsilon1, epsilon0, n_boot=500, truncpscore=c(0,1))
{
  fsa_result = fsa_point(X, Z, Y, epsilon1, epsilon0)
  n_boot = n_boot
  truncpscore = truncpscore

  ## nonparametric bootstrap
  n_sample   = length(Z)
  X          = as.matrix(X)
  boot_est   = replicate(n_boot,
                         {id_boot = sample(1:n_sample, n_sample, replace = TRUE)
                         fsa_point(X[id_boot, ], Z[id_boot], Y[id_boot], epsilon1, epsilon0, truncpscore)})
  # boot_point = apply(data.frame(boot_est), 1, mean)
  boot_se    = apply(data.frame(boot_est), 1, sd)

  res = cbind(fsa_result, boot_se)
  colnames(res) = c("est", "boot_se")
  rownames(res) = c("pred", "ht", "haj", "dr")
  res = data.frame(res)
  res$p_value = 2 * pnorm(-abs(res$est / res$boot_se))
  res$ci_lb = res$est + qnorm(0.025) * res$boot_se
  res$ci_ub = res$est + qnorm(0.975) * res$boot_se
  return(res)
}


#' Your title goes here: example: sensitivity boot
#'
#' Your description goes here: example: Creates a plot of the crayon colors in \code{\link{brocolors}}
#'
#' @param X an R array
#' @param Y an R array
#' @return A number.
#' @export
#' @examples
#' fsa_ate()

fsa_ate <- function(X, Z, Y, eps1_list, eps0_list, n_boot=500, truncpscore=c(0,1)) {
  pred <- data.frame()
  ht <- data.frame()
  haj <- data.frame()
  dr <- data.frame()
  for (epsilon1 in eps1_list){
    for (epsilon0 in eps0_list){
      res = fsa_boot(X, Z, Y, epsilon1, epsilon0)
      res$eps1 = epsilon1
      res$eps0 = epsilon0
      pred <- rbind(pred, res[1,])
      ht <- rbind(ht, res[2,])
      haj <- rbind(haj, res[3,])
      dr <- rbind(dr, res[4,])
    }
  }
  rownames(pred) <- NULL
  rownames(ht) <- NULL
  rownames(haj) <- NULL
  rownames(dr) <- NULL
  res_list <- list(pred = pred,
                   ht = ht,
                   haj = haj,
                   dr = dr)
  return(res_list)
}

#' Your title goes here: example: sensitivity boot
#'
#' Your description goes here: example: Creates a plot of the crayon colors in \code{\link{brocolors}}
#'
#' @param X an R array
#' @param Y an R array
#' @return A number.
#' @export
#' @examples
#' plot_contour()

plot_contour <- function(ate_result, caption='', estimator="dr", value="est"){
  if (caption == ''){
    caption <- "Contour plots of the estimated ATE with different values of sensitivity parameters"
  }
  eps1_list <- ate_result$eps1_list
  eps0_list <- ate_result$eps0_list
  res_wide = ate_result[[estimator]] %>%
    select(all_of(value), "eps1", "eps0") %>%
    spread(key="eps0", value=value)
  rownames(res_wide) <- eps1_list
  res_wide <- res_wide[,-1]
  z <- as.matrix(res_wide)
  contour(eps1_list, eps0_list, z, labcex = 1,
          plot.title = title(main = caption, xlab="epsilon1",  ylab="epsilon0"))
}


#' Your title goes here: example: sensitivity boot
#'
#' Your description goes here: example: Creates a plot of the crayon colors in \code{\link{brocolors}}
#'
#' @param X an R array
#' @param Y an R array
#' @return A number.
#' @examples
#' fsa_att_point()

fsa_att_point = function(X, Z, Y, epsilon0, truncpscore=c(0,1)){
  # pscore model
  pscore_logit = glm(Z ~ as.matrix(X), family='binomial')
  pscore = predict(pscore_logit, type='response')
  pscore = pmax(truncpscore[1], pmin(truncpscore[2], pscore))
  odd = pscore / (1 - pscore)

  # outcome model
  mu1_lm = lm(Y ~ as.matrix(X), weights = Z)
  mu1 = predict(mu1_lm)
  mu0_lm = lm(Y ~ as.matrix(X), weights = 1-Z)
  mu0 = predict(mu0_lm)

  n1 = sum(Z)
  mu_t1 = sum(Z * Y) / n1

  # 1. pred estimator
  pred_t0 = sum(Z * epsilon0 * mu0) / n1
  pred = mu_t1 - pred_t0

  # 2. ht estimator
  ht_t0 = sum(epsilon0 * odd * (1 - Z) * Y) / n1
  ht = mu_t1 - ht_t0

  # 3. hajek estimator
  haj_t0 = sum(epsilon0 * odd * (1 - Z) * Y) / sum(odd * (1 - Z))
  haj = mu_t1 - haj_t0

  # 4. dr estimator
  dr_t0 = sum(Z * epsilon0 * mu0 + epsilon0 * odd * (1 - Z) * (Y - mu0)) / n1
  dr = mu_t1 - dr_t0

  c(pred, ht, haj, dr)
}


#' Your title goes here: example: sensitivity boot
#'
#' Your description goes here: example: Creates a plot of the crayon colors in \code{\link{brocolors}}
#'
#' @param X an R array
#' @param Y an R array
#' @return A number.
#' @examples
#' fsa_att_boot()

fsa_att_boot = function(X, Z, Y, epsilon0, n_boot=500, truncpscore=c(0,1))
{
  fsa_result = fsa_att_point(X, Z, Y, epsilon0)
  n_boot = n_boot
  truncpscore = truncpscore

  ## nonparametric bootstrap
  n_sample   = length(Z)
  X          = as.matrix(X)
  boot_est   = replicate(n_boot,
                         {id_boot = sample(1:n_sample, n_sample, replace = TRUE)
                         fsa_att_point(X[id_boot, ], Z[id_boot], Y[id_boot], epsilon0, truncpscore)})
  # boot_point = apply(data.frame(boot_est), 1, mean)
  boot_se    = apply(data.frame(boot_est), 1, sd)

  res = cbind(fsa_result, boot_se)
  colnames(res) = c("est", "boot_se")
  rownames(res) = c("pred", "ht", "haj", "dr")
  res = data.frame(res)
  res$p_value = 2 * pnorm(-abs(res$est / res$boot_se))
  res$ci_lb = res$est + qnorm(0.025) * res$boot_se
  res$ci_ub = res$est + qnorm(0.975) * res$boot_se
  return(res)
}


#' Your title goes here: example: sensitivity boot
#'
#' Your description goes here: example: Creates a plot of the crayon colors in \code{\link{brocolors}}
#'
#' @param X an R array
#' @param Y an R array
#' @return A number.
#' @export
#' @examples
#' fsa_att()

fsa_att <- function(X, Z, Y, eps0_list, n_boot=500, truncpscore=c(0,1)) {
  pred <- data.frame()
  ht <- data.frame()
  haj <- data.frame()
  dr <- data.frame()
  for (epsilon0 in eps0_list){
    res = fsa_att_boot(X, Z, Y, epsilon0)
    res$eps0 = epsilon0
    pred <- rbind(pred, res[1,])
    ht <- rbind(ht, res[2,])
    haj <- rbind(haj, res[3,])
    dr <- rbind(dr, res[4,])
  }
  rownames(pred) <- NULL
  rownames(ht) <- NULL
  rownames(haj) <- NULL
  rownames(dr) <- NULL
  res_list <- list(pred = pred,
                   ht = ht,
                   haj = haj,
                   dr = dr)
  return(res_list)
}

