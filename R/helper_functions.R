#' Linearly interpolate a time-varying variable at specific time point
#'
#' @export
#'
#' @param x a vector of a time-varying values for which to interpolate
#' @param t_ind a vector of times corresponding to the values in x
#' @param t a single time point at which to return a value
#' @param method a string indicating what method to use for
#'  interpolation. One of: "linear"
#' @param search a string indicating what method to use for
#'  finding the correct indices within t_ind. One of: "bisection",
#'  "interpolation"
#'
#' @returns
#' the numeric value of the time-varying variable interpolated at time t
#'
csm_get_at_t <- function(x, t_ind, t,
                         method = "linear",
                         search = c("interpolation", "bisection",
                                    "bruteforce", "i=t+1")){
  search <- search[1]

  if(search == "i=t+1"){
    return(x[t+1])
  }

  method <- method[1]

  t_ind_len <- length(t_ind)
  if(length(t) != 1) stop("Length of argument t not equal to 1.")
  method <- method[1]
  if(method != "linear") stop("method must be one of: linear.")
  search <- search[1]
  if(!any(search != c("interpolation", "bisection", "bruteforce",
                      "i=t+1"))){
    stop("method must be one of: linear.")
  }
  if(t > t_ind[t_ind_len]){
    stop("t is greater than last value of t_ind.")
  }
  if(t < t_ind[1]) stop("t is less than last value of t_ind.")

  if(search == "i=t+1"){
    return(x[t+1])
  }

  i <- floor(t) |> max(1)
  j <- ceiling(t) |> min(length(t_ind)) |> max(1)

  if(t_ind[i] == t){
    return(x[i])
  }else{
    if(search == "bisection"){
      if(t_ind[i] > t) i <- 1
      if(t_ind[j] < t) j <- length(t_ind)
      if(t_ind[i] == t) return(x[i])
      if(t_ind[j] == t) return(x[j])
      while(i < j && (i + 1) != j){
        ij_cand <- floor((i+j)/2)
        if(t_ind[ij_cand] == t){
          return(x[ij_cand])
        }else if(t_ind[ij_cand] < t){
          i <- ij_cand
        }else{
          j <- ij_cand
        }
      }
    }else if(search == "interpolation"){
      i <- 1
      j <- length(t_ind)
      while(i < j && (i + 1) != j){
        if(t_ind[i] > t) i <- 1
        if(t_ind[j] < t) j <- length(t_ind)
        m <- (j - i)/(t_ind[j] - t_ind[i])
        b <- j - m*t_ind[j]
        ij <- m*t+b
        i <- floor(ij)
        if(t_ind[i] == t) return(x[i])
        j <- ceiling(ij)
      }
    }else if(search == "bruteforce"){
      i <- which.min(abs(t_ind - t))
      if(t_ind[i] == t){
        return(x[i])
      }else if(t_ind[i] > t){
        j = i
        i = j - 1
      }else{
        j = i + 1
      }
    }
  }

  # i = which.min(abs(t_ind - t))
  # if(t_ind[i] == t){
  #   j = i
  # }else if(t_ind[i] > t){
  #   j = i
  #   i = j - 1
  # }else{
  #   j = i + 1
  # }

  if(!(i > 0)) stop("i must be greater than 0.")
  if(!(j <= length(x))) stop("j must be less than or equal to length(x).")

  if(t_ind[i] == t){
    return(x[i])
  }else{
    return(x[i]*(t_ind[j]-t) + x[j]*(t-t_ind[i]))/(t_ind[j] - t_ind[i])
  }

}

#' Modified Arrhenius function
#'
#' @export
#'
#' @param Tt temperature in Celsius
#' @param ko reaction rate at the optimum temperature (To)
#' @param H deactivation energy parameter
#' @param E activation energy parameter
#' @param To optimum temperature in Celsius
#'
#' @returns
#' a numeric value of the reaction rate at temperature Tt
#'
csm_mod_arr <- function(Tt, ko, H, E, To){
  R <- 8.314
  unname(ko*(H*exp(E/R*(1/(To+273.15)-1/(Tt + 273.15))))/(H - E*(1-exp(H/R*(1/(To+273.15)-1/(Tt + 273.15))))))
}

#' Fraction of active enzymes based on modified Arrhenius function
#'
#' This function computes the fraction of active enzymes according to the
#' the modified Arrhenius function. The fraction of denatured enzymes can be
#' calculated by subtracting this function from 1.
#'
#' @export
#'
#' @param Tt temperature in Celsius
#' @param H deactivation energy parameter
#' @param E activation energy parameter
#' @param To optimum temperature in Celsius
#'
#' @returns
#' a numeric value of the fraction of active enzymes rate at temperature Tt
#'
csm_arr_fr_active <- function(Tt, H, E, To){
  R <- 8.314
  unname((H-E)/(H - E*(1-exp(H/R*(1/(To+273.15)-1/(Tt + 273.15))))))
}

#' Hill equation for up-regulation
#'
#' @export
#'
#' @param L a numeric value providing the ligand concentration
#' @param K a numeric value providing the ligand concentration at half occupation
#' @param n a numeric value providing the Hill coefficient
#'
#' @returns
#' a numeric value between 0 and 1
#'
csm_hill_up_reg <- function(L, K, n){
  L_n <- L**n
  L_n/(L_n + K**n)
}

#' Hill equation for down-regulation
#'
#' @export
#'
#' @param L a numeric value providing the ligand concentration
#' @param K a numeric value providing the ligand concentration at half occupation
#' @param n a numeric value providing the Hill coefficient
#'
#' @returns
#' a numeric value between 0 and 1
#'
csm_hill_down_reg <- function(L, K, n){
  K_n <- K**n
  K_n/(K_n + L**n)
}

#' Hill coefficient for up-regulation
#'
#' @export
#'
#' @param x a numeric value providing the value of the Hill equation at L
#' @param L a numeric value providing the ligand concentration
#' @param K a numeric value providing the ligand concentration at half occupation
#'
#' @returns
#' a numeric value
#'
csm_calc_n_up <- function(x, L, K){
  log((1-x)/x)/log(K/L)
}

#' Hill coefficient for down-regulation
#'
#' @export
#'
#' @param x a numeric value providing the value of the Hill equation at L
#' @param L a numeric value providing the ligand concentration
#' @param K a numeric value providing the ligand concentration at half occupation
#'
#' @returns
#' a numeric value
#'
csm_calc_n_down <- function(x, L, K){
  log((1-x)/x)/log((L/K))
}

#' Logistic switch function
#'
#' @export
#'
#' @param x a numeric value providing the input variable
#' @param r a numeric value providing the rate parameter of the logistic function
#' @param x0 a numeric value providing the mid-point parameter of the logistic function
#'
#' @returns
#' a numeric value between 0 and 1
#'
csm_logistic <- function(x, r, x0){
  1./(1.+exp(-r*(x - x0)))
}

#' Calculate logistic rate parameter
#'
#' @export
#'
#' @param x a numeric value providing a target input value
#' @param f_x a numeric value providing the target logistic function value at x
#' @param x0 a numeric value providing the mid-point parameter of the logistic function
#'
#' @returns
#' a numeric value
#'
csm_calc_logistic_r <- function(x, f_x, x0){
  -log(1/f_x - 1)/(x - x0)
}


#' Calculate smooth approximation to a bounded variable
#'
#' @export
#'
#' @param x a numeric input value
#' @param min_bound a numeric value of the minimum bound
#' @param max_bound a numeric value of the maximum bound
#'
#' @returns
#' a numeric value
#'
csm_smooth_bounds <- function(x, min_val, max_val){
  x0 <- (max_val + min_val)/2
  r <- csm_calc_logistic_r(x = min_val, f_x = 0.01, x0 = x0)
  min_val + (max_val - min_val)*csm_logistic(x, r, x0)
}

#' Calculate the LogSumExp smooth approximation to a bounded variable
#'
#' @export
#'
#' @param x1 a numeric input value
#' @param x2 a numeric input value
#' @param alpha a numeric value providing the sharpness of the approximation function
#'
#' @returns
#' a numeric value
#'
csm_log_sum_exp <- function(x, y, alpha = 1){
  log(sum(exp(alpha*x)+exp(alpha*y)))/alpha
}

#' Calculate the LogSumExp smooth approximation to a bounded variable
#'
#' @export
#'
#' @param x1 a numeric input value
#' @param x2 a numeric input value
#' @param alpha a numeric value providing the sharpness of the approximation function
#'
#' @returns
#' a numeric value
#'
csm_log_sum_exp <- function(x, y, alpha = 1){
  stopifnot(length(x) == 1)
  stopifnot(length(y) == 1)
  stopifnot(length(alpha) == 1)
  max_val <- max(x,y)*alpha
  (max_val + log((exp(alpha*x-max_val)+exp(alpha*y-max_val))))/alpha
}


#' Calculate the Boltzmann operator smooth approximation to the maximum function
#'
#' @export
#'
#' @param x1 a scalar numeric input value
#' @param x2 a scalar numeric input value
#' @param alpha a scalar numeric value providing the sharpness of the approximation function
#'
#' @returns
#' a numeric value
#'
csm_boltz_max <- function(x1, x2, alpha = 1){
  stopifnot(length(x1) == 1)
  stopifnot(length(x2) == 1)
  stopifnot(length(alpha) == 1)
  ex1 <- exp(alpha*x1)
  ex2 <- exp(alpha*x2)
  (x1*ex1+x2*ex2)/(ex1+ex2)
}


#' @export
csm_hs_max <- function(x1, x2, epsilon = 1e-6){
(x1 + x2 + sqrt((x1-x2)**2 + epsilon))/2
}

#' @export
csm_hs_min <- function(x1, x2, epsilon = 1e-6){
  (x1 + x2 - sqrt((x1-x2)**2 + epsilon))/2
}
