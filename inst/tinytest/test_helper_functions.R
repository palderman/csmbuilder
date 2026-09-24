library(tinytest)

#############################
# Modified Arrhenius function
#############################

prms <- c(ko = 1,
          H = 317470.,
          E = 91450.,
          To = 25)

k <- csmbuilder::csm_mod_arr(Tt = -273:100,
                             ko = prms["ko"],
                             H = prms["H"],
                             E = prms["E"],
                             To = prms["To"])

expect_true(all(k >= 0))
expect_true(all(k <= 1))

expect_equal(csmbuilder::csm_mod_arr(Tt = prms["To"],
                                    ko = prms["ko"],
                                    H = prms["H"],
                                    E = prms["E"],
                                    To = prms["To"]),
                                    1.)

expect_equivalent(
  csmbuilder::csm_mod_arr(Tt = c(-273, 273),
                          ko = prms["ko"],
                          H = prms["H"],
                          E = prms["E"],
                          To = prms["To"]),
             c(0, 0))

###########################
# Arrhenius fraction active
###########################

f_a <- csmbuilder::csm_arr_fr_active(Tt = -273:100,
                                     H = prms["H"],
                                     E = prms["E"],
                                     To = prms["To"])

expect_true(all(f_a >= 0))
expect_true(all(f_a <= 1))

Thalf <- 1/(1/(prms["To"]+273.15)-log(1+(2*(prms["H"]-prms["E"])-prms["H"])/prms["E"])*8.314/prms["H"]) - 273.15

expect_equal(
  csmbuilder::csm_arr_fr_active(Tt = Thalf,
                                H = prms["H"],
                                E = prms["E"],
                                To = prms["To"]),
  1/2)


expect_equivalent(
  csmbuilder::csm_arr_fr_active(Tt = c(-273, 273),
                                H = prms["H"],
                                E = prms["E"],
                                To = prms["To"]),
  c(1., 0.))

###################
# csm_logistic()
###################

x <- seq(-10, 10, by = 0.01)

r <- -log(1/0.01 - 1)/(0 - 0.5)

k <- csmbuilder::csm_logistic(x, r = r, x0 = 0.5)

expect_true(all(k >= 0))
expect_true(all(k <= 1))
expect_true(all(k[x<0.5] < 0.5))
expect_true(all(k[x>0.5] > 0.5))

expect_equal(
  csmbuilder::csm_logistic(0.5, r = r, x0 = 0.5),
  0.5
)

expect_equal(
  csmbuilder::csm_logistic(c(0, 1), r = r, x0 = 0.5),
  c(0.01, 0.99)
)

r <- -log(1/0.99 - 1)/(0 - 0.5)

k <- csmbuilder::csm_logistic(x, r = r, x0 = 0.5)

expect_true(all(k >= 0))
expect_true(all(k <= 1))
expect_true(all(k[x<0.5] > 0.5))
expect_true(all(k[x>0.5] < 0.5))

expect_equal(
  csmbuilder::csm_logistic(0.5, r = r, x0 = 0.5),
  0.5
)

expect_equal(
  csmbuilder::csm_logistic(c(0, 1), r = r, x0 = 0.5),
  c(0.99, 0.01)
)

#######################
# csm_calc_logistic_r()
#######################

x <- seq(-10, 10, by = 0.01)

expect_true(
  csmbuilder::csm_calc_logistic_r(0, 0.01, 0.5) > 0
)

expect_true(
  csmbuilder::csm_calc_logistic_r(1, 0.01, 0.5) < 0
)

r <- csmbuilder::csm_calc_logistic_r(0, 0.01, 0.5)

k <- csmbuilder::csm_logistic(x, r = r, x0 = 0.5)

expect_true(all(k >= 0))
expect_true(all(k <= 1))
expect_true(all(k[x<0.5] < 0.5))
expect_true(all(k[x>0.5] > 0.5))

expect_equal(
  csmbuilder::csm_logistic(c(0, 1), r = r, x0 = 0.5),
  c(0.01, 0.99)
)

r <- csmbuilder::csm_calc_logistic_r(0, 0.99, 0.5)

k <- csmbuilder::csm_logistic(x, r = r, x0 = 0.5)

expect_true(all(k >= 0))
expect_true(all(k <= 1))
expect_true(all(k[x<0.5] > 0.5))
expect_true(all(k[x>0.5] < 0.5))

expect_equal(
  csmbuilder::csm_logistic(c(0, 1), r = r, x0 = 0.5),
  c(0.99, 0.01)
)


#######################
# csm_smooth_bounds()
#######################

x <- seq(0, 1, by = 0.01)

min_val <- 0.3
max_val <- 0.7

k <- csmbuilder::csm_smooth_bounds(x, min_val, max_val)

expect_true(all(k >= min_val))
expect_true(all(k <= max_val))
expect_true(all(k[x<(max_val+min_val)/2] < (max_val+min_val)/2))
expect_true(all(k[x>(max_val+min_val)/2] > (max_val+min_val)/2))

expect_equal(
  csmbuilder::csm_smooth_bounds((max_val+min_val)/2, min_val, max_val),
  (max_val+min_val)/2
)

###################
# csm_log_sum_exp()
###################

csmbuilder::csm_log_sum_exp(1e-16, 1e16)

