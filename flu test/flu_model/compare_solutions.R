# Load results
load("MOL_results.RData")       
S_mol <- Sm; V_mol <- Vm; E_mol <- Em; I_mol <- Im; R_mol <- Rm

load("Explicit_results.RData")  
S_exp <- Se; V_exp <- Ve; E_exp <- Ee; I_exp <- Ie; R_exp <- Re

load("FVM_results.RData")       
S_FVM <- Sf; V_FVM <- Vf; E_FVM <- Ef; I_FVM <- If; R_FVM <- Rf

load("Crank-Nicolson_results.RData")       
S_CN <- Sc; V_CN <- Vc; E_CN <- Ec; I_CN <- Ic; R_CN <- Rc


load("results_strang.RData")  
S_str <- Ss; V_str <- Vs; E_str <- Es; I_str <- Is; R_str <- Rs

load("Spectra_results.RData")  
S_sp <- Sp; V_sp <- Vp; E_sp <- Ep; I_sp <- Ip; R_sp <- Rp


compare_solutions <- function(S1, V1, E1, I1, R1, S2, V2, E2, I2, R2, show_summary = TRUE) {
  calc_error <- function(u1, u2) {
     abs(u1 - u2)  
  }
  
  err_S <- calc_error(S1, S2)
  err_V <- calc_error(V1, V2)
  err_E <- calc_error(E1, E2)
  err_I <- calc_error(I1, I2)
  err_R <- calc_error(R1, R2)
  
  if (show_summary) {
    mean_errors <- c(
      S = mean(err_S, na.rm = TRUE),
      V = mean(err_V, na.rm = TRUE),
      E = mean(err_E, na.rm = TRUE),
      I = mean(err_I, na.rm = TRUE),
      R = mean(err_R, na.rm = TRUE)
    )
    
    range_errors <- rbind(
      S = range(err_S, na.rm = TRUE),
      V = range(err_V, na.rm = TRUE),
      E = range(err_E, na.rm = TRUE),
      I = range(err_I, na.rm = TRUE),
      R = range(err_R, na.rm = TRUE)
    )
    
    cat("Mean Errors between methods:\n")
    print(round(mean_errors, 4))
    
    cat("\nRange of Errors (min, max):\n")
    print(round(range_errors, 6))
  }
  
  return(list(
    S_error_matrix = err_S,
    V_error_matrix = err_V,
    E_error_matrix = err_E,
    I_error_matrix = err_I,
    R_error_matrix = err_R
  ))
}
###############################################################################
# IMPORTANT NOTE ON SOLUTION COMPARISONS:
#
# All comparison blocks below assume that each method (exp, RK4, CN, str) was
# run with THE SAME INITIAL CONDITIONS (IP). For valid comparisons between methods:
#
# 1. ALL methods must use either IP=1 OR IP=2 consistently
# 2. NEVER mix different IP values across methods
# 3. The reference solution (S_mol, V_mol, etc.) must use the same IP as all others
#
# Failure to use consistent initial conditions will make comparisons invalid
# as the methods would be solving different problems with different starting points.
###############################################################################

errors <- compare_solutions(S_mol, V_mol, E_mol, I_mol, R_mol,
                            S_exp, V_exp, E_exp, I_exp, R_exp)
print(errors)

errors <- compare_solutions(S_mol, V_mol, E_mol, I_mol, R_mol,
                            S_FVM, V_FVM, E_FVM, I_FVM, R_FVM)
print(errors)

errors <- compare_solutions(S_mol, V_mol, E_mol, I_mol, R_mol,
                            S_CN, V_CN, E_CN, I_CN, R_CN)
print(errors)

errors <- compare_solutions(S_mol, V_mol, E_mol, I_mol, R_mol,
                            S_str, V_str, E_str, I_str, R_str)
print(errors)

errors <- compare_solutions(S_mol, V_mol, E_mol, I_mol, R_mol,
                            S_sp, V_sp, E_sp, I_sp, R_sp)
print(errors)