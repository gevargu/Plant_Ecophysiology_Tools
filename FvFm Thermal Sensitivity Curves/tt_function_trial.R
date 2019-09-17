#'
#' # Thermal tolerance curve analysis package trials
#'
#'

trial <- read.csv(file = "Thermal_tolerance_sample.csv")
trial

#'
#' ## Load the function
#'

source(file = "scripts/functions/T_tol_func_V2.R")

#'
#' ## Fit the thermal tolerance curves and obtain parameters:
#'
#' ### Not saving the graphs as jpeg outside the R environment
#' 
#' 
Thermal_output <- Thermal_tol(fv_fm = fv.fm_aft,
                              temperature = tmp,
                              species = code,
                              sp_levels = c("BURSIM","BROALI"),
                              individual = ind,
                              data = trial,
                              print_graph = FALSE,
                              graph_dir = "graphs/thermal_sp_id/")

Thermal_output

#' 
#' ### Condense data together and visualize values per species
#' 

Thermal_resp <- do.call(what = "rbind",Thermal_output)

plot(Thermal_resp$ctmax~Thermal_resp$species,main="C50")
plot(Thermal_resp$tcrit~Thermal_resp$species,main="Tcrit")
