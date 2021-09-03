library(igraph)
library(purrr)
library(brainGraph)
library(progress)
library(DirectedClustering)
library(DescTools)
library(aricode)
library(rslurm)


run_wrapper <- function(seed, alpha, nNodes){
  
  source("/sfs/qumulo/qhome/ycp6wm/Thresholding_Sims/bin/rand_corrmat.R")
  source("/sfs/qumulo/qhome/ycp6wm/Thresholding_Sims/bin/helper_functions.R")
  set.seed(seed)
  nets <- gen_corr_mats(5, nNodes,
                        function(d){
                          return(rep(alpha, d-1))})
  
  out_ge <-  apply_GE(nets)
  out_mod <- apply_modularity(nets)
  out_swp  <- apply_SWP(nets)
  
  out = merge(out_ge, out_mod, by = c("ID", "thres", "thres_type"))
  out = merge(out, out_swp, by = c("ID", "thres", "thres_type"))
  out$seed = seed
  out$alpha = alpha
  out$nNodes = nNodes
  return(out)
}



design_frame = data.frame(seed = 1:500, alpha = 1, nNodes = 100)

sjob <- slurm_apply(run_wrapper, design_frame, jobname = "100n_run", nodes=nrow(design_frame), cpus_per_node = 1, submit = TRUE,
                    slurm_options = c(account = "netlab", partition = "standard",time = "5:00:00"), preschedule_cores = F)

save(sjob, file = "100n_run_sjob.Rdata")


design_frame = data.frame(seed = 1:500, alpha = 1, nNodes = 200)

sjob <- slurm_apply(run_wrapper, design_frame, jobname = "200n_run", nodes=nrow(design_frame), cpus_per_node = 1, submit = TRUE,
                    slurm_options = c(account = "netlab", partition = "standard",time = "5:00:00"), preschedule_cores = F)

save(sjob, file = "200n_run_sjob.Rdata")
