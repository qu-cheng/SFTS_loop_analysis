# ======================================================
# Scenario grid: 3 L:N:A ratios x 3 kappa x 3 tick abundance = 27
#
#   L:N:A      100:10:1 / 200:10:1 / 52:5:1
#   kappa      0.15 / 0.05 / 0.5
#   abundance  low / medium / high
# =====================================================

ncore = parallel::detectCores() - 2
cl = makeCluster(ncore)
registerDoParallel(cl)

print(unique(all.dat$Var))


#================= scenario grid ==================

# larval : nymphal : adult weights, keyed by the tag used in the file name

lna.list = list(
  "100" = c(L = 100, N = 10, A = 1),
  "200" = c(L = 200, N = 10, A = 1),
  "52"  = c(L = 52,  N = 5,  A = 1))

kappa.list = c(0.15, 0.05, 0.5)

sce.list = c("low", "medium", "high")

n.draw  = 1000   # draws per location-year
var.idx = 1:38   # 38 location-years

skip.existing = TRUE



worker.objs = c("all.dat", "tick.scaled", "base.para",
                "compute_R0", "compute_K_discrete",
                "loopanalysis.Adams", "loopanalysis.Guneralp",
                "decomp.new", "loop.random.new", "calculate_elasticity",
                "FindCycles", "chain.to.df")


run_scenario = function(lna, kappa, sce, n.j = n.draw, var.idx = 1:38) {
  
  out.file = str_c("./", lna, "_", kappa, "_", sce, "_sample1000.rds")
  
  if (skip.existing && file.exists(out.file)) {
    message("skip (already there): ", out.file)
    return(invisible(out.file))
  }
  
  alpha_L = alpha_N = alpha_A = kappa
  
  w = lna.list[[lna]] / sum(lna.list[[lna]])
  p.L = unname(w["L"])
  p.N = unname(w["N"])
  p.A = unname(w["A"])
  
  sce.res = vector("list", length(var.idx))
  
  for (k in seq_along(var.idx)) {
    
    var_k = unique(all.dat$Var)[var.idx[k]]
    
    message("Var ", var.idx[k], " | L:N:A ", lna, " | kappa ", kappa, " | ", sce)
    
    res = foreach(
      j = 1:n.j,
      .packages = c("tidyverse", "igraph", "CVXR"),
      .export = worker.objs,
      .errorhandling = "pass"
    ) %dopar% {
      
      tryCatch({
        
        abund = tick.scaled[j, sce]
        
        current.dat = all.dat %>%
          filter(Var == var_k) %>%
          mutate(total = sum(Nymph + Larva + Adult),
                 sum_N = sum(Nymph),
                 sum_L = sum(Larva),
                 sum_A = sum(Adult),
                 wN = p.N,
                 wL = p.L,
                 wA = p.A,
                 target_N = total * wN,
                 target_L = total * wL,
                 target_A = total * wA,
                 kN = target_N / sum_N,
                 kL = target_L / sum_L,
                 kA = target_A / sum_A,
                 Nymph = Nymph * kN,
                 Larva = Larva * kL,
                 Adult = Adult * kA) %>%
          select(1:5) %>%
          mutate(no.tick = sum(Nymph + Larva + Adult, na.rm = TRUE),
                 scale.para = no.tick / (abund * 365)) %>%
          mutate(across(c(Adult, Nymph, Larva), ~ .x / scale.para), scale.tick = sum(Nymph + Larva + Adult, na.rm = TRUE)) %>%
          mutate(a_N = Nymph / sum(Nymph),
                 a_L = Larva / sum(Larva),
                 a_A = Adult / sum(Adult)) %>%
          mutate(Pr_L0 = (alpha_L / (alpha_L + Larva))^alpha_L,
                 C_LL = sum(ifelse(Pr_L0 == 1, 0,
                                   a_L * ((Larva / (1 - Pr_L0)) - 1)))) %>%
          mutate(Pr_N0 = (alpha_N / (alpha_N + Nymph))^alpha_N,
                 C_NN = sum(ifelse(Pr_N0 == 1, 0,
                                   a_N * ((Nymph / (1 - Pr_N0)) - 1)))) %>%
          mutate(Pr_A0 = (alpha_A / (alpha_A + Adult))^alpha_A,
                 C_AA = sum(ifelse(Pr_A0 == 1, 0,
                                   a_A * ((Adult / (1 - Pr_A0)) - 1)))) %>%
          mutate(C_LN = sum(a_N * Larva),
                 C_LA = sum(a_A * Larva),
                 C_NL = sum(a_L * Nymph),
                 C_NA = sum(a_A * Nymph),
                 C_AL = sum(a_L * Adult),
                 C_AN = sum(a_N * Adult))
        
        Res.origin = compute_R0(current.dat)
        R0 = Res.origin[[1]]
        K.mat = Res.origin[[2]]
        
        
        E.mat = matrix(0, 7, 7)
        
        for (ei in 1:7)
          for (ej in 1:7) {
            delta.kij = K.mat[ei, ej] * 0.001
            if (delta.kij == 0) next
            
            K.mat.pert = K.mat
            K.mat.pert[ei, ej] = K.mat[ei, ej] + delta.kij
            
            delta.R0 = Re(eigen(K.mat.pert)$values[1]) - R0
            E.mat[ei, ej] = K.mat[ei, ej] / R0 * delta.R0 / delta.kij
          }
        
        E.mat[is.na(E.mat)] = 0
        gemat = E.mat
        
        result_adams = as_tibble(loopanalysis.Adams(gemat))
        result_Guneralp = as_tibble(loopanalysis.Guneralp(gemat))
        result_package = decomp.new(gemat)
        
        loop_info = list()
        
        for (i in seq_along(result_package)) {
          tmp = calculate_elasticity(result_package[[i]])
          if (!is.null(tmp)) {
            loop_info[[i]] = data.frame(
              loop = gsub(" ", "", tmp$path),
              length = tmp$length,
              loop.package.elas = tmp$elasticity
            )
          }
        }
        
        current.loop.res = bind_rows(loop_info) %>%
          mutate(
            R0.origin = R0,
            sce = sce,
            j = j,
            Var = var_k
          ) %>%
          left_join(result_Guneralp, by = "loop") %>%
          left_join(result_adams[, c(1, 2, 5, 6)])
        
        current.Kmat = K.mat %>%
          as.numeric() %>%
          matrix(nrow = 1) %>%
          as.data.frame() %>%
          set_names(paste0("K", rep(1:7, times = 7), rep(1:7, each = 7))) %>%
          mutate(R0.origin = R0, sce = sce, j = j, Var = var_k)
        
        current.e.mat = E.mat %>%
          as.numeric() %>%
          matrix(nrow = 1) %>%
          as.data.frame() %>%
          set_names(paste0("e", rep(1:7, times = 7), rep(1:7, each = 7))) %>%
          mutate(R0.origin = R0, sce = sce, j = j, Var = var_k)
        
        list(
          loop = current.loop.res,
          Kmat = current.Kmat,
          e.mat = current.e.mat
        )
        
        
      }, error = function(e) {
        NULL
      })
      
    }
    
    n.fail = sum(vapply(res, is.null, logical(1)))
    
    if (n.fail > 0) {
      warning(sprintf("%d of %d draws returned nothing in Var %s (%s, kappa %s, %s)",
                      n.fail, n.j, var_k, lna, kappa, sce), call. = FALSE)
    }
    
    sce.res[[k]] = list(
      loop = bind_rows(lapply(res, `[[`, "loop")),
      Kmat = bind_rows(lapply(res, `[[`, "Kmat")),
      e.mat = bind_rows(lapply(res, `[[`, "e.mat"))
    )
  }
  
  write_rds(sce.res, out.file)
  
  message("written: ", out.file, " (", length(sce.res), " Var)")
  
  invisible(out.file)
}


#=========== run the 27 scenarios =============

for (lna in names(lna.list)) {
  for (kappa in kappa.list) {
    for (sce in sce.list) {
      run_scenario(lna, kappa, sce)
    }
  }
}

stopCluster(cl)
