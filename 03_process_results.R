DIR_RES   = Sys.getenv("LOOP_RES_DIR", "./loop res")
OUT_STEM  = Sys.getenv("LOOP_SLIM_STEM", "./output/all.loop.slim")
OUT_EIJ   = Sys.getenv("LOOP_EIJ_LONG_OUT",  "./output/all.eij.wide.rds")

# ---- batching: 27 scenarios -> 3 parts of 9 ----------------------------------
BATCH      = 9
N_SPLIT    = 3
BATCH_EIJ  = 9

# -----------  the 27-scenario grid  -------------
# abundance ratio (L:N:A) x dilution rate kappa x tick abundance level
sce.list   = c("low", "medium", "high")
lna.list   = c("100:10:1", "200:10:1", "52:5:1")
kappa.list = c(0.15, 0.05, 0.5)

# file names only encode L (the larval part of the ratio), kappa and sce, e.g. 100_0.15_low_sample1000.rds
grid = expand_grid(lna = lna.list, kappa = kappa.list, sce = sce.list) %>%
  mutate(L    = str_remove(lna, ":.*"),
         file = file.path(DIR_RES, str_c(L, "_", kappa, "_", sce, "_sample1000.rds")),
         ok   = file.exists(file))

print(grid[, c("lna", "kappa", "sce", "ok")], n = 27)


#-------------------------------------------------------
# pass 1: loop table
# read one scenario -> keep the 7 columns needed downstream -> flush to disk every BATCH scenarios (one .rds per part)
#-------------------------------------------------------

# columns kept: loop identity, its elasticity point estimate and interval,R0 of the same draw, tick abundance level and location-year.
#-------------------------------------------------------

slim.cols = c("loop", "loop.elas", "loop.elas.low", "loop.elas.high",
              "R0.origin", "sce", "Var")

pieces = vector("list", BATCH)
k      = 0L
part   = 0L

for (i in seq_len(nrow(grid))) {

  g = grid[i, ]
  message(sprintf("%2d/%d  %s  [loop]", i, nrow(grid), basename(g$file)))

  x = read_rds(g$file)

  # one scenario = 38 location-years; bind them into a single long tableand tag the two scenario dimensions that are only in the file name
  piece = x %>%
    lapply(`[[`, "loop") %>%
    bind_rows() %>%
    select(all_of(slim.cols)) %>%
    mutate(ratio = g$lna,
           kappa = g$kappa)

  rm(x)

  k = k + 1L
  pieces[[k]] = piece
  rm(piece)

  # flush every BATCH scenarios, or at the end of the grid
  if (k == BATCH || i == nrow(grid)) {

    part = part + 1L
    nb   = rbindlist(pieces[seq_len(k)], use.names = TRUE)
    pieces = vector("list", BATCH)
    k = 0L
    invisible(gc())

    out.file = sprintf("%s_%d.rds", OUT_STEM, part)

    write_rds(nb, out.file)

    cat(sprintf("\npart %d: %s rows, lna = %s\n",
                part, format(nrow(nb), big.mark = ","),
                paste(unique(nb$ratio), collapse = " / ")))
    cat("written:", out.file,
        " (", round(file.size(out.file) / 1024^3, 2), " GB)\n")

    # sanity check: all 38 location-years and all 3 abundance levels present
    stopifnot(n_distinct(nb$Var) == 38,
              all(c("low", "medium", "high") %in% nb$sce))

    rm(nb); invisible(gc())
  }
}

cat("\nloop splits written:", part, "\n")

#-------------------------------------------------------
# pass 2: e_ij table
# same loop, but keep the whole elasticity matrix (e11:e77) of each draw.
#-------------------------------------------------------

pieces = vector("list", BATCH_EIJ)
k      = 0L
out    = NULL

for (i in seq_len(nrow(grid))) {

  g = grid[i, ]

  message(sprintf("%2d/%d  %s  [e.mat]", i, nrow(grid), basename(g$file)))

  x = read_rds(g$file)

  # e.mat is the 7x7 elasticity matrix per draw; lna and kappa come from
  # the file name, sce and Var are already columns of the table
  piece = x %>%
    lapply(`[[`, "e.mat") %>%
    bind_rows() %>%
    mutate(lna   = g$lna,
           kappa = g$kappa)

  rm(x)

  k = k + 1L
  pieces[[k]] = piece
  rm(piece)

  if (k == BATCH_EIJ || i == nrow(grid)) {
    nb  = rbindlist(pieces[seq_len(k)], use.names = TRUE)
    out = if (is.null(out)) nb else rbindlist(list(out, nb), use.names = TRUE)
    pieces = vector("list", BATCH_EIJ)
    k = 0L
    invisible(gc())
  }
}

all.eij.res = out
rm(out); invisible(gc())

# -------- checks on the e_ij table ---------
cat("\n===== all.eij.res =====\n")
cat("rows:", format(nrow(all.eij.res), big.mark = ","), " cols:", ncol(all.eij.res), "\n")
print(names(all.eij.res))
cat("location-years:", n_distinct(all.eij.res$Var), "\n")


eij = str_c("e", rep(1:7, each = 7), rep(1:7, times = 7))
cat("e11:e77 all present:", all(eij %in% names(all.eij.res)), "\n")

# draws surviving per scenario: 38 location-years x 1000 draws = 38000 max; anything below that is a failed draw dropped by the tryCatch in 02
chk = as_tibble(all.eij.res) %>% count(lna, kappa, sce)
cat("\neij rows per scenario (max 38 x 1000 = 38000):\n")
print(as.data.frame(chk), row.names = FALSE)

n.miss = 38000 * 27 - sum(chk$n)
cat("\ntotal rows:", format(sum(chk$n), big.mark = ","),
    " missing draws:", n.miss,
    " (", round(100 * n.miss / (38000 * 27), 2), "% )\n")

stopifnot(n_distinct(all.eij.res$Var) == 38,
          all(c("low", "medium", "high") %in% all.eij.res$sce),
          all(eij %in% names(all.eij.res)),
          nrow(chk) == 27,
          all(chk$n <= 38000))

print(head(all.eij.res[, c("Var", "sce", "j", "R0.origin", "e11", "e21", "lna", "kappa")], 3))

write_rds(all.eij.res, OUT_EIJ)

