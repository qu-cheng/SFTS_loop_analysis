library(readr)
library(tidyverse)
library(data.table)
library(lhs)
library(truncnorm)
library(EnvStats)
library(igraph)
library(igraph)
library(CVXR)
library(doParallel)
library(foreach)
library(cowplot)
library(ggsci)
library(ggtext)

all.dat = read_rds("./data/all_tick_dat.rds") %>% 
  rename("Var" = "clean.var")

alpha_L = alpha_N = alpha_A = 0.15

#===== parameters that DO NOT VARY with scenario =====

base.params = c("Sl", "Sn", "Sa", "E", "Ra",
                "theta",
                "Hcn", "Hcs",
                "Pl", "Pn", "Pa",
                "Dl","Dn", "Da",
                "Ql", "Qn", "Qa",
                "I", "Cs")

# sampling from LHS
base.unscaled = randomLHS(n = 1000, k = length(base.params))

colnames(base.unscaled) = base.params

base.scaled = base.unscaled

# Survival rate
base.scaled[, "Sl"] = qlnorm(base.unscaled[,"Sl"], meanlog = -3.1, sdlog = 0.4)
base.scaled[, "Sn"] = qlnorm(base.unscaled[,"Sn"], meanlog = -2.4, sdlog = 0.4)
base.scaled[, "Sa"] = qlnorm(base.unscaled[,"Sa"], meanlog = -2.4, sdlog = 0.4)

# Mean eggs
base.scaled[, "E"] = qtruncnorm(base.unscaled[,"E"], 0, Inf, mean = 1400, sd = 612)

# Efficiency from female adult to egg
base.scaled[, "Ra"] = 0.0343 + base.unscaled[, "Ra"]*0.9657

# Transmission efficiency from tick to tick
base.scaled[, "theta"] = 0.33 + base.unscaled[, "theta"]*0.67

# Efficiency from competent host to larva/nymph/adult
base.scaled[, "Pl"] = qtruncnorm(base.unscaled[,"Pl"], 0, Inf, 14/120, 0.03)
base.scaled[, "Pn"] = qtruncnorm(base.unscaled[,"Pn"], 0, Inf, 12/30, 0.0894)
base.scaled[, "Pa"] = qtruncnorm(base.unscaled[,"Pa"], 0, Inf, 6/13, 0.1382642)

# days of attachment
base.scaled[, "Dl"] = qtruncnorm(base.unscaled[,"Dl"], 0, Inf, 4.1, 0.5)
base.scaled[, "Dn"] = qtruncnorm(base.unscaled[,"Dn"], 0, Inf, 3.6, 0.9)
base.scaled[, "Da"] = qtruncnorm(base.unscaled[,"Da"], 0, Inf, 5.7, 1)

# Efficiency from larva/nymph/adult to competent host
base.scaled[, "Ql"] = qlnormTrunc(base.unscaled[,"Ql"], meanlog = -1.478410, sdlog = 0.77, min = 0, max = 1) 
base.scaled[, "Qn"] = qlnormTrunc(base.unscaled[,"Qn"], meanlog = -1.478410, sdlog = 0.77, min = 0, max = 1)
base.scaled[, "Qa"] = qlnormTrunc(base.unscaled[,"Qa"], meanlog = -0.6753073, sdlog = 0.41, min = 0, max = 1)

# Systemic infection duration (days)
base.scaled[, "I"] = 5 + base.unscaled[, "I"]*4  # (5-9 days)

# Mean fraction of ticks on a host feeding close enough for cofeeding
base.scaled[, "Cs"] = 0.1 + base.unscaled[, "Cs"]*0.4

# Fraction of blood meals on hosts competent for systemic transmission and non-systemic transmission
base.scaled[, "Hcn"] = 0.1 + base.unscaled[, "Hcn"]*0.8
base.scaled[, "Hcs"] = 0.1 + base.unscaled[, "Hcs"]*0.8

base.para = colMeans(base.scaled)

# base.para %>% write_rds("./data/base_para_.rds")

#================= parameters that VARY with scenarios ================

tick.unscaled = randomLHS(n = 1000, k = 3)

colnames(tick.unscaled) = c("low", "medium", "high")

tick.scaled = tick.unscaled

# Low, total abundance of 0.01~ 2, Liaoning
tick.scaled[, "low"] = 0.01 + tick.unscaled[, "low"] * (2 - 0.01)

# Medium, total abundance of 2~ 5, Shandong
tick.scaled[, "medium"] = 2 + tick.unscaled[, "medium"] * (5 - 2)

# High, total abundance of 5~12, Hubei
tick.scaled[, "high"] = 5 + tick.unscaled[, "low"] * (12 - 5)

# tick.scaled %>% write_rds("./data/1000_sets_tick_abundance.rds")

compute_R0 = function(dat){
  
  K.mat = matrix(0, 7, 7)
  
  base.dat = base.para %>%
    t %>% 
    as.data.frame() 
  
  attach(base.dat)
  
  K.mat[1,1] = Sl*Sn*Sa*E*Ra
  K.mat[1,2] = Sn*Sa*E*Ra
  K.mat[1,3] = Sa*E*Ra
  K.mat[1,4] = E*Ra
  
  cofeeding.dat = dat %>% 
    dplyr::select(C_LL:C_AN, -Pr_N0, -Pr_A0) %>% 
    distinct()
  
  attach(cbind(base.dat, cofeeding.dat))
  
  K.mat[2,1] = (Sl*C_LL + Sl*Sn*C_LN + Sl*Sn*Sa*C_LA)*theta*Cs*Hcn
  K.mat[2,2] = (Sn*C_LN + Sn*Sa*C_LA)*theta*Cs*Hcn
  K.mat[2,3] = (Sa*C_LA)*theta*Cs*Hcn
  
  K.mat[3,1] = (Sl*C_NL + Sl*Sn*C_NN + Sl*Sn*Sa*C_NA)*theta*Cs*Hcn
  K.mat[3,2] = (Sn*C_NN + Sn*Sa*C_NA)*theta*Cs*Hcn
  K.mat[3,3] = (Sa*C_NA)*theta*Cs*Hcn
  
  K.mat[4,1] = (Sl*C_AL + Sl*Sn*C_AN + Sl*Sn*Sa*C_AA)*theta*Cs*Hcn
  K.mat[4,2] = (Sn*C_AN + Sn*Sa*C_AA)*theta*Cs*Hcn
  K.mat[4,3] = (Sa*C_AA)*theta*Cs*Hcn
  
  K.mat[2,5] = compute_K_discrete(dat, "L", "L", Pl, Dl, I)
  K.mat[2,6] = compute_K_discrete(dat, "N", "L", Pl, Dl, I)
  K.mat[2,7] = compute_K_discrete(dat, "A", "L", Pl, Dl, I)
  
  
  K.mat[3,5] = compute_K_discrete(dat, "L", "N", Pn, Dn, I)
  K.mat[3,6] = compute_K_discrete(dat, "N", "N", Pn, Dn, I)
  K.mat[3,7] = compute_K_discrete(dat, "A", "N", Pn, Dn, I)
  
  K.mat[4,5] = compute_K_discrete(dat, "L", "A", Pa, Da, I)
  K.mat[4,6] = compute_K_discrete(dat, "N", "A", Pa, Da, I)
  K.mat[4,7] = compute_K_discrete(dat, "A", "A", Pa, Da, I)
  
  K.mat[5,1] = Sl*Ql*Hcn
  K.mat[6,1] = Sl*Sn*Qn*Hcn
  K.mat[6,2] = Sn*Qn*Hcn
  K.mat[7,1] = Sl*Sn*Sa*Qa*Hcn
  K.mat[7,2] = Sn*Sa*Qa*Hcn
  K.mat[7,3] = Sa*Qa*Hcn 
  
  R0 = Re(eigen(K.mat)$values[1])
  
  return(list(R0, K.mat))
}

compute_K_discrete = function(df, stage_from, stage_to, P, D, I) {
  
  a_col = paste0("a_", stage_from)
  z_col = switch(stage_to,
                 "L" = "Larva",
                 "N" = "Nymph",
                 "A" = "Adult")
  
  wrap = function(x) ((x - 1) %% 365) + 1
  
  K_val = 0
  for (t in 1:365) {
    a_val = df[[a_col]][t]
    inner_sum = sum(df[[z_col]][wrap(t:(t + I))])
    K_val = K_val + a_val * P / D * inner_sum
  }
  
  return(K_val)
}

# Guneralp
loopanalysis.Guneralp = function(gemat)
{
  g = graph_from_adjacency_matrix(t(gemat), mode = "directed", weighted = TRUE)
  # plot(g, layout=layout_as_star)
  
  Cycles = FindCycles(g)
  
  result = data.frame(loop = Cycles, length = NA, character.elas = NA, loop.elas = NA)
  
  current.gemat = gemat
  for(i in 1:nrow(result))
  {
    current.loop = result$loop[i]
    result$length[i] = length(unlist(regmatches(current.loop, gregexpr("\\b[[:digit:].]+\\b", current.loop)))) - 1
    
    current.result = chain.to.df(current.loop) %>%
      mutate(weight = NA,
             from = as.numeric(from),
             to = as.numeric(to))
    
    for(j in 1:nrow(current.result))
      current.result$weight[j] = current.gemat[current.result$to[j], current.result$from[j]]
    
    current.weight = min(current.result$weight)
    
    mat.reduce = matrix(0, nrow = nrow(gemat), ncol = ncol(gemat))
    for(j in 1:nrow(current.result))
      mat.reduce[current.result$to[j], current.result$from[j]] = current.weight
    
    if(min(current.gemat - mat.reduce) >= 0){
      current.gemat = current.gemat - mat.reduce
      result$character.elas[i] = current.weight
    }
  }
  
  result$loop.elas = result$character.elas * result$length
  result[result$character.elas != 0,]
}

chain.to.df = function(chain)
{
  chain.element = unlist(regmatches(chain, gregexpr("\\b[[:digit:].]+\\b", chain)))
  n = length(chain.element)
  
  data.frame(from = chain.element[1:(n-1)], to = chain.element[2:n])
}

## More efficient version
FindCycles = function(g) {
  Cycles = NULL
  for(v1 in V(g)) {
    if(degree(g, v1, mode="in") == 0) { next }
    GoodNeighbors = neighbors(g, v1, mode="out")
    GoodNeighbors = GoodNeighbors[GoodNeighbors > v1]
    for(v2 in GoodNeighbors) {
      TempCyc = lapply(all_simple_paths(g, v2,v1, mode="out"), function(p) c(v1,p))
      TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
      Cycles  = c(Cycles, TempCyc)
    }
  }
  
  Cycles = unlist(lapply(Cycles, function(x) paste(x, collapse = "->")))
  
  self_loop_indices = which(is.loop(g))
  self_loops = get.edgelist(g)[self_loop_indices, ]
  
  Cycles = c(Cycles, apply(self_loops, 1, function(x) paste(x, collapse = "->")))
  Cycles = Cycles[order(nchar(Cycles), Cycles)] # reorder Cycles from shortest to longest
  
  Cycles
}


# AdamS
loopanalysis.Adams = function(gemat, etol = 1)
{
  g = graph_from_adjacency_matrix(t(gemat), mode = "directed", weighted = TRUE)
  edge.list = get.edgelist(g)
  
  Cycles = FindCycles(g)
  
  Lambda = matrix(0, nrow = nrow(edge.list), ncol = length(Cycles))
  rownames(Lambda) = apply(edge.list, 1, function(x) paste(x, collapse = ","))
  colnames(Lambda) = Cycles
  
  result = data.frame(loop = Cycles, length = NA, character.elas.low = NA, character.elas.high = NA)
  
  for(i in 1:length(Cycles)){
    current.chain = Cycles[i]
    result$length[i] = length(unlist(regmatches(current.chain, gregexpr("\\b[[:digit:].]+\\b", current.chain)))) - 1
    
    current.chain = chain.to.df(current.chain) %>%
      mutate(from = as.numeric(from), to = as.numeric(to)) %>%
      as.matrix()
    
    # 合并两个矩阵
    combined_matrix = rbind(edge.list, current.chain)
    
    # 找到重复的行
    duplicated_rows = duplicated(combined_matrix, fromLast = TRUE)
    
    # 提取 matrix1 中重复的行索引
    row_indices = which(duplicated_rows[1:nrow(edge.list)])
    
    Lambda[row_indices, i] = 1
  }
  
  
  for(i in 1:length(Cycles))
  {
    # CVXR::installed_solvers()
    v = Variable(ncol(Lambda))
    objective = v[i]
    constraints = list(v >= 0, Lambda %*% v <= E(g)$weight, Lambda %*% v >= etol*E(g)$weight)
    
    prob = Problem(Maximize(objective), constraints)
    solution = solve(prob)
    result$character.elas.high[i] = solution$getValue(v)[i]
    
    prob = Problem(Minimize(objective), constraints)
    solution = solve(prob)
    result$character.elas.low[i] = solution$getValue(v)[i]
  }
  
  result %>%
    mutate(character.elas.low = ifelse(character.elas.low < 0, 0, character.elas.low),
           character.elas.high = ifelse(character.elas.high <0, 0, character.elas.high),
           loop.elas.low = character.elas.low*length,
           loop.elas.high = character.elas.high*length)
}


decomp.new = function(gemat){
  results = list()
  remainder = list()
  loops = list()
  i = 0
  results = loop.random.new(gemat)
  while (is.list(results) == TRUE) {
    i = i + 1
    gemat = results$newmat
    loops[[i]] = results$loopmat
    remainder = results$newmat
    results = loop.random.new(gemat)
  }
  loops[[length(loops) + 1]] = remainder
  return(loops)
}


loop.random.new = function (gemat) 
{
  nodenum = dim(gemat)[1]
  series = 1:nodenum
  gemat[abs(gemat) < 1e-10] = 0
  
  for (i in series) {
    if (gemat[i, i] > 0) {
      thisloop = c(i, i)
      weight = gemat[i, i]
      thismat = matrix(0, ncol = nodenum, nrow = nodenum)
      thismat[i, i] = weight
      return(list(loop = thisloop, weights = weight, loopmat = thismat, 
                  newmat = gemat - thismat))
      break
    }
    else {
      unchecked = series
      thisnode = i
      father = thisnode
      out.index = which(gemat[thisnode, ] > 0)
      out.index = out.index[which(out.index != thisnode)]
      in.index = which(gemat[, thisnode] > 0)
      in.index = in.index[which(in.index != thisnode)]
      thisloop = vector()
      weight = vector()
      flist = list()
      if (length(out.index) > 0 & length(in.index) > 0) {
        thisloop = c(thisloop, thisnode)
        flist[[length(flist) + 1]] = list(ansc = father, 
                                           desc = out.index)
        rid = sample(1:length(out.index), 1, replace = FALSE)
        thisnode = flist[[length(flist)]]$desc[rid]
        father = thisnode
        flist[[length(flist)]]$desc = flist[[length(flist)]]$desc[-rid]
        flag = 1
        while (gemat[thisnode, i] == 0 & flag == 1) {
          out.index = which(gemat[thisnode, ] > 0)
          out.index = out.index[which(out.index != thisnode)]
          if (length(out.index) > 0 & thisnode %in% unchecked == 
              TRUE) {
            thisloop = c(thisloop, thisnode)
            unchecked = unchecked[-thisnode]
            weight = c(weight, gemat[thisloop[length(thisloop) - 
                                                 1], thisloop[length(thisloop)]])
            flist[[length(flist) + 1]] = list(ansc = father, 
                                               desc = out.index)
            rid = sample(1:length(out.index), 1, replace = FALSE)
            thisnode = flist[[length(flist)]]$desc[rid]
            father = thisnode
            flist[[length(flist)]]$desc = flist[[length(flist)]]$desc[-rid]
            flag = 1
          }
          else {
            backone = length(flist)
            while (length(flist[[backone]]$desc) == 0) {
              flist[[backone]] = NULL
              thisloop = thisloop[-length(thisloop)]
              weight = weight[-length(weight)]
              backone = backone - 1
              if (backone == 0) {
                flag = 0
                break
              }
            }
            if (flag == 1) {
              rid = sample(1:length(flist[[backone]]$desc), 1, replace = FALSE)
              thisnode = flist[[backone]]$desc[rid]
              father = thisnode
              flist[[backone]]$desc = flist[[backone]]$desc[-rid]
            }
          }
          if (length(flist) > 0) {
            flag = 1
          }
          else {
            flag = 0
            break
          }
        }
        if (flag == 1) {
          thisloop = c(thisloop, thisnode, i)
          weight = c(weight, gemat[thisloop[length(thisloop) - 
                                               2], thisnode])
          weight = c(weight, gemat[thisnode, i])
          thismat = matrix(0, ncol = nodenum, nrow = nodenum)
          for (i in 1:(length(thisloop) - 1)) {
            thismat[thisloop[i], thisloop[i + 1]] = min(weight)
          }
          return(list(loop = thisloop, weights = weight, 
                      loopmat = thismat, newmat = gemat - thismat))
          break
        }
      }
    }
  }
}


calculate_elasticity = function(loop_matrix) {
  #转置矩阵
  loop_matrix = t(loop_matrix)
  # 创建有向图
  g = graph_from_adjacency_matrix(loop_matrix, mode = "directed", weighted = TRUE)
  
  # 提取所有边及其权重
  edges = get.edgelist(g)
  weights = E(g)$weight
  
  # 如果无边，返回空
  if (length(weights) == 0) return(NULL)
  
  # 处理自环（单节点环路）
  if (nrow(edges) == 1 && edges[1,1] == edges[1,2]) {
    path = paste(edges[1,1], edges[1,2], sep = " -> ")
    return(list(
      path = path,
      weight = min(weights),
      length = 1,
      elasticity = min(weights) * 1
    ))
  }
  
  # 寻找简单环路路径（假设每个矩阵只含一个环路）
  # 通过深度优先搜索找到闭合路径
  start_node = edges[1,1]
  current_node = start_node
  path_nodes = start_node
  visited_edges = logical(length(weights))
  
  for (i in seq_along(weights)) {
    idx = which(edges[,1] == current_node & !visited_edges)
    if (length(idx) == 0) break
    
    next_node = edges[idx[1], 2]
    path_nodes = c(path_nodes, next_node)
    visited_edges[idx[1]] = TRUE
    current_node = next_node
    
    # 检查是否回到起点形成环路
    if (next_node == start_node) break
  }
  
  # 验证路径闭合
  if (path_nodes[length(path_nodes)] != start_node) {
    return(NULL) # 无法形成闭合环路
  }
  
  # 生成路径字符串
  path = paste( path_nodes, collapse = " -> ")
  
  # 计算指标
  valid_weights = weights[visited_edges]
  loop_weight = min(valid_weights)
  loop_length = length(valid_weights)
  elasticity = loop_weight * loop_length
  
  return(list(
    path = path,
    weight = loop_weight,
    length = loop_length,
    elasticity = elasticity
  ))
}



