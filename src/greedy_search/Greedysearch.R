## BGe score
logBGe <- function(Data, incidence, v=1,a=nrow(Data)+2,
                   mu=numeric(nrow(Data)), T_0=diag(0.5,nrow(Data),nrow(Data))){
  n <- nrow(Data)
  m <- ncol(Data)
  T_m <- T_0 + (m-1)* cov(t(Data)) + ((v*m)/(v+m))* (mu - rowMeans(Data))%*%t(mu - rowMeans(Data))
  c_function <- function(N,A){
    fact <- numeric(N)
    for (i in 1:N){
      fact[i] <- -lgamma((A+1-i)/2)
    }
    product <- sum(fact) -(A*N/2)*log(2)- (N*(N-1)/4)*log(pi)
    return(product)}
  P_local_num <- numeric(n)
  P_local_den <- numeric(n)
  for (j in 1:n) {
    n_nodes <- which(incidence[,j]==1)
    P_local_num[j] <- (-(length(n_nodes)+1)*m/2)*log(2*pi) + ((length(n_nodes)+1)/2)*log(v/(v+m)) + c_function((length(n_nodes)+1),a)-c_function((length(n_nodes)+1),a+m)+ (a/2)*log(det(as.matrix(T_0[sort(c(n_nodes,j)),sort(c(n_nodes,j))])))+ (-(a+m)/2)*log(det(as.matrix(T_m[sort(c(n_nodes,j)),sort(c(n_nodes,j))])))
    if(sum(incidence[,j])>0){
      P_local_den[j] <- (-(length(n_nodes))*m/2)*log(2*pi) + (length(n_nodes)/2)*log(v/(v+m)) + c_function(length(n_nodes),a)- c_function(length(n_nodes),a+m)+ (a/2)*log(det(as.matrix(T_0[n_nodes,n_nodes])))+ (-(a+m)/2)*log(det(as.matrix(T_m[n_nodes,n_nodes])))
    }
    else{
      P_local_den[j] <- 0
    }
  }
  log_bge <- sum(P_local_num) - sum(P_local_den)
  return(log_bge)
}
##
# a: Incidence to Ancestor matrix
incidence.to.ancestor <- function(I) {
  ones <- matrix(rep(1, nrow(I)**2), nrow = nrow(I), byrow = T)
  E <- t(expm::expm(I) - diag(nrow(I)))
  pmin(ceiling(E), ones)
}
# b: output a list of all edges that can be deleted
# without introducing a cycle (which is all the edges of the graph)
valid.deletions <- function(I){
  which(I!=0, arr.ind = T)
  1
}
# c: outputs a list of all edges that can be added
# without introducing a cycle
valid.additions <- function(I){
  ones <- matrix(rep(1, nrow(I)**2), nrow = nrow(I), byrow = T)
  identity <- diag(nrow(I))
  A <- incidence.to.ancestor(I)
  which (ones - identity - I -A == 1, arr.ind = T)
}
# d: outputs a list of all edges that can be reversed
# without introducing a cycle
valid.reversals <- function(I){
  A <- incidence.to.ancestor(I)
  which(I - t(A%*%t(I))== 1, arr.ind = T)
}
# e: outputs a list of all single edge operations
# that give a valid DAG
valid.ops <- function(I) {
  deletions <- valid.deletions(I) %>% as_tibble() %>%
    mutate(type = 1)
  additions <- valid.additions(I) %>% as_tibble() %>%
    mutate(type = 2)
  reversals <- valid.reversals(I) %>% as_tibble() %>%
    mutate(type = 3)
  bind_rows(deletions, additions, reversals)
}
# f: outputs a list of the incidence matrices of
# all valid neighbour graphs.
ops.to.incidence <- function(row, col, ops, I){
  I[row, col] <- case_when (
    ops == 1 ~ 0, # deletion, reversal
    ops == 2 ~ 1,
    ops == 3 ~ 0# addition
  )
  if (ops == 3)
    I[col, row] <- 0 # reversal
  I
}
incidence.from.ops <- function(I){
  valids <- valid.ops(I)
  incidences <- pmap(list(valids$row, valids$col, valids$type),
                     ops.to.incidence, I = I )
  valids$I <- incidences
  valids
}
#### And from Marco,
cpdag <- function(incidence){
  z <- order.edges(incidence)
  new_mat <- cbind(z,numeric(nrow(z))) # edges, parents, children, order, zeros
  n_mat <- new_mat[order(new_mat[,4]),] # sort the edges by its order
  vec <- numeric(nrow(z))
  while(any(vec==0)){ # while there are unlabeled edges l.3
    if (length(vec)>1){ # if there are at least 2 edges
      first <- which(n_mat[,5]==0)[1] # first EDGE that ist labeled "unknown" (0) l.4
      parent1 <- n_mat[first,2] # x parent NODE
      child1 <- n_mat[first,3] # y child NODE
      comp1 <- n_mat[which(n_mat[,3]==parent1 & n_mat[,5]==1),2] # w NODES that have an edge incident into the parent labeled compelled)
    }
    if (length(vec)==1){
      first <- which(n_mat[5]==0) # first edge that ist labeled "unknown" (0)
      parent1 <- n_mat[2] # x parent
      child1 <- n_mat[3] # y child
      comp1 <- numeric(0)
    }
    for (j in comp1){ # l.5
      if (incidence[j,child1]==0){ # if w is not a parent of the child l.6
        n_mat[first,5] <- 1 # label x -> y compelled l.7
        n_mat[which(n_mat[,3]==child1),5] <- 1 # label every edge incident into y compelled l.7
        vec[first] <- 1
        vec[which(n_mat[,3]==child1)] <- 1
        break
      }
      if (incidence[j,child1]!=0) {
        n_mat[which(n_mat[,2]==j & n_mat[,3]==child1),5] <- 1 # label w -> y compelled l.10
        vec[which(n_mat[,2]==j & n_mat[,3]==child1)] <- 1
      }
    }
    if (length(vec)>1){
      if(n_mat[first,5]==0){
        moep <- n_mat[which(n_mat[,3]==child1 & n_mat[,2]!=parent1),2] # other parents of the child
        if(length(moep)>0){ # l.11
          for(o in moep){
            if(incidence[o,parent1]==0){
              vec[first] <- 1
              vec[which(n_mat[,3]==child1 & n_mat[,5]==0)] <- 1
              n_mat[first,5] <- 1 # label x -> y compelled
              n_mat[which(n_mat[,3]==child1 & n_mat[,5]==0),5] <- 1 # label all "unknown" edges incident into y compelled
              break
            }
            if(all(incidence[moep,parent1]!=0)){
              vec[first] <- -1
              vec[which(n_mat[,3]==child1 & n_mat[,5]==0)] <- -1
              n_mat[first,5] <- -1 # label x -> y reversible
              n_mat[which(n_mat[,3]==child1 & n_mat[,5]==0),5] <- -1 # label all "unknown" edges incident into y reversible
            }
          }
        }
        3
        if(length(moep)==0){
          vec[first] <- -1
          vec[which(n_mat[,3]==child1 & n_mat[,5]==0)] <- -1
          n_mat[first,5] <- -1 # label x -> y reversible
          n_mat[which(n_mat[,3]==child1 & n_mat[,5]==0),5] <- -1 # label all "unknown" edges incident into y reversible
        }
      }
    }
    if (length(vec)==1){
      n_mat[5] <- -1 # label x -> y reversible
      vec <- -1
    }
  }
  return(n_mat)
}
top_order <- function(incidence){
  n <- length(incidence[1,])
  Order <- numeric(n)
  fan_in <- numeric(n)
  no_fan_in <- numeric(0)
  m <- 1
  for (p in 1:n){ # number of parent nodes at the beginning
    fan_in[p] <- sum(incidence[,p])
  }
  no_fan_in <- which(fan_in==0)
  while (length(which(Order==0))>0){ # as long as there is a node without an order
    fan_in[which(incidence[no_fan_in[1],]==1)] <- fan_in[which(incidence[no_fan_in[1],]==1)] - 1
    no_fan_in <- c(no_fan_in, c(which(incidence[no_fan_in[1],]==1),which(fan_in==0))[duplicated(c(which(incidence[no_fan_in[1],]==1),which(fan_in==0)))])
    Order[m] <- no_fan_in[1]
    no_fan_in <- no_fan_in[-1]
    m <- m+1
  }
  return(Order)
}
order.edges <- function(incidence){
  top.order <- top_order(incidence)
  n <- length(top.order)
  edges <- which(incidence!=0)
  children <- child(edges,n)
  parents <- parent(edges,n)
  m <- length(edges)
  ordered_edges <- numeric(m)
  incidence_n <- incidence
  tog <- matrix(c(edges,parents,children,ordered_edges),ncol=4, byrow=FALSE)
  k <- 1
  while(any(tog[,4]==0)){
    node1 <- top.order[which(colSums(incidence_n[,top.order])>0)][1] # first node in top. order that has at least one parent
    par1<- tog[which(tog[,3]==node1),2] # find the parents of first child in the top. order that has an unordered edge incident into it
    4
    g <- par1[which(par1>0)]
    f1 <- numeric(length(g))
    for (i in 1:length(g)){
      f1[i] <- which(top.order==g[i])
    }
    par2 <- g[which.max(f1)] # find the highest ordered node that has an edge leading into node1
    tog[which(tog[,2]==par2 & tog[,3]==node1),4] <- k
    k <- k + 1
    incidence_n[tog[which(tog[,2]==par2 & tog[,3]==node1),1]] <- 0 # delete the edge in the "incidence" matrix
    tog[which(tog[,2]==par2 & tog[,3]==node1),2] <- 0
  }
  to <- matrix(c(edges,parents,children,tog[,4]),ncol=4,byrow=FALSE)
  return(to) # return the whole matrix, the order is the fourth column
}
############################################################
child <- function(edges,n){
  # input: the numbers of the edges in the incidence matrix and the number of nodes
  p <- ceiling(edges/n)
  return(p)
}
##############################################################
parent <- function(edges,n){
  ch <- edges + n - child(edges,n)*n
  return(ch)
}

greedy_search <- function(obs, init.I){
    best.neighbor <- init.I
    old_score <- logBGe(obs, init.I)
    i <- 0
    while (abs(old_score)>0){
      i <- i+1
      N <- incidence.from.ops(best.neighbor)
      scores <- unlist(map(N$I, logBGe, Data = obs))
      if (max(scores)[1] > old_score)
        best.neighbor <- N$I[which(scores == max(scores))][[1]]
      if (max(scores)[1] == old_score)
        break
      old_score <- max(scores)[1]
    }
    best.neighbor
  }
