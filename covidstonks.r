# Covidstonks code
# Library installation and loading
packages = c("xml2",
             "rvest",
             "RCurl",
             "httr",
             "quantmod",
             "tidyquant",
             "stringr",
             "huge",
             "qgraph",
             "CINNA",
             "ggplot2",
             "igraph",
             "tidyverse",
             "Matrix",
             "doParallel",
             "glasso",
             "sm",
             "devtools",
             "foreach",
             "loggle")
install.packages(packages[packages != "loggle"])
install_github(repo="jlyang1990/loggle")
lapply(packages, require, character.only=TRUE)
g_stat = function(graph, decay_closeness = 0.5, decay_prestige = 0.5, include_prestige = FALSE){
  # Centrality measures
  # For any graph, you need to define the following objects:
  g <- graph.adjacency(graph)
  g_mat <- as.matrix(as_adjacency_matrix(g))
  # The as.matrix has to go twice if you have a sparse matrix
  
  # Degree centrality : number of edges / number of total nodes
  g_degree <- degree(g) # This stores the number of edges for every node
  g_degree_centrality <- g_degree/(dim(g_mat)[1] - 1) # Normalized degrees
  g_degree_ranking <- c(1:dim(g_mat)[1])[order(-g_degree_centrality)]
  # To identify who these nodes are in stock data, just do:
  stocks_degree_ranking <- sp500.ticker.all[g_degree_ranking]
  g_degree_data <- list(measure = g_degree, centrality = g_degree_centrality,
                        ranking = g_degree_ranking,
                        stocks_ranking = stocks_degree_ranking)
  
  # (Decay) Closeness centrality : Inverse avg. distance between other nodes
  delta <- 0.5 # This parameter determines the weight (or decay) of closeness
  # delta approaching 1 : largeness of node's component
  # delta approaching 0 : proportional to degree centrality
  g_distance <- distances(g) # This stores distance between i and j
  g_decay_centrality <- c()
  for(x in c(1:dim(g_mat)[1])){
    g_decay_centrality[x] <- sum(delta^g_distance[x,])
  }
  g_decay_ranking <- c(1:dim(g_mat)[1])[order(-g_decay_centrality)]
  # To identify who these nodes are in stock data, just do:
  stocks_decay_ranking <- sp500.ticker.all[g_decay_ranking]
  g_decay_data <- list(measure = g_distance, centrality = g_decay_centrality, 
                       ranking  = g_decay_ranking,
                       stocks_ranking = stocks_decay_ranking)
  
  # Betweenness centrality : number of paths that node lies on
  g_betweenness <- betweenness(g)
  g_betweenness_ranking <- c(1:dim(g_mat)[1])[order(-g_betweenness)]
  # To identify who these nodes are in stock data, just do:
  stocks_betweenness_ranking <- sp500.ticker.all[g_betweenness_ranking]
  g_betweenness_data <- list(measure = g_betweenness, ranking = g_betweenness_ranking,
                             stocks_ranking = stocks_betweenness_ranking)
  
  # Bonacich prestige centrality : importance of neighbors
  if(include_prestige == TRUE){
    g_prestige <- power_centrality(g, exponent = 0.5)
    # The exponent parameter sets the decay of prestige over distances
    g_prestige_ranking <- c(1:dim(g_mat)[1])[order(-g_prestige)]
    # To identify who these nodes are in stock data, just do:
    stocks_prestige_ranking <- sp500.ticker.all[g_prestige_ranking]
    g_prestige_data <- list(measure = g_prestige, ranking = g_prestige_ranking,
                            stocks_ranking = stocks_prestige_ranking)
  }
  
  # Sparsity
  g_sparsity <- sum(g_mat)/(dim(g_mat)[1]*(dim(g_mat)[1]-1))
  
  # Influence
  library(CINNA)
  maxcomp <- giant_component_extract(g) 
  maxcompedgenum <- dim(maxcomp[[2]])[1] #num of edges
  maxcompnodenum <- length(V(maxcomp[[1]])) #num of vertices 
  nodesize_time <- length(V(g))
  
  maxcompinfluence_edges <- maxcompedgenum / sum(g_mat)
  maxcompinfluence_nodes <- maxcompnodenum / nodesize_time
  g_influence_data <- list(maxcomp_graph = maxcomp, maxcomp_edges = maxcompedgenum, 
                           maxcomp_nodes = maxcompnodenum, 
                           PTE = maxcompinfluence_edges, PTN = maxcompinfluence_nodes)
  g_data <- list(degree = g_degree_data, decay_closeness = g_decay_data, 
                 betweenness = g_betweenness_data,
                 prestige = if(include_prestige == TRUE){g_prestige_data}else{"Not computed."}, 
                 sparsity = g_sparsity, 
                 giant_influence = g_influence_data)
  return(g_data)
}
# Scraping
sp500list_page <- "https://en.wikipedia.org/wiki/List_of_S%26P_500_companies" 
sp500.file <- read_html(sp500list_page)
sp500.tables <- html_nodes(sp500.file, "table")
sp500.df <- html_table(sp500.tables[1], fill = TRUE)[[1]]
sp500.ticker.all <- sp500.df$Symbol
sp500.ticker.all <- gsub('\\.','-',sp500.ticker.all)
sp500.name.all <- sp500.df$Security
sp500.sector.all <- sp500.df$'GICS Sector'
# Data acquisition
sp500.quantdat <- tq_get(sp500.ticker.all, from = '1981-03-20',
                         to = '2021-02-28', do.parallel = TRUE) # 1981 to 2021
saveRDS(sp500.quantdat, file = "sp500_data.rds")
# Data cleaning
transformation_type = "unbalanced" # Options: "balanced", "unbalanced", or "censored". All other values make the code not run.
# Note, doing unbalanced can always be reverted to being balanced or censored so there is no loss of data in getting the unbalanced transformation
# The only issue is that it is more time consuming than may be desired.

# Balanced data transformation (Only data with observations present for all dates is included)
if(transformation_type == "balanced"){
  print("Transforming data... (balanced)")
  sp500.all <- data.frame(date = unique(sp500.quantdat$date)) # Creates a dataframe with all of the dates
  for(i in sp500.ticker.all){ # For every stock in the scraped list
    if(length(subset(sp500.quantdat$close, sp500.quantdat$symbol == i)) == 
       length(sp500.all$date)){ # If there's as many observations correspondent to this stock as there are dates (so the data is complete)
      sp500.all$new <- subset(sp500.quantdat$close, sp500.quantdat$symbol == i) # Then, make a new column for this stock with all the observations
      colnames(sp500.all)[colnames(sp500.all) == "new"] <- i } # Then finally, rename this column according to the stock name
  }
  rownames(sp500.all) <- sp500.all$date # Once this is complete, set the rownames to be the dates themselves
  sp500.all <- subset(sp500.all, select = -c(date)) # And to avoid redundancy, get rid of the column with all the dates.
}
# Unbalanced data transformation (All observations included with NAs for missing observations, larger files than censored and slower.)
if(transformation_type == "unbalanced"){
  print("Transforming data... (unbalanced)")
  sp500.all <- data.frame(date = unique(sp500.quantdat$date)) # Creates a dataframe with all of the dates
  for(i in sp500.ticker.all){ # For every stock in the scraped list
    if(length(subset(sp500.quantdat$close, sp500.quantdat$symbol == i)) == 
       length(sp500.all$date)){ # If there's as many observations correspondent to this stock as there are dates (so the data is complete)
      sp500.all$new <- subset(sp500.quantdat$close, sp500.quantdat$symbol == i) # Then, make a new column for this stock with all the observations
      colnames(sp500.all)[colnames(sp500.all) == "new"] <- i } # Then rename this column according to the stock name
    else{ # If there are more dates than observations for this stock (so the data is incomplete)
      lastmissingindex <- which(sp500.all$date == min(subset(sp500.quantdat$date, sp500.quantdat$symbol == i))) - 1 # Find the first date that an observation appears for this stock
      sp500.all$new <- NA # Create a new column with nothing but NAs
      sp500.all$new[lastmissingindex:length(sp500.all$new)] <- subset(sp500.quantdat$close, sp500.quantdat$symbol == i) # From this date onwards, put in the data
      colnames(sp500.all)[colnames(sp500.all) == "new"] <- i # Rename this column to stock name
    }
  }
  rownames(sp500.all) <- sp500.all$date # Once this is complete, set the rownames to be the dates themselves
  sp500.all <- subset(sp500.all, select = -c(date)) # And to avoid redundancy, get rid of the column with all the dates.
}
# Censored data transformation (All observations included with 1s for missing observations, smaller files than unbalanced and faster.)
if(transformation_type == "censored"){
  print("Transforming data... (censored)")
  sp500.all <- data.frame(date = unique(sp500.quantdat$date)) # Creates a dataframe with all of the dates
  for(i in sp500.ticker.all){ # For every stock in the scraped list
    if(length(subset(sp500.quantdat$close, sp500.quantdat$symbol == i)) == 
       length(sp500.all$date)){ # If there's as many observations correspondent to this stock as there are dates (so the data is complete)
      sp500.all$new <- subset(sp500.quantdat$close, sp500.quantdat$symbol == i) # Then, make a new column for this stock with all the observations
      colnames(sp500.all)[colnames(sp500.all) == "new"] <- i } # Then rename this column according to the stock name
    else{ # If there are more dates than observations for this stock (so the data is incomplete)
      lastmissingindex <- which(sp500.all$date == min(subset(sp500.quantdat$date, sp500.quantdat$symbol == i))) - 1 # Find the first date that an observation appears for this stock
      sp500.all$new <- 1 # Create a new column with nothing but 1s
      sp500.all$new[lastmissingindex:length(sp500.all$new)] <- subset(sp500.quantdat$close, sp500.quantdat$symbol == i) # From this date onwards, put in the data
      colnames(sp500.all)[colnames(sp500.all) == "new"] <- i # Rename this column to stock name
    }
  }
  rownames(sp500.all) <- sp500.all$date # Once this is complete, set the rownames to be the dates themselves
  sp500.all <- subset(sp500.all, select = -c(date)) # And to avoid redundancy, get rid of the column with all the dates.
}
saveRDS(sp500.all, file = "sp500_transformed.rds")
# Log returns transformation
X <- log(sp500.all[-1,]/sp500.all[-dim(sp500.all)[1],])
rownames(X) = rownames(sp500.all[-1,])
N <- nrow(X)
p <- ncol(X)
# Date matrix
X.date <- as.Date(rownames(X), "%Y-%m-%d")
year = as.numeric(format(X.date, "%Y"))
month = as.numeric(format(X.date, "%m"))
day = as.numeric(format(X.date, "%d"))
quarter = c(rep(0,length(rownames(X))))
X.date <- cbind(year, month, day, quarter)
X.date <- as_tibble(X.date)

for (i in 1:length(X.date$month)){
  if (X.date$month[i] < 10) {
    X.date$quarter[i] = 3
  }
  if (X.date$month[i] < 7){
    X.date$quarter[i] = 2
  }
  if (X.date$month[i] < 4){
    X.date$quarter[i] = 1
  }
  if (X.date$month[i] > 9){
    X.date$quarter[i] = 4
  }
}
# Loggle time stamps
yearly_loggle_stamps <- c()
quarterly_loggle_stamps <- c()
for(y in unique(X.date$year)){
  yearly_loggle_stamps[which(unique(X.date$year) == y)] <- max(rownames(X.date)[X.date$year == y])
  for(q in unique(X.date$quarter[X.date$year == y])){
    quarterly_loggle_stamps[q + 4*(which(y == unique(X.date$year)) - 1)] <- max(rownames(X.date)[X.date$year == y & X.date$quarter == q])
  }
}
yearly_loggle_stamps
quarterly_loggle_stamps
# Estimation (big boy!)
numCores <- parallel::detectCores()
ts <- proc.time()
year_loggle <- loggle.cv(X, yearly_loggle_stamps,
                         h.list = c(0.1,0.15),
                         d.list = c(0,0.001,0.01,
                                    0.025,0.05,
                                    0.075,0.1,0.15,
                                    0.2,0.25,0.3,1),
                         lambda.list = 10^seq(-2,0,0.1),
                         num.thread = numCores)
te <- proc.time()
sprintf("Time used for year loggle.cv: %.2fs", (te-ts)[3])
ts <- proc.time()
saveRDS(year_loggle, file="year_loggle.rds",compress=FALSE)
te <- proc.time()
sprintf("Time used to save year loggle.cv: %.2fs", (te-ts)[3])
ts <- proc.time()
quarter_loggle <- loggle.cv(X, quarterly_loggle_stamps,
                         h.list = c(0.1,0.15),
                         d.list = c(0,0.001,0.01,
                                    0.025,0.05,
                                    0.075,0.1,0.15,
                                    0.2,0.25,0.3,1),
                         lambda.list = 10^seq(-2,0,0.1),
                         num.thread = numCores)
te <- proc.time()
sprintf("Time used for quarterly loggle.cv: %.2fs", (te-ts)[3])
ts <- proc.time()
saveRDS(quarter_loggle, file="quarter_loggle.rds",compress=FALSE)
te <- proc.time()
sprintf("Time used to save quarterly loggle.cv: %.2fs", (te-ts)[3])
