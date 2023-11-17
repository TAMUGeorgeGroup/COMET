# Create an environment to store global variables
my_package_env <- new.env()

#' Read necessary files and parameters for pipeline to run
#' This function reads all the necessary files to run this pipeline
#' @param tables.dir directory with a csv file that had the address
#' for the data and metadata
#' @param input.data.dir this is the directory where the input file
#' should be saved in
#'
#' @return nothing, stores necessary parameters for the model in global
#' variables
#' @importFrom readxl read_excel
#' @importFrom utils read.csv
#' @export
start_pipeline<-function(tables.dir, input.data.dir){
  # Read files and store in the environment
  my_package_env$EMT.genes <- read.csv(system.file("extdata", "EMTGeneListForMLR.csv", package = "COMET"))
  my_package_env$EMT.genes$name <- my_package_env$EMT.genes

  hold.xl <- readxl::read_excel(system.file("extdata", "EM_gene_signature_cellLine_KS.xlsx", package = "COMET"), col_names = FALSE)

  my_package_env$EMT.sig <-data.frame(hold.xl)

  my_package_env$color.scheme <- c("#191935", "#1B1B3A", "#2F2246", "#422951", "#693668", "#e3f6f5", "#bae8e8", "#2c698d")
  my_package_env$emt.color.scheme <- c('#CDF0EA', '#F7DBF0', '#BEAEE2') # Mesenchymal, Hybrid, Epithelial
  my_package_env$emt.color.scheme.bold <- c("#24A19C","#D96098",  "#BEAEE2")

  my_package_env$flow.dat <- read.csv(system.file("extdata", "markov_flow.csv", package = "COMET"))
  my_package_env$flow.dat[1:7, ] <- my_package_env$flow.dat

  my_package_env$inc_num <- 0.04
  my_package_env$CTMC.scale <- 3/20
}




#' This function runs the data driven pipeline for inferring trajectories
#' @param data.inputs the input datasheet stored in a csv file in the tables dir
#' @param tables.dir directory with a csv file that had the address
#' for the data and metadata
#' @param input.data.dir this is the directory where the input file
#' should be saved in
#' @param cutoff of highly variable EMT genes to be considered
#'
#' @return the inferred trajectories, also saves them within the
#' COMET_populated_files directory
#' @importFrom data.table fread
#' @importFrom phateR library.size.normalize
#' @importFrom stats princomp
#' @importFrom Seurat VariableFeatures
#' @importFrom reshape2 melt
#' @importFrom utils write.csv
#' @importFrom umap umap
#' @importFrom stats princomp kmeans
#' @importFrom Rmagic magic
#' @importFrom utils head
#' @export
run_pipeline<- function(data.inputs, tables.dir, input.data.dir, cutoff){
  .<-model.cluster<-NULL
  .<-cluster<-NULL
  .<-clust.1.lab<-NULL
  .<-clust.2.lab<-NULL
  .<-clust.3.lab<-NULL
  .<-meta.data.Time<-NULL
  for(k in 1:nrow(data.inputs)){
    for(turn in 1:10){
    data.inputs[k,]->data.input #For demonstration purposes we only read the first file
    data <- fread(file.path(paste0(input.data.dir,data.input$DataPath)))
    data$V1->gene.names #Extract gene names
    meta.data <- fread(file.path(paste0(input.data.dir,data.input$MetaData)))

    sort(unique(meta.data$Time))->time.points
    time.points->timepoints

    #Step I. library size normalization
    rownames(data)<-gene.names #Set the gene names
    data$V1 <- NULL
    t(data)->data #Transpose data
    colnames(data)<- gene.names
    rownames(data)<-as.character(seq(1, length(rownames(data))))#Change rownames to numbers so we can feed to MAGIC

    #Filter cells having at least 10 expressed genes
    keep_cols <- colSums(data > 0) > 10
    data.filtered.not.norm <- data[,keep_cols]

    #library size normalize
    data.filtered <- library.size.normalize(data)

    #Find highly variable genes
    as.data.frame(library.size.normalize(data.filtered.not.norm))->data.df
    data.df[colnames(data.df) %in% my_package_env$EMT.genes]->data.df

    #Choose a specific cutoff of highly variable genes
    Seurat::CreateSeuratObject(t(data.df))->seur
    Seurat::FindVariableFeatures(seur, selection.method='vst', nfeatures=200)->seur
    top.var.genes <- utils::head(Seurat::VariableFeatures(seur), cutoff)

    data_MAGIC <- Rmagic::magic(data.filtered, genes=top.var.genes,
                                knn=15)
    #Add time to MAGIC data
    as.data.frame(data_MAGIC)->data.magic
    data.magic$time <- meta.data$Time


    for(i in 1:length(timepoints)){
      meta.data$color[meta.data$Time==timepoints[i]]<-my_package_env$color.scheme[i]
      meta.data$Time[meta.data$Time==timepoints[i]]<-timepoints[i]
    }

    umapped <- umap::umap(data.magic[1:(length(colnames(data.magic))-1)])
    meta.data$Time->my.colors

    for(i in 1:length(timepoints)){
      my.colors[my.colors==timepoints[i]]<-my_package_env$color.scheme[i]
    }
    #We are considering only 3 states with 1 hybrid
    states.no <- 3
    stats::kmeans(umapped$layout, states.no)->model

    data.frame(as.data.frame(umapped$layout)$V1, as.data.frame(umapped$layout)$V2, model$cluster)->df.model

    for(state in 1:states.no){
      name <- paste0("df.model.", state)
      df.model %>% filter(model.cluster==state)->df.model.hold
      assign(name, df.model.hold)
    }

    pca_x <- princomp(df.model[1:2])
    x_cluster <- data.frame(pca_x$scores,model$cluster, meta.data$color, meta.data$Time)
    x_cluster[order(as.numeric(x_cluster$meta.data.Time)),]->x_cluster

    model$cluster->meta.data$cluster
    data.magic$cluster<-meta.data$cluster
    top100 <- utils::head(VariableFeatures(seur), 100)
    ####Labeling the clusters with KS test
    data.df$cluster<-data.magic$cluster

    #Run MAGIC
    data_MAGIC <- Rmagic::magic(data.filtered, genes=top100,
                                knn=15)
    as.data.frame(data_MAGIC)->data.magic
    data.magic$time <- meta.data$Time
    data.magic$cluster<-meta.data$cluster

    x_cluster$ks.scores<-0

    for(state in 1:states.no){
      name <- paste0("clust.", state, ".lab")
      (data.magic %>% filter(cluster==state))[1:(dim(data.magic)[2]-2)]->exp.mat
      colnames(exp.mat)->genes
      KS.label.me(exp.mat, genes, top100)->clust.lab
      assign(name, clust.lab)
      x_cluster[x_cluster$model.cluster==state,]$ks.scores<-clust.lab
    }

    #you would need to modify this if adding more states - this is hardcoded for 3 states
    clust.labs <- data.frame(type=c(rep("1", length(clust.1.lab)), rep("2", length(clust.2.lab)), rep("3", length(clust.3.lab))),
                             value = c(clust.1.lab, clust.2.lab, clust.3.lab))



    as.matrix(as.numeric(clust.labs$value))->my.mat
    colnames(my.mat)<-"KS.score"
    rownames(my.mat)<-as.character(seq(1, length(my.mat[,1])))
    cbind(my.mat, as.numeric(clust.labs$type))->my.mat
    colnames(my.mat)<-c("KS.score", "Cluster")
    colnames(clust.labs)<-c("cluster", "score")
    clust.labs[order(clust.labs$score),]->clust.labs

    #again hardcoded for three states modify accordingly
    x_cluster[x_cluster$model.cluster==unique(clust.labs$cluster)[1],]$model.cluster<-"Epithelial"
    x_cluster[x_cluster$model.cluster==unique(clust.labs$cluster)[2],]$model.cluster<-"Hybrid"
    x_cluster[x_cluster$model.cluster==unique(clust.labs$cluster)[3],]$model.cluster<-"Mesenchymal"

    #x_cluster$meta.data.Time <- factor(x_cluster$meta.data.Time,                 # Relevel group factor
    #                                   levels = as.numeric(x_cluster$meta.data.Time))
    clusterM.day <- NA
    clusterH.day <- NA
    clusterE.day <- NA
    for(i in 1:length(unique(x_cluster$meta.data.Time))){
      timepoint <- unique(x_cluster$meta.data.Time)[i]
      x_cluster %>% filter(meta.data.Time==timepoint)->data.time
      nrow(data.time %>% filter(model.cluster=="Mesenchymal"))/nrow(data.time)->clust.M.perc
      nrow(data.time %>% filter(model.cluster=="Hybrid"))/nrow(data.time)->clust.H.perc
      nrow(data.time %>% filter(model.cluster=="Epithelial"))/nrow(data.time)->clust.E.perc
      clusterM.day <- c(clusterM.day, clust.M.perc)
      clusterH.day <- c(clusterH.day, clust.H.perc)
      clusterE.day <- c(clusterE.day, clust.E.perc)
    }
    clusterM.day[2:length(clusterM.day)]->clusterM.day
    clusterH.day[2:length(clusterH.day)]->clusterH.day
    clusterE.day[2:length(clusterE.day)]->clusterE.day
    data.frame(clusterM.day, clusterH.day, clusterE.day, unique(x_cluster$meta.data.Time))->df.cluster
    colnames(df.cluster)<-c("Mesenchymal", "Hybrid", "Epithelial", "time")
    melt(df.cluster, "time")->melted.clust
    if(!dir.exists("COMET_populated_files")){
      dir.create("COMET_populated_files")
    }
    write.csv(melted.clust, paste0("COMET_populated_files/", data.input$Sample, "_", turn, "_", cutoff, ".csv"))
    saveRDS(melted.clust, paste0("COMET_populated_files/", data.input$Sample, "_", turn, "_", cutoff, ".Rds"))
    }
  }

  return(df.cluster)
}




#' This function perform Kolmogorov Smirnov scoring on data, credit given to
#' Priyanka Chakraborty as the code was adapted from her work and modified
#' @param exp.mat gene expression matrix
#' @param genes genes
#' @param topgenes receives top 200 highly variable genes for scoring
#'
#' @return KS score
#' @importFrom stats ks.test
#' @export
KS.label.me <- function(exp.mat, genes, topgenes){
  common.sig <- intersect(my_package_env$EMT.sig[,1],genes)
  EMT.exp.idx <- match(common.sig, genes)
  exp.mat[,EMT.exp.idx]->EMT.exp
  topgenes[topgenes %in% my_package_env$EMT.sig[my_package_env$EMT.sig$...2=="Mes",][,1]]->Mes.genes
  topgenes[topgenes %in% my_package_env$EMT.sig[my_package_env$EMT.sig$...2=="Epi",][,1]]->Epi.genes
  Epi.idx <- colnames(EMT.exp[colnames(EMT.exp) %in% Epi.genes])
  Mes.idx <- colnames(EMT.exp[colnames(EMT.exp) %in% Mes.genes])


  data.frame(matrix(0, nrow <- dim(EMT.exp)[1], ncol <- 6))->ks
  for(i in 1:nrow(EMT.exp)){
    ks.test(as.numeric(EMT.exp[i,Mes.idx]),      as.numeric(EMT.exp[i,Epi.idx]))$statistic->ks$X1[i]
    ks.test(as.numeric(EMT.exp[i,Mes.idx]), as.numeric(EMT.exp[i,Epi.idx]))$p.value->ks$X2[i]
    ks.test(as.numeric(EMT.exp[i,Mes.idx]), as.numeric(EMT.exp[i,Epi.idx]), alternative="greater")$statistic->ks$X3[i]
    ks.test(as.numeric(EMT.exp[i,Mes.idx]), as.numeric(EMT.exp[i,Epi.idx]), alternative="greater")$p.value->ks$X4[i]
    ks.test(as.numeric(EMT.exp[i,Epi.idx]), as.numeric(EMT.exp[i,Mes.idx]), alternative="greater")$statistic->ks$X5[i]
    ks.test(as.numeric(EMT.exp[i,Epi.idx]), as.numeric(EMT.exp[i,Mes.idx]), alternative="greater")$p.value->ks$X6[i]
  }
  final.score <- matrix(0,nrow = nrow(ks),ncol = 1)
  for(i in 1:dim(ks)[1]){
    if(ks$X4[i]<0.05){
      final.score[i,1] <- -1*ks$X3[i]
    }else if(ks$X6[i]<0.05){
      final.score[i,1] <- ks$X5[i]
    }else{
      if(ks$X5[i] == max(c(ks$X3[i], ks$X5[i]))){
        final.score[i,1] <- max(c(ks$X3[i], ks$X5[i]))
      }else{
        final.score[i,1] <- -1 * max(c(ks$X3[i], ks$X5[i]))
      }
    }
  }
  return(final.score[,1])
}






#' This function is to be ran for all cutoffs, purpose is to find the optimal
#' number of EMT genes to minimize the DTW distance between the flow cytometry
#' trajectories and data
#' @param data.inputs the input datasheet stored in a csv file in the tables dir
#' @param tables.dir directory with a csv file that had the address
#' for the data and metadata
#' @param input.data.dir this is the directory where the input file
#' should be saved in
#'
#' @return does not return, saves the files in the COMET_populated_files dir
#' @export
generate_pipeline_files <- function(data.inputs, tables.dir, input.data.dir){

  for(cutoff in seq(5, 100, 5)){
  run_pipeline(data.inputs, tables.dir, input.data.dir, cutoff)

  }

}





#' This function calculates confidence intervals for every sample over 10
#' runs
#' @param data.inputs the input datasheet stored in a csv file in the tables dir
#'
#' @return nothing, saves results within the Confidence_Interval_Calculations dir
#' @importFrom stats setNames
#' @importFrom utils read.csv
#' @export
calculate_conf_intervals<-function(data.inputs){
  setNames(data.frame(matrix(ncol = 3, nrow = 0)),c("time","variable", "value"))->binded


  for(k in 1:nrow(data.inputs)){
    data.inputs[k,]->data.input
    for(cutoff in seq(5, 100, 5)){
      setNames(data.frame(matrix(ncol = 3, nrow = 0)),c("time","variable", "value"))->binded
      for(p in 1:10){

        read.csv(paste0("COMET_populated_files/",data.input$Sample,"_", p,"_", cutoff, ".csv"))->data
        rbind(binded, data)->binded

      }
      if(!dir.exists("Confidence_Interval_Calculations")){
        dir.create("Confidence_Interval_Calculations")
      }
      summ.binded <- Rmisc::summarySE(binded, measurevar = "value", groupvars = c("time", "variable"), conf.interval = 0.95)
      saveRDS(summ.binded, paste0("Confidence_Interval_Calculations/", data.input$Sample, "_", cutoff, ".Rds"))
    }

  }

}




#' This function calculates the DTW distance bewteen the inferred trajectories
#' for every cutoff and the flow cytometry data
#' @param data.inputs the input datasheet stored in a csv file in the tables dir
#' @param MET.range range of time for MET to take place
#'
#' @return nothing, saves the matrix in the DTW_Matrix dir
#' @importFrom dplyr %>% filter
#' @export
DTW_calculate <- function(data.inputs,  MET.range){
  for(k in 1:nrow(data.inputs)){
    data.inputs[k,]->data.input
  #initialize matrix
  matmat<-matrix(nrow=length(seq(5, 100, 5)), ncol=(3+1))
  h<-1
  for(cutoff in seq(5, 100, 5)){

    df <- readRDS(paste0("Confidence_Interval_Calculations/", data.input$Sample, "_", cutoff, ".Rds"))
    .<-time<-NULL
    .<-variable<-NULL
    #Exclude data from MET range
    df %>% filter(!(time %in% MET.range))->df
    #iterate over states and find DTW
    for(state in unique(df$variable)){

      df %>% filter(variable==state)->df.state

      dtw::dtw(df.state$value, my_package_env$flow.dat[state])$distance->dtw.d

      matmat[h, which(unique(df$variable)==state)]<-dtw.d
    }
    h+1->h
  }
  rownames(matmat)<-seq(5, 100, 5)
  #Save in matrix
  sum.mats<-0
  for(state in 1:3){

    matmat[,state]+sum.mats->sum.mats

  }
  #Find total DTW distance
  sum.mats->matmat[,(4)]
  #matmat[,1]+matmat[,2]+matmat[,3]->matmat[,4]
  if(!dir.exists("DTW_Matrix")){
    dir.create("DTW_Matrix")
  }
  saveRDS(matmat, paste0("DTW_Matrix/", data.input$Sample, "_DTW_Matrix.Rds"))
  }
}



#' This function fits optimal CTMC trajectories to timecourse data
#' @param data.inputs the input datasheet stored in a csv file in the tables dir
#' @param MET.range range of time for MET to take place
#'
#' @return the final dataframe with the optimal trajectories fitted to data
#' @export
fit.all.data <- function(data.inputs, MET.range){
  for(k in 1:nrow(data.inputs)){
    data.inputs[k,]->data.input
    find.optimal.cutoff(data.input)->opt.cutoff

    return(fit.CTMC(data.input, MET.range, opt.cutoff))
  }
}

#' Find the optimal cutoff
#' Finds the best number of highly variable EMT genes
#' @param data.input input data to use
#'
#' @return optimal cutoff of highly variable genes
#' @export
find.optimal.cutoff <- function(data.input){
  #Read the DTW matrix
  matmat <- readRDS(paste0("DTW_Matrix/", data.input$Sample, "_DTW_Matrix.Rds"))

  #Find the optimal cutoff
  opt.cutoff <- as.numeric(names(which(min(matmat[,4]) == matmat[,4])))

  return(opt.cutoff)
}


#' This function optimally fits 3 CTMC models to data (1st phase, 2nd phase,
#' MET range)
#' @param data.input input data to use
#' @param MET.range range where MET takes place
#' @param opt.cutoff optimal cutoff of highly variable genes
#'
#' @return final trajectories, lambda_E, mu, and lambda_M respectively
#' @importFrom dplyr %>% filter
#' @export
fit.CTMC <- function(data.input, MET.range, opt.cutoff){

  #Read the corresponding data for the optimal cutoff
  opt.vals <- readRDS(paste0("Confidence_Interval_Calculations/", data.input$Sample, "_", opt.cutoff, ".Rds"))
  .<-time<-NULL
  .<-variable<-NULL
  #Filter for EMT range
  opt.vals.emt <- opt.vals %>% filter(!(time %in% MET.range))

  #Find where the hybrid peaks
  opt.vals.emt %>% filter(variable == "Hybrid")->opt.vals.emt.H

  #For the first phase
  opt.vals.emt.H[which.max(opt.vals.emt.H$value),]->ind.first.phase

  #values at the point where hybrid peaks
  pi<- c((opt.vals.emt %>% filter(variable == "Epithelial") %>% filter(time ==ind.first.phase$time))$value,
  (opt.vals.emt %>% filter(variable == "Hybrid") %>% filter(time ==ind.first.phase$time))$value,
  (opt.vals.emt %>% filter(variable == "Mesenchymal") %>% filter(time ==ind.first.phase$time))$value)

  #Initial values
  E_1 <- (opt.vals.emt %>% filter(variable == "Epithelial") %>% filter(time ==min(opt.vals.emt$time)))$value
  H_1 <- (opt.vals.emt %>% filter(variable == "Hybrid") %>% filter(time ==min(opt.vals.emt$time)))$value
  M_1 <- (opt.vals.emt %>% filter(variable == "Mesenchymal") %>% filter(time ==min(opt.vals.emt$time)))$value


  p0 <- c(E_1, H_1, M_1)

  #Define these parameters for first Markov chain
  lambda_E <- pi[1]/pi[2]
  mu <- 1
  lambda_M<- pi[3]/pi[2]

  #If the epithelial value drops to zero during EMT, make sure to add a small number to avoid errors.
  if(((opt.vals.emt %>% filter(variable == "Epithelial")%>% filter(time ==max(opt.vals.emt$time))))$value==0){
    M_sc <- (((opt.vals.emt %>% filter(variable == "Epithelial")%>% filter(time ==max(opt.vals.emt$time))))$value+1e-2)/((opt.vals.emt %>% filter(variable == "Mesenchymal")%>% filter(time ==max(opt.vals.emt$time))))$value
    Mu_sc <- (((opt.vals.emt %>% filter(variable == "Epithelial")%>% filter(time ==max(opt.vals.emt$time))))$value+1e-2)/((opt.vals.emt %>% filter(variable == "Hybrid")%>% filter(time ==max(opt.vals.emt$time))))$value
  }else{
    M_sc <- (((opt.vals.emt %>% filter(variable == "Epithelial")%>% filter(time ==max(opt.vals.emt$time))))$value)/((opt.vals.emt %>% filter(variable == "Mesenchymal")%>% filter(time ==max(opt.vals.emt$time))))$value
    Mu_sc <- (((opt.vals.emt %>% filter(variable == "Epithelial")%>% filter(time ==max(opt.vals.emt$time))))$value)/((opt.vals.emt %>% filter(variable == "Hybrid")%>% filter(time ==max(opt.vals.emt$time))))$value
  }

  #Extract timepoints
  timepoints <- unique(opt.vals.emt$time)

  #Set alpha
  alph<-1
  (opt.vals %>% filter(variable == "Epithelial"))$value->E_cad
  (opt.vals %>% filter(variable == "Hybrid"))$value->hybrid
  (opt.vals %>% filter(variable == "Mesenchymal"))$value->ZEB

  #Find alpha that minimizes mse_total using the Nelder Mead algorithm
  fminres <- pracma::fminsearch(find.min.alpha, x0 = c(0.0001, 0.001),E_cad = E_cad,hybrid =  hybrid,
                                ZEB = ZEB, M_sc = M_sc, Mu_sc = Mu_sc, eq = which.max(opt.vals.emt.H$value) - 1,
                                ref_eq_day = 5, timepoints = timepoints, method="Nelder-Mead")
  alph <- fminres$xmin[1]
  #alph<-0.0145
  #Fitting first Markov chain
  E_st <- c()
  H_st <- c()
  M_st <- c()

  E_st[1]<-E_1
  H_st[1]<-H_1
  M_st[1]<-M_1

  #In case Hybrid peaks at first day, set to NULL
  if(ind.first.phase$time==0){
    E_st<-NULL
    H_st<-NULL
    M_st<-NULL
    loop_to <- length(timepoints)/my_package_env$CTMC.scale
  }else{
  for(jj in seq(1, ind.first.phase$time/my_package_env$CTMC.scale)){
    p <- transition_matrix(lambda_E, lambda_M, mu, mu, my_package_env$inc_num*(jj-1));
    vect = p0 %*% p
    E_st[jj] <- vect[,1]
    H_st[jj] <- vect[,2]
    M_st[jj] <- vect[,3]
    p0 <- vect
    loop_to <- (max(timepoints)-ind.first.phase$time)/(my_package_env$CTMC.scale*ind.first.phase$time);
  }
  }

  #Second phase CTMC
  E_f<-c()
  H_f<-c()
  M_f<-c()

  lst.f <- run.CTMC(alph, loop_to*my_package_env$CTMC.scale, 1/M_sc, 1/Mu_sc, p0)
  lst.f[[1]]-> E_f
  lst.f[[2]]-> H_f
  lst.f[[3]]-> M_f
  lst.f[[4]]->p0
  #MET range Markov chain
  #alpha_rev <- 1
  trans_lambda_M <- alph*M_sc
  trans_Mu <- alph*Mu_sc
  trans_lambda_E <- alph

  fminresrev <- pracma::fminsearch(find.min.alpha, x0 = c(1, 0.1),E_cad = E_cad,hybrid =  hybrid,
                                ZEB = ZEB, M_sc = 1/(ZEB[1]/E_cad[1]), Mu_sc = 1/(hybrid[1]/E_cad[1]), eq = 4,
                                ref_eq_day = 8, timepoints = c(max(timepoints), MET.range), method="Nelder-Mead")
  alph_rev <- fminresrev$xmin[1]

  E_end<-c()
  H_end<-c()
  M_end<-c()

  lst.end <- run.CTMC(alph_rev, length(MET.range), (ZEB[1]/E_cad[1]), (hybrid[1]/E_cad[1]), p0)
  lst.end[[1]]-> E_end
  lst.end[[2]]-> H_end
  lst.end[[3]]-> M_end

  E_final <- c(E_st, E_f[-1], E_end)
  H_final <- c(H_st, H_f[-1], H_end)
  M_final <- c(M_st, M_f[-1], M_end)

  #Step II. Fit CTMC
  time_0 <- seq(0, unique(opt.vals.emt$time)[which.max(opt.vals.emt.H$value)], length.out = unique(opt.vals.emt$time)[which.max(opt.vals.emt.H$value)]/my_package_env$CTMC.scale)
  time_1 <- seq(unique(opt.vals.emt$time)[which.max(opt.vals.emt.H$value)], max(opt.vals.emt$time), length.out = loop_to)
  time_2 <- seq(max(opt.vals.emt$time), max(MET.range), length.out =length(MET.range)/my_package_env$CTMC.scale)
  #Should add a conditional to take into account when no MET range is considered.
  time_2[-1]->time_2
  data.frame(time = c(time_0, time_1[-1], time_2), E_final, H_final, M_final)->final.df
  return(list(final.df, trans_lambda_E, trans_Mu, trans_lambda_M))
}




#' This function generates trajectories for a CTMC model given parameters
#' @param alph_fun alpha parameter used in function
#' @param time.range specific timerange to generate trajectory
#' @param M_sc_fun ratio of M/E at steady state
#' @param Mu_sc_fun ratio of H/E at steady state
#' @param p0_fun initial state vector
#'
#' @return trajectories along with the resulting p vector
#' @export
run.CTMC <- function(alph_fun, time.range, M_sc_fun, Mu_sc_fun, p0_fun){
  E_end_fun<-c()
  H_end_fun<-c()
  M_end_fun<-c()

  for (j in seq(1, time.range/my_package_env$CTMC.scale)){
    lambda_E_fun = alph_fun;
    lambda_M_fun = alph_fun*M_sc_fun;
    mu_fun = alph_fun*Mu_sc_fun;
    p_fun <- transition_matrix(lambda_E_fun, lambda_M_fun, mu_fun, mu_fun, my_package_env$inc_num*(j-1))
    vect_fun = p0_fun %*% p_fun;
    E_end_fun[j] = vect_fun[,1];
    H_end_fun[j] = vect_fun[,2];
    M_end_fun[j] = vect_fun[,3];
    p0_fun <- vect_fun;
  }
  return(list(E_end_fun, H_end_fun, M_end_fun, p0_fun))
}


#' This function finds the mse_total given a certain alpha
#' @param alpha parameter alpha
#' @param E_cad Epithelial percentage, just named E_cad
#' @param hybrid hybrid percentage
#' @param ZEB Mesenchymal percentage, just named ZEB
#' @param M_sc ratio of M/E at steady state
#' @param Mu_sc ratio of H/E at steady state
#' @param eq what timepoint to start
#' @param ref_eq_day steady state timepoint
#' @param timepoints total timepoints
#'
#' @return mse_total
#' @export
find.min.alpha <- function(alpha, E_cad, hybrid, ZEB, M_sc, Mu_sc, eq, ref_eq_day, timepoints){

  p0 <- c(E_cad[eq+1], hybrid[eq+1], ZEB[eq+1])

  E<-0
  H<-0
  M<-0

  for(k in 1:100){
    lambda_E <- alpha
    lambda_M <- alpha/M_sc
    mu <- alpha/Mu_sc

    p <- transition_matrix(lambda_E, lambda_M, mu, mu, my_package_env$inc_num*(k-1));
    vect = p0 %*% p
    E[k] = vect[,1]
    H[k] = vect[,2]
    M[k] = vect[,3]

  }

  scaling_sys <- floor(timepoints/max(timepoints)*(length(E)-1)+1)

  mse_E <- mean((E[scaling_sys]-(E_cad[(eq+1):ref_eq_day])^2))/(ref_eq_day - eq)
  mse_H <- mean((H[scaling_sys]-(hybrid[(eq+1):ref_eq_day])^2))/(ref_eq_day  - eq)
  mse_M <- mean((M[scaling_sys]-(ZEB[(eq+1):ref_eq_day])^2))/(ref_eq_day  - eq)
  mse_total = mse_E + mse_M + mse_H
  mse_total = abs(mse_total)
  if(is.nan(mse_total)){
    mse_total=Inf
  }

  print(mse_total)

}



#' This function finds the transition matrix given the parameters for the
#' generator matrix
#' @param lambda_E Rate of transition from H to E
#' @param lambda_M Rate of transition from H to M
#' @param mu_E Rate of transition from E to H
#' @param mu_M Rate of transition from M to H
#' @param t time
#'
#' @return the probability transition matrix
#' @importFrom expm expm
#' @export
transition_matrix <- function(lambda_E, lambda_M, mu_E, mu_M, t){

  Q <- t(matrix(c(-mu_E, mu_E, 0, lambda_E, -(lambda_E+lambda_M), lambda_M, 0, mu_M, -mu_M), nrow = 3, ncol = 3))

  return(expm(Q*t))
}

#   ````````````````````````````````````````````
#   /dsoooooooooooooooooooooooooooooooooooooydd-
#   /mdo///////////////////////////////////sh+d-
#   /mmmmmmmmmmmmmmmmmmmmmo:hmmmmmmmmmmmmmmm.`d-
#   /mmmmmmm+/++++++++hmmmo oy+++++++++smmmm.`d.
#   +mmmmmmm.        `ymmmo os`        /mmmm+.d.
#   /hhhhhhh.        `ymmmo os`        :hhhhhyh.
#   ````````         `ymmmo os`        ````````
#                    `ymmmo os`
#      -ooooo/`      `ymmmo os` `+ooo+`    `/oooo.
#      :hmmmdo`      `ymmmo os` .sdmmm+`   /dmmdy.
#      .ymhmd/`      `ymmms os`  `ymddd+` :dmhmd.
#     `smd:omd:      `ymmms os`  `ymy+md/:dmsomd.
#    `ommo//hmh-     `ymmms os`  `ymy.omdhmy.omd.
#    /mmdddhddmy.    `ymmms os`  `ymy``smmy.`omd.
# ``:dms-....:dms``  `ymmms os`  .ymy. .sh- `omd.
# :ydmmh.   `+dmmys` `ymmms +s` .sdmds. .- `ohmmy-
# -+++++.   `:++++/` `ymmms os` .+ooo+`    `/oooo-
#                    `ymmms os`
#                ````.hmmms os.````
#                `+hhhhdmmms odhhyy/
#                `ommmmmmmmhymmmy.y+
#                `smmmmmmmmmmmmms.y+
#                `oddddddddddddddhh/
#                ``````````````````-+







