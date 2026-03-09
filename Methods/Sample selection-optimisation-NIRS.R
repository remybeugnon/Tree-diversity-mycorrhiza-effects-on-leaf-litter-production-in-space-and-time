rm(list = ls())
#### > 0. PACKAGES AND FUNCTIONS ####
library(tidyverse)
library(magrittr)
library(stringr)
library(purrr)
library(FactoMineR)
library(factoextra)
library(ClusterR)
np.split.Multiblock.DX <- function (file = NULL, wd = NULL, export = FALSE){
  # intro ####
  if(is.null(wd) == TRUE) wd <- getwd()
  if((("samples" %in% dir(wd))==FALSE)&(export == TRUE)) dir.create("samples")
  if (is.null(file))  stop("You must provide a file name")
  
  # start import ####
  lines <- readLines(file)
  block_pat <- "##BLOCKS=(.*)"
  bc <- grep(block_pat, lines)
  bc <- sub(block_pat, "\\1", lines[bc])
  bc <- as.integer(bc)
  st_pat <- "##TITLE="
  st <- grep(st_pat, lines)
  st <- st[-1]
  end_pat <- "##END="
  end <- grep(end_pat, lines)
  end <- end[-length(end)]
  if (length(st) != length(end)) stop("Block starts and stops did not match")
  nb <- length(st)
  if (nb != bc) stop("Block count in link block did not match the number of blocks found")
  
  # block decomposition ####
  blocks <- vector("list", nb)
  for (i in 1:nb) blocks[[i]] <- (lines[st[i]:end[i]])
  
  fnames <- sapply(st,function(X,vec,sub_ch) gsub(sub_ch, "\\1", vec[X]), vec = lines, sub_ch = st_pat)
  fnames <- str_squish(fnames) %>% str_replace_all(" ","_") %>% str_replace_all("\\+","p") %>% str_replace_all("\\-","m")
  
  if (anyDuplicated(fnames))warning("Duplicated sample names found.\n\t\tLater samples will overwrite earlier samples\n\t\tunless you edit the original multiblock file.")
  names(blocks) <- paste0("S_",fnames)
  if(export == TRUE) for (i in 1:nb) writeLines(blocks[[i]], paste0(wd,"/samples/",fnames[i],".dx"))
  invisible(blocks)
}
np.read.block <- function (jdx){
  sstt <- grep( "^\\s*##XYDATA\\s*=\\s*\\(X\\+\\+\\(Y\\.\\.Y\\)\\)$", jdx)
  send <- grep("^\\s*##END\\s*=", jdx)
  
  metadata <- jdx[1:(sstt[1]-1)]
  
  Format <- c("metadata", "XYY")
  FirstLine <- c(1, sstt)
  LastLine <- c(sstt[1] - 1, send)
  
  DF <- data.frame(Format, FirstLine, LastLine, stringsAsFactors = FALSE)
  keep_lines <- (1+DF[2,2]):(DF[2,3]-2)
  fmr <- sapply(jdx[keep_lines],str_split,pattern = " ", simplify = TRUE) %>% sapply(as.numeric) %>% t()
  rownames(fmr) <- NULL
  
  VL <- list()
  VL$Dataguide <- DF
  VL$Metadata <- metadata
  VL$spectra <-  cbind(fmr[,1],rowMeans(fmr[,-1])) %>% as.data.frame()
  colnames(VL$spectra) <- c("wl","int")
  
  fmr <- metadata[grep("YFACTOR",metadata)] %>% str_split(" ")
  VL$spectra$int <- as.numeric(fmr[[1]][2])*VL$spectra$int
  
  fmr <- (grep("CONCENTRATIONS",metadata)+1):(grep("DELTAX",metadata)-1)
  VL$concentration <- str_remove_all(metadata[fmr],"\\(") %>%
    str_remove("\\)") %>% str_split(",",simplify = TRUE)
  fmr <- paste0(VL$concentration[,1],"_",VL$concentration[,3])
  VL$concentration <- as.numeric(VL$concentration[,2])
  names(VL$concentration) <- fmr
  return(VL)
}
np.export.date <- function(spl){
  spl$name <- str_split(spl$Metadata[grep("##TITLE",spl$Metadata)],"= ",simplify = TRUE)[,2] %>% str_trim() %>% str_squish()
  spl$date <- str_split(spl$Metadata[grep("##DATE",spl$Metadata)],"= ",simplify = TRUE)[,2] %>% str_trim() %>% str_squish()
  spl$time <- str_split(spl$Metadata[grep("##TIME",spl$Metadata)],"= ",simplify = TRUE)[,2] %>% str_trim() %>% str_squish()
  return(spl)
}
np.import <- function(name = "sample.dx"){
  np.split.Multiblock.DX(file = name) %>% lapply(np.read.block) %>% lapply(np.export.date)
}

#### > 1. Data ####
df = np.import('1-data/1-2-NIRS/MyDiv_subset1.dx')
df.1 = 
  df |> 
  purrr::map_df(
    .f = ~{
      .x$spectra |> 
        dplyr::mutate(sample = paste0(.x$name, "-",.x$time))
    }
  )

df.large = 
  pivot_wider(df.1, 
              names_from = wl,
              values_from = int)


#### > 2. Sample selection #####
#### >> 2.1 PCA ####
n.tot.1 = 100
n.tot.2 = 252
library(furrr)

plan(multisession, workers = 6)

df.optim = 
  future_map(
    .x = 1:1000,
    .progress = T,
    .f = ~{
      set.seed(.x)
      # PCA measure
      pca.nirs = PCA(df.large[,-1], 
                     scale.unit = T,
                     ncp = 10, 
                     graph = F)
      
      # Variables effects
      pca.var = get_pca_var(pca.nirs)
      pca.var.c = pca.var$contrib |> 
        data.frame()
      
      pca.var.c$id = 
        unlist(rownames(pca.var.c)) 
      
      pca.var.c = pivot_longer(pca.var.c, 
                               1:10, 
                               names_to = 'dim', 
                               values_to = 'loading')
      
      #### >> 2.2 kmeans ####
      pca.ind = get_pca_ind(pca.nirs)$coord[,1:3]
      wss <- function(k) {
        kmeans(pca.ind, k, nstart = 10)$tot.withinss
      }
      k.values <- 1:15
      wss_values <- map_dbl(k.values, wss)
      
      nb.k.selected = 4
      
      k.means.nirs = kmeans(pca.ind, nb.k.selected, nstart = 10)
      
      cluster.ind = 
        pca.ind |>
        data.frame() |> 
        mutate(cluster = k.means.nirs$cluster |> 
                 unname() |> 
                 factor())
      
      #### >> 2.3 Second kmeans for 100 samples ####
      cluster.2.ind = 
        map_df(.x = 1: nb.k.selected, 
               .f = ~{
                 d = cluster.ind |> 
                   filter(cluster == .x)
                 
                 k.sub = kmeans(d[1:2], 
                                round(n.tot.1/nb.k.selected), 
                                nstart = 10)
                 c.ind = 
                   d |>
                   data.frame() |> 
                   mutate(c.2 = k.sub$cluster |> unname())
               })
      
      cluster.2.ind = 
        cluster.2.ind |>
        mutate(cc = paste0(cluster, '-', c.2))
      
      # Random selection in cluster
      df.selection.100 = 
        map_df(
          .x = unique(cluster.2.ind$cc),
          .f = ~{
            cluster.2.ind |> 
              filter(cc == .x) |> 
              sample_n(1)
          })
      
      # Joint data 
      df.100 = 
        df.large |>
        select(sample) |> 
        add_column(as.data.frame(pca.ind)) |>
        right_join(df.selection.100 |> 
                     select(-c.2), 
                   by = c('Dim.1', 'Dim.2', 'Dim.3'))
      
      #### >> 2.4 selections 250 ####
      cluster.3.ind = 
        map_df(.x = 1: nb.k.selected, 
               .f = ~{
                 d = cluster.ind |> 
                   filter(cluster == .x)
                 
                 k.sub = kmeans(d[1:2], 
                                round(n.tot.2/nb.k.selected), 
                                nstart = 10)
                 c.ind = 
                   d |>
                   data.frame() |> 
                   mutate(c.3 = k.sub$cluster |> unname())
               })
      
      cluster.3.ind = 
        cluster.3.ind |>
        mutate(cc = paste0(cluster, '-', c.3))
      
      # Random selection in cluster
      df.selection.250 = 
        map_df(
          .x = unique(cluster.3.ind$cc),
          .f = ~{
            d = cluster.3.ind |> 
              filter(cc == .x)
            d.j = inner_join(d, df.100 |>
                               select('Dim.1', 'Dim.2', 'Dim.3'), 
                             by = c('Dim.1', 'Dim.2', 'Dim.3'))
            if(nrow(d.j)>0){d.j}else{sample_n(d,1)}
          })
      
      # Joint data 
      df.250 = 
        df.large |>
        select(sample) |> 
        add_column(as.data.frame(pca.ind)) |>
        right_join(df.selection.250 |> 
                     select(-c.3), 
                   by = c('Dim.1', 'Dim.2', 'Dim.3'))
      
      #### > 3. Selection control ####
      #### >> 3.1 Control match selection 100 and 250 ####
      # sum(df.100$sample %in% df.250$sample)
      
      #### >> 3.2 Check groups balance ####
      
      # Should not be selected: '15-AE2-Ca-Mar'
      a = unique(df.large$sample)
      a = a[grepl('15-AE2-Ca-Mar', a)]
      a %in% df.100$sample
      a %in% df.250$sample
      
      df.large$species = 
        df.large$sample |> 
        str_split('-') |> 
        map_chr(.f = ~{.x[3]})
      
      df.large$plot = 
        df.large$sample |> 
        str_split('-') |> 
        map_chr(.f = ~{paste(.x[1],'-',.x[2])})
      
      df.large$month = 
        df.large$sample |> 
        str_split('-') |> 
        map_chr(.f = ~{.x[4]})
      
      df.100$species = 
        df.100$sample |> 
        str_split('-') |> 
        map_chr(.f = ~{.x[3]})
      
      df.100$plot = 
        df.100$sample |> 
        str_split('-') |> 
        map_chr(.f = ~{paste(.x[1],'-',.x[2])})
      
      df.100$month = 
        df.100$sample |> 
        str_split('-') |> 
        map_chr(.f = ~{.x[4]})
      
      df.design = readxl::read_xlsx(path ='1-data/1-0-design.xlsx') |> 
        select(plot = plotName, TSR = Tree_sp_rich, myc = Myc_type) |> 
        mutate(treatment = paste0(TSR, '-', myc)) |> 
        mutate(plot = str_replace_all(plot, '\\.', ' - ')) |> 
        group_by(plot, treatment) |> 
        summarise()
      
      df.100 = df.100 |> 
        left_join(df.design, by = 'plot')
      
      df.large = df.large |> 
        left_join(df.design, by = 'plot')
      
      dist.species.tot = 
        df.large |> 
        group_by(species) |> 
        summarise(n = n())
      
      dist.plot.tot = 
        df.large |> 
        group_by(treatment) |> 
        summarise(n = n())
      
      dist.month.tot = 
        df.large |> 
        group_by(month) |> 
        summarise(n = n())
      
      dist.species = 
        df.100 |> 
        group_by(species) |> 
        summarise(n = n())
      
      dist.plot = 
        df.100 |> 
        group_by(treatment) |> 
        summarise(n = n())
      
      dist.month = 
        df.100 |> 
        group_by(month) |> 
        summarise(n = n())
      
      L = list(
        dist.m = left_join(dist.month, dist.month.tot, by = 'month') |> mutate(d = n.x - n.y),
        dist.s = left_join(dist.species, dist.species.tot, by = 'species') |> mutate(d = n.x - n.y),
        dist.p = left_join(dist.plot, dist.plot.tot, by = 'treatment') |> mutate(d = n.x - n.y)
      )
      c(ini = .x, 
        to.rem = sum(a %in% df.100$sample, a %in% df.250$sample),
        m = mean(L$dist.m$d),
        s = mean(L$dist.s$d),
        p = mean(L$dist.p$d))
    }
  ) |> 
  map_df(.f = ~.x)

df.optim.1 = df.optim |> 
  mutate(m.s = (m-mean(m))/sd(m)) |> 
  mutate(s.s = (s-mean(s))/sd(s)) |> 
  mutate(p.s = (p-mean(p))/sd(p))

df.optim.1$val = rowSums(abs(df.optim.1[,c('m.s','s.s','p.s')]), na.rm = T)

write_csv(df.optim.1, '1-data/1-2-NIRS/optim.csv')

df.optim.1 = read_csv('1-data/1-2-NIRS/optim.csv')
df.optim.1 = df.optim.1 |> filter(to.rem == 0)
selected.ini = df.optim.1[df.optim.1$val == min(df.optim.1$val), 1]
selected.ini
df.optim.1[selected.ini$ini,] |> data.frame()

