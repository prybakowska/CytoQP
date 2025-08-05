.plot_spill <- function(data_df, custom_colors = c("white", "lightcoral", "red2", "darkred"),
                        title="Spillover Matrix", ...){
  if(!T %in%  unname(rowSums(data_df)<=100)){
    data_df <- data_df[,-unname(which(colMeans(data_df)==0))]

  }
else{
  col_selected <- unname(which(colSums(data_df)<=100))
  row_selected <- unname(which(rowSums(data_df)<=100))
  data_df <- data_df[-row_selected,-col_selected]

}


  nrow_grid = nrow(data_df)
  ncol_grid = ncol(data_df)
  scale_data <- range(data_df[data_df != 100])
  data_df <- as.data.frame(as.table(data_df))

  plot_sm <- ggplot(data_df, aes(x = Var2, y = Var1, fill = Freq)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = c(custom_colors[1:3], "grey", custom_colors[4]), limits = scale_data, guide = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5, size = 22)) +
    labs(title = title,
         x = NULL,
         y = NULL)+
    geom_hline(yintercept = seq(0.5, nrow_grid + 0.5), color = "black", size = 0.5) +
    geom_vline(xintercept = seq(0.5, ncol_grid + 0.5), color = "black", size = 0.5) +
    geom_text(aes(label = ifelse(Freq != 0 & Freq != 100, as.character(Freq), "")),
              color = "black", size = 2) +
    coord_fixed(ratio = 1)

  return(plot_sm)
}

.catalyst_compensation <- function(input,choosen_channels = NULL, isotope_list = NULL, ...){
  if(is.null(choosen_channels)){
    choosen_channels = flowCore::colnames(input)[grep("Di",
                                                      flowCore::colnames(input))]

  }
  if (is.null(isotope_list)){
    isotope_list <- .fixed_isotypes(input)
  }

  bc_key <- as.numeric(gsub("\\D", "", choosen_channels))
  bc_key <- bc_key[order(bc_key)]
  sce <- prepData(input)
  sce_bead <- assignPrelim(sce, bc_key)
  sce_bead <- applyCutoffs(estCutoffs(sce_bead))
  sce_bead <- computeSpillmat(sce_bead)
  return(metadata(sce_bead)$spillover_matrix)
}



.fixed_isotypes <- function(ff){
  isotope_list <- CATALYST::isotope_list
  if(class(ff)[1]=="flowframe"){
  ### fixation of isotope list ###
  marker_names <- flowCore::colnames(ff)[grep("Di",
                                              flowCore::colnames(ff))]
  }else{
    marker_names <- colnames(ff)
  }
  isotope_list_input <- unique(stringr::str_extract(marker_names, "[A-Za-z]+"))

  diff_isotopes <- setdiff(isotope_list_input,names(isotope_list))
  for (isotopes in diff_isotopes) {
    isotope_list[isotopes] <- as.numeric(gsub("\\D", "", marker_names[grepl(isotopes,marker_names)]))

  }


  return(isotope_list)
}


#' Plot Spillover Matrix
#'
#' @description Creates Spillover matrix and plots.
#'
#' @param input Matrix, dataframe or FCS file path.
#' @param choosen_channels You can select channels using character vector which you used in experiment (e.g. c("Cd110Di","In115Di", "Pr141Di")).
#' @param transformation Default is TRUE. If you want to use transformed data, it fixes automatically according to arcsinh.
#' @param isotope_list If you have isotope which is not including in CATALYST isotope list, you can assign but it fixes automatically.
#'
#' @return returns ggplot object
#' @export
#'
#' @examples
#' plt <- plotspillmat("231115_CompBeads_Aq1_4C_01.FCS", choosen_channels = c("Cd110Di","In115Di", "Pr141Di",
#'                                                                            "Nd142Di", "Nd143Di", "Nd144Di",
#'                                                                            "Nd145Di", "Nd146Di", "Dy164Di",
#'                                                                            "Er167Di", "Er168Di", "Yb172Di"))
#' ggsave(plt, "Spillover matrix")
#'
plotspillmat <- function(input, choosen_channels = NULL,
                         transformation =T, isotope_list = NULL,
                         savespillmat = NULL,
                         ...){
  if(is.matrix(input) | is.data.frame(input) ){
    sm <- round(input * 100, 1)
    plot_data <- .plot_spill(sm, ...)
  }
  else if (flowCore::isFCSfile(input)){
    input <- read.FCS(input, transformation = F, truncate_max_range = F)
    if(isFALSE(transformation)){
      input <- flowCore::transform(input, flowCore::transformList(flowCore::colnames(input)[grep("Di",
                                                                                            flowCore::colnames(input))], CytoNorm::cytofTransform.reverse ))

    }

    if (is.null(choosen_channels)){
        sm <- .catalyst_compensation(input, ...)
        sm <- sm[,-unname(which(colMeans(sm)==0))]
      }
    else {
        sm <- .catalyst_compensation(input,choosen_channels, ...)
    }
    if (!is.null(savespillmat)){
      write.csv(sm, savespillmat)
    }
      sm <- round(sm * 100, 1)
      plot_data <- .plot_spill(sm, ...)

      }else {
    stop("'spillover_matrix' must be matrix, DataFrame or FCS file")
  }

return(plot_data)

}





#' Data Compensation
#'
#' @description Compensates the data using CATALYST package.
#'
#' @param samples_input Flowframe, character vector with the paths of fcs files.
#' @param compensation_input Matrix, dataframe or FCS file path.
#' @param choosen_channels You can select channels using character vector which you used in experiment (e.g. c("Cd110Di","In115Di", "Pr141Di")).
#' @param transformation Default is TRUE. If you want to use transformed data, it fixes automatically according to arcsinh.(It is only for flowframe variable not fcs files.)
#' @param isotope_list If you have isotope which is not including in CATALYST isotope list, you can assign but it fixes automatically.
#'
#' @return returns ggplot object
#' @export
#'
#' @examples
#' samples <- list.files("./", pattern = ".fcs", full.names = T)
#' compensation(samples, "231115_CompBeads_Aq1_4C_01.FCS", choosen_channels = c("Cd110Di","In115Di", "Pr141Di",
#'                                                                              "Nd142Di", "Nd143Di", "Nd144Di",
#'                                                                              "Nd145Di", "Nd146Di", "Dy164Di",
#'                                                                              "Er167Di", "Er168Di", "Yb172Di"))
#'
compensation <- function(samples_input, compensation_input,
                         choosen_channels = NULL, transformation =T,
                         isotope_list = NULL, out_dir = NULL,
                         suffix_old = "cleaned", suffix_new = "compensated",...){
  if(!is.vector(samples_input) && class(samples_input) != "flowFrame"){
    stop("Please sample input must be in vector or Flowframe format.")
  }
  else{
    if(is.matrix(compensation_input) | is.data.frame(compensation_input)){
      isotope_list = .fixed_isotypes(compensation_input)
      sm = compensation_input
    }else if(flowCore::isFCSfile(compensation_input)){
      if (is.null(out_dir)){
        out_dir <- "Compensated_files"
      }
      if(!dir.exists(out_dir)){dir.create(out_dir)}
      bead_file <- read.FCS(compensation_input, transformation = F, truncate_max_range = F)
      if(is.null(isotope_list)){
        isotope_list <- .fixed_isotypes(bead_file)
      }
      sm <- .catalyst_compensation(bead_file,choosen_channels, ...)
    }else{
      stop("'spillover_matrix' must be matrix, DataFrame or fcs file")
    }
    isotope_list[["BCKG"]] <- c(190) #Added because got error that compensation matrix is not valid
    if(is.vector(samples_input)){
      for (file in samples_input) {
        sample_ff <- read.FCS(file, transformation = F, truncate_max_range = F)
        desc <- unname(sample_ff@parameters$desc)
        sample_ff@parameters$desc <- stringr::str_split_fixed(desc,"_",2)[,1]

        colnames(sample_ff@exprs) <- unname(colnames(sample_ff@exprs))

        sce_sample <- prepData(sample_ff)
        sample_comp <- CATALYST::compCytof(sce_sample, sm, overwrite = F, isotope_list = isotope_list)
        out_fcs <- sce2fcs(sample_comp, assay = "compexprs")
        out_name <- sprintf(paste0(out_dir, "/%s"),gsub(suffix_old,suffix_new,basename(file)))
        out_fcs <- flowCore::transform(out_fcs, flowCore::transformList(flowCore::colnames(out_fcs)[grep("Di",
                                                                                    flowCore::colnames(out_fcs))], CytoNorm::cytofTransform.reverse ))
        out_fcs@parameters$desc <- desc
        write.FCS(out_fcs,filename = out_name)
      }
    }else{

      if(isFALSE(transformation)){
        input <- flowCore::transform(samples_input, flowCore::transformList(flowCore::colnames(samples_input)[grep("Di",
                                                                                                   flowCore::colnames(samples_input))], CytoNorm::cytofTransform.reverse ))

      }

      sce_sample <- prepData(samples_input)
      sample_comp <- CATALYST::compCytof(sce_sample, sm, overwrite = F, isotope_list = isotope_list)
      out_fcs <- sce2fcs(sample_comp, assay = "compexprs")
      out_fcs <- flowCore::transform(out_fcs, flowCore::transformList(flowCore::colnames(out_fcs)[grep("Di",
                                                                                                       flowCore::colnames(out_fcs))], CytoNorm::cytofTransform.reverse ))
      return(out_fcs)
    }

  }

}



#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# plt <- plotspillmat("231115_CompBeads_Aq1_4C_01.FCS", choosen_channels = c("Cd110Di","In115Di", "Pr141Di",
#                                                                                               "Nd142Di", "Nd143Di", "Nd144Di",
#                                                                                               "Nd145Di", "Nd146Di", "Dy164Di",
#                                                                                               "Er167Di", "Er168Di", "Yb172Di"))
# plt
#
# ff <- read.FCS("231115_CompBeads_Aq1_4C_01.FCS", transformation = F, truncate_max_range = F)
# sm <- .catalyst_compensation(ff, choosen_channels = c("Cd110Di","In115Di", "Pr141Di",
#                                                       "Nd142Di", "Nd143Di", "Nd144Di",
#                                                       "Nd145Di", "Nd146Di", "Dy164Di",
#                                                       "Er167Di", "Er168Di", "Yb172Di"))

# k <- plotspillmat(sm, title = "test")
#
# k
#
# plt2 <- plotspillmat(sm, title="plt2")
#
# plt2
#
#
# samples <- list.files("./", pattern = ".fcs", full.names = T)
#
# compensation(samples, "231115_CompBeads_Aq1_4C_01.FCS", choosen_channels = c("Cd110Di","In115Di", "Pr141Di",
#                                                                              "Nd142Di", "Nd143Di", "Nd144Di",
#                                                                              "Nd145Di", "Nd146Di", "Dy164Di",
#                                                                              "Er167Di", "Er168Di", "Yb172Di"))
#
#
#
# test <- read.FCS("231115_CompBeads_BD1_01.fcs", transformation = F, truncate_max_range = F)
#
#
# test2 <- compensation(test, sm, choosen_channels = c("Cd110Di","In115Di", "Pr141Di",
#                                                      "Nd142Di", "Nd143Di", "Nd144Di",
#                                                      "Nd145Di", "Nd146Di", "Dy164Di",
#                                                      "Er167Di", "Er168Di", "Yb172Di"))
#
