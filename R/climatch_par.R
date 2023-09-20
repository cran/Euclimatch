#' Run climatch in parallel
#'
#' @param recipient List of data.frames of the recipient regions
#' @param source List of dataf.rames of the source regions
#' @param biovar Vector of the columns (climate variables) to use, default all columns
#' @param globvar Vector of the global variance of each variable
#' @param ncores The number of cores to use in parallel
#' @param type Choose between "perc" (default) or "mean" passed to climatch_sum() and "vec" passes to climatch_vec()
#' @param threshold The climatch score (0-10) to use in calculating the percentage match, which is the number of grid cells within the recipient region with a climatch >= the threshold (default is 6).
#' @return "perc" and "mean" returns data.frame of climatch within recipients (rows) for each source represented in columns, "vec" returns data.frame of climatch of a recipient (each column corresponds to grid cell), to sources (corresponding to rows)
#'
#' @importFrom foreach %:% %dopar%
#' @import doParallel
#' @import RcppParallel
#'
#' @usage climatch_par(recipient,
#' @usage              source,
#' @usage              globvar,
#' @usage              biovar = 1:length(globvar),
#' @usage              ncores,
#' @usage              type = "perc",
#' @usage              threshold = 6)
#'
#' @examples
#' # Dummy data
#' i1 <- data.frame("clim1" = 1:10, "clim2" = 9:18) # Fake source climate data
#' i <- list(i1, i1) # list the source dataframes
#' j1 <- data.frame("clim1" = 11:20, "clim2" = 16:25) # Fake recipient climate data
#' j <- list(j1, j1) # list the recipient dataframes
#' variance <- c(60, 80) # Fake global variance
#'
#' # Climate matching
#' climatch_par(recipient = j, source = i, globvar = variance, ncores = 1, type = "vec")
#' @export
climatch_par <- function(recipient, source, globvar, biovar = 1:length(globvar), ncores = 1, type = "perc", threshold = 6) {

  if(class(recipient)[[1]]=="data.frame"){
    recipient <- list(recipient)
  }

  if(class(source)[[1]]=="data.frame"){
    source <- list(source)
  }

  recipient_n <- length(recipient)
  source_n <- length(source)

  climatch_pairwise <- data.frame(nrow = recipient_n, ncol = source_n)

  # Set i and j as global variables
  i <- NULL
  j <- NULL

  # Set up parallel processing
  nc <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(nc)

  if(type == 'perc'){
    # Parallel foreach loop - each row will represent recipient regions and columns for sources
    climatch_pairwise <- foreach::foreach(i = 1:source_n, .combine = "cbind", .inorder = TRUE, .packages = c("Euclimatch")) %:%
      foreach::foreach(j = 1:recipient_n, .combine = "c", .inorder = TRUE, .packages = c("Euclimatch")) %dopar% {
        Euclimatch::climatch_sum(recipient = recipient[[j]][, biovar, drop = FALSE], source = source[[i]][, biovar, drop = FALSE], globvar = globvar, type = type, threshold = threshold)
      }
  }

  if(type == 'mean'){
    # Parallel foreach loop - each row will represent recipient regions and columns for sources
    climatch_pairwise <- foreach::foreach(i = 1:source_n, .combine = "cbind", .inorder = TRUE, .packages = c("Euclimatch")) %:%
      foreach::foreach(j = 1:recipient_n, .combine = "c", .inorder = TRUE, .packages = c("Euclimatch")) %dopar% {
        Euclimatch::climatch_sum(recipient = recipient[[j]][, biovar, drop = FALSE], source = source[[i]][, biovar, drop = FALSE], globvar = globvar, type = type)
      }
  }

  if(type == 'vec'){
    recipient_df <- recipient[[1]]
    recipient_n <- nrow(recipient_df)
    source_n <- length(source)
    climatch_pairwise <- data.frame(nrow = source_n, ncol = recipient_df)
    # Single foreach loop where each is climatch vector is a row representing match to each source, each column in the grid cell or data point in the recipient region
    climatch_pairwise <- foreach::foreach(i = 1:source_n, .combine = "rbind", .inorder = TRUE, .packages = c("Euclimatch")) %dopar% {
      Euclimatch::climatch_vec(recipient = recipient_df[, biovar, drop = FALSE], source = source[[i]][, biovar, drop = FALSE], globvar = globvar)
    }
  }

  parallel::stopCluster(nc)  # Stop cluster
  return(climatch_pairwise)
}

