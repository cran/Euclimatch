
#' Climatch Summary Score
#'
#' A summarized climatch score within the recipient region to a source region. Provides the percentage of climate data points in the recipient region equal to or above a specified score (default is 6), or the mean climatch score across the whole recipient region. Note no floor function is used as in Crombie et al. (2008)
#'
#' @param recipient A data.frame or list of data.frames of climatic variables for the recipient region.
#' @param source A data.frame or list of data.frames of climatic variables for the source region.
#' @param globvar A vector of the global variance of each climate variable.
#' @param type Specifies the type of summary score to use. "perc" (default) specifies a percent climatch score representing the number of grid cells above or equal to a given value specified with the 'thershold' argument. "mean" provides the mean climatch score across the recipient region.
#' @param threshold The climatch score to use in calculating the percentage match, which is the number of grid cells within the recipient region with a climatch >= the threshold (default is 6).
#'
#' @return A numeric value, vector or data.frame of the percentage of climatch scores within the recipient region(s) >= a threshold value, or the mean climatch score across the region(s).
#'
#' @usage climatch_sum(recipient, source, globvar, type = "perc", threshold = 6)
#'
#' @references Predicting invasiveness of species in trade: climate match, trophic guild and fecundity influence establishment and impact of non-native freshwater fishes"<doi:10.1111/ddi.12391>
#'
#' @examples
#' i <- as.data.frame(matrix(runif(n=180, min=1, max=20), nrow=60)) # Fake source climate data
#' j <- as.data.frame(matrix(runif(n=300, min=10, max=40), nrow=100)) # Fake recipient data
#' variance <- c(600, 800, 450) # Fake global variance
#'
#' climatch_sum(recipient = j, source = i, globvar = variance, type = "perc", threshold = 6)
#'
#' @export

climatch_sum <- function(recipient, source, globvar, type = "perc", threshold = 6){

  # If recipient is list of data.frames and source is a single dataframe
  if(class(recipient)[[1]]=="list" & class(source)[[1]]=="data.frame"){
    # Run the climatch algorithm for the vector of climatch scores
    match_vec <- numeric()
    for(i in 1:length(recipient)){
      match <- climatch_vec(recipient = recipient[[i]], source = source, globvar = globvar)
      # Summarize the climatch vector as a percent match above the threshold
      if(type == "perc"){
        match_vec[i] <- (sum(match >= threshold) / nrow(recipient[[i]])) * 100
      }

      # Summarize the climatch vector as the mean score
      if(type == "mean"){
        match_vec[i] <- mean(match)
      }
    }
    return(match_vec)
  }
  #If recipient is a data.frame and source a list
  if(class(recipient)[[1]]=="data.frame" & class(source)[[1]]=="list"){
    # Run the climatch algorithm for the vector of climatch scores
    match_vec <- numeric()
    for(i in 1:length(source)){
      match <- climatch_vec(recipient = recipient, source = source[[i]], globvar = globvar)
      # Summarize the climatch vector as a percent match above the threshold
      if(type == "perc"){
        match_vec[i] <- (sum(match >= threshold) / nrow(recipient)) * 100
      }

      # Summarize the climatch vector as the mean score
      if(type == "mean"){
        match_vec[i] <- mean(match)
      }
    }
    return(match_vec)
  }

  # If both recipient and sources are lists
  if(class(recipient)[[1]]=="list" & class(source)[[1]]=="list"){
    # Run the climatch algorithm for the vector of climatch scores

    match_df <- data.frame()
    for(j in 1:length(recipient)){
      for(i in 1:length(source)){
        match_j <- climatch_vec(recipient = recipient[[j]], source = source[[i]], globvar = globvar)
        # Summarize the climatch vector as a percent match above the threshold
        if(type == "perc"){
          match_df[j,i] <- (sum(match_j >= threshold) / nrow(recipient[[j]])) * 100
        }

        if(type == "mean"){
          match_df[j,i] <- mean(match_j)
        }
      }
    }
    return(match_df)
  }

  if(class(recipient)[[1]]=="data.frame" & class(source)[[1]]=="data.frame"){
  # Run the climatch algorithm for the vector of climatch scores
  match_vec <- climatch_vec(recipient = recipient, source = source, globvar = globvar)

  # Summarize the climatch vector as a percent match above the threshold
  if(type == "perc"){
    match_score <- (sum(match_vec >= threshold) / nrow(recipient)) * 100
  }

  # Summarize the climatch vector as the mean score
  if(type == "mean"){
    match_score <- mean(match_vec)
  }


  return(match_score)
  }
}

