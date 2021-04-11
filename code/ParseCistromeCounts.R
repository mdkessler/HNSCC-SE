ParseCistromeCounts <- function(df, f, shade = "10k",
                             count_min = 0){
  # function to parse counts from df from Ilya
  # This function parses the counts into a 2x2x2 matrix that is
  # ready to be provided as input to Sasha's fisher like 2x2x2 test (FLINT)
  
  # df is assumed to have column Transcription.factor,
  # which will be used as the output name for each list entry
   
  # count_min allows filtering out count tables where any
  # of the counts are higher lower than count_min
  
  # mandate shade is 10k or 100k
  stopifnot(shade == "10k" | shade == "100k")
  
  # mandate f is an epigenomic feature  == SE | E | P
  stopifnot(f == "Promoter" | f == "Enhancer" | f == "SE")
  
  # save TF vector before subsetting down
  TFs <- as.character(df$Transcription.factor)
  
  # subset df down to only cols with f in the name
  df <- df[,grepl(f, colnames(df))]
    
  # output list to store 8 value vectors
  output <- list()
  for (i in 1:nrow(df)){
    if (shade == "10k"){
      vals <- c(df[i,10],
                df[i,11],
                df[i,13],
                df[i,14],
                df[i,2],
                df[i,3],
                df[i,5],
                df[i,6])
    } else{
      vals <- c(df[i,10],
                df[i,11],
                df[i,15],
                df[i,16],
                df[i,2],
                df[i,3],
                df[i,7],
                df[i,8])
    }
    
    # filter out if any counts < count_min
    if (sum(vals < count_min) >= 1){
      next
    }
    
    # add vals to output list
    TF <- TFs[i]
    output[[TF]] <- vals
  }
  # return output list
  return(output)
}
