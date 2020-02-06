# For all proteins, do pairwise comparison based on the group variance estimated from linear mixed model
# data: long-format input data with columns: RUN, LABEL, SUBJECT, GROUP, ABUNDANCE, FEATURE
# MS: data type, which can be "FG.NormalizedMS1PeakArea", "FG.NormalizedMS2PeakArea", c("FG.NormalizedMS1PeakArea", "FG.NormalizedMS2PeakArea")
# summarization: summarize precursor intensity to protein intensity by summarization method. Default is `median`. The other option is `mean`.
# adj.method: methods for adjusting the pvalue. Default is `BH``
# Return a data frame whose columns are Protein, Compare, nFeatures, log2FC, SE, DF, pvalue, adjusted.pvalue. Compare indicates the comparison between which two groups.nFeatures is the number of precursors belonging to each protein
protein.linear <- function(data, MS, summarization = "median", adj.method = "BH") {
  
  # Input checking
  if(length(MS) == 1){
    if(!(MS == "FG.NormalizedMS1PeakArea" |
       MS == "FG.NormalizedMS2PeakArea")){
      stop("MS must be one of \"FG.NormalizedMS1PeakArea\", \"FG.NormalizedMS2PeakArea\" or 
           c(\"FG.NormalizedMS1PeakArea\", \"FG.NormalizedMS2PeakArea\")!")
    }
  } else{
    if(length(MS) == 2){
      if(!all(MS == c("FG.NormalizedMS1PeakArea", "FG.NormalizedMS2PeakArea"))){
        stop("MS must be one of \"FG.NormalizedMS1PeakArea\", \"FG.NormalizedMS2PeakArea\" or 
             c(\"FG.NormalizedMS1PeakArea\", \"FG.NormalizedMS2PeakArea\")!")
      }
    } else{
      stop("MS must be one of \"FG.NormalizedMS1PeakArea\", \"FG.NormalizedMS2PeakArea\" or 
           c(\"FG.NormalizedMS1PeakArea\", \"FG.NormalizedMS2PeakArea\")!")
    }
  } 
  
  
  # do normalization for each precursor across all the runs across MS1 and MS2
  if(length(MS) == 2){
    feature.median <- data[ ,.(meidan = median(ABUNDANCE)), by = .(PROTEIN, FEATURE, LABEL)]
    data <- merge(data, feature.median, all.x=TRUE)
    data$ABUNDANCE <- data$ABUNDANCE - data$meidan
    data[ ,meidan:=NULL]
  }
  
  # Extract the required data 
  data <- data[LABEL %in% MS]
  
  # Summarize peptide intensities to protein intensity by mean or median
  if(summarization == "median"){
    prot.int <- data[, .(PROT.INT = median(ABUNDANCE, na.rm = T)), by=.(LABEL, GROUP, SUBJECT, RUN, PROTEIN)]
  } else{
    if(summarization == "mean"){
      prot.int <- data[, .(PROT.INT = mean(ABUNDANCE, na.rm = T)), by=.(LABEL, GROUP, SUBJECT, RUN, PROTEIN)]
    } else{
      stop("summarization must be through median or mean!", call.=FALSE)
    }
  }
  
  proteins <- as.character(unique(data$PROTEIN))
  groups <- unique(data$GROUP)
  nlabel <- length(unique(data$LABEL))
  pvalue <- NULL 
  log2FC <- NULL
  SE <- NULL
  DF <- NULL
  prots <- NULL
  comparison <- NULL
  
  if(length(MS) > 1){ # both MS1 and MS2
    for(i in 1:length(proteins)) {
      message("Protein ", proteins[i], ": ",  i)
      # Select one protein
      sub <- prot.int[PROTEIN == proteins[i]]
      
      # Retired code: Two-way anova model
      # model <- try(lm(PROT.INT ~ 1 + GROUP + LABEL, data=sub), TRUE)
      #if(!inherits(model, "try-error")){
        #group.size <- sub[ , .(df = .N ) , by = .(GROUP)]
        ### Get estimated fold change from mixed model
        #coeff <- coef(model)
        #group.coeff <- coeff[grep("GROUP", names(coeff))]
        #label.coeff <- coeff[grep("LABEL", names(coeff))]
        #baseline <- (nlabel*coeff["(Intercept)"] + sum(label.coeff))/nlabel # calculate group mean for the first group
        #group.coeff <- group.coeff + baseline
        #group.coeff <- c(baseline, group.coeff)
        ### Find the group name for baseline
        #names(group.coeff) <- gsub("GROUP", "", names(group.coeff))
        #names(group.coeff)[names(group.coeff)=="(Intercept)"] <- setdiff(as.character(groups), names(group.coeff))
        ### Estimate the group variance (MSE) from linear model
        #k = length(model$coefficients) - 1 #Subtract one to ignore intercept
        #n = length(model$residuals) # total df
        #df = n-(1+k) # residual df
        #SSE = sum(model$residuals**2)
        #MSE = SSE/(n-(1+k))
      
      ## if there is only one run in the data, then train one-way anova
      rand.model <- try(lmer(PROT.INT ~ 1 + GROUP + LABEL + (1 | GROUP:RUN), data=sub), TRUE)
      
      if(!inherits(rand.model, "try-error")){  
        group.size <- sub[ , .(df = .N ) , by = .(GROUP)]
        ## Get estimated fold change from mixed model
        coeff <- fixed.effects(rand.model)
        group.coeff <- coeff[grep("GROUP", names(coeff))]
        label.coeff <- coeff[grep("LABEL", names(coeff))]
        baseline <- (nlabel*coeff["(Intercept)"] + sum(label.coeff))/nlabel # calculate group mean for the first group
        group.coeff <- group.coeff + baseline
        group.coeff <- c(baseline, group.coeff)
        ## Find the group name for baseline
        names(group.coeff) <- gsub("GROUP", "", names(group.coeff))
        names(group.coeff)[names(group.coeff)=="(Intercept)"] <- setdiff(as.character(groups), names(group.coeff))
        
        # Estimate the group variance from fixed model
        fixed.model <- lm(PROT.INT ~ 1 + GROUP + LABEL + GROUP:RUN, data=sub)
        av <- anova(fixed.model)
        MSE <- av$'Mean Sq'[3]
        df <- av$Df[3]
        
        for(j in 1:(length(groups)-1)){
          for(k in (j+1):length(groups)){
            g1 <- groups[j]
            g2 <- groups[k]
            prots <- c(prots, proteins[i])
            comparison <- c(comparison, paste(g1, g2, sep="-"))
            
            ## test whether the protein abundance are signigicantly different between two geiven groups
            # Estimate fold change
            diff<- group.coeff[g1] - group.coeff[g2]
            # Estimate group variance
            variance<-MSE*sum(1/group.size[GROUP==g1,df] + 1/group.size[GROUP==g2, df])
            # Calculate t statistic
            t<-diff/sqrt(variance)
            # Calculate p-value
            pv <- 2*pt(-abs(t),df = df)
            
            pvalue <- c(pvalue, pv)
            log2FC <- c(log2FC, diff)
            SE <- c(SE, sqrt(variance))
            DF <- c(DF, df)
          }
        } # end for pairwise comparison
        
      } else { # if unable to fit the linear model
        for(j in 1:(length(groups)-1)){
          for(k in (j+1):length(groups)){
            # Assign NA to the testing results
            g1 <- groups[j]
            g2 <- groups[k]
            prots <- c(prots, proteins[i])
            comparison <- c(comparison, paste(g1, g2, sep="-"))
            pvalue <- c(pvalue, NA)
            log2FC <- c(log2FC, NA)
            SE <- c(SE, NA)
            DF <- c(DF, NA)
          }
        } # end for pairwise comparison
        
      } # end if model error for one protein
      
    } # end for all the proteins
    
  } else { # single MS1 or MS2 data
    for(i in 1:length(proteins)) {
      # message("Protein ", proteins[i], ": ",  i)
      # Select one protein
      sub <- prot.int[PROTEIN == proteins[i]]
      # Train linear model
      model <- try(lm(PROT.INT ~ 1 + GROUP, data=sub), TRUE)
      if(!inherits(model, "try-error")){
        group.size <- sub[ , .(df = .N ) , by = .(GROUP)]
        ## Get estimated fold change from mixed model
        coeff <- coef(model)
        group.coeff <- coeff[grep("GROUP", names(coeff))]
        baseline <- coeff["(Intercept)"] # calculate group mean for the first group
        group.coeff <- group.coeff + baseline
        group.coeff <- c(baseline, group.coeff)
        ## Find the group name for baseline
        names(group.coeff) <- gsub("GROUP", "", names(group.coeff))
        names(group.coeff)[names(group.coeff)=="(Intercept)"] <- setdiff(as.character(groups), names(group.coeff))
        ## Estimate the group variance (MSE) from linear model
        k = length(model$coefficients) - 1 #Subtract one to ignore intercept
        n = length(model$residuals) # total df
        df = n-(1+k) # residual df
        SSE = sum(model$residuals**2)
        MSE = SSE/(n-(1+k))
        
        message("Protein ", proteins[i], ": ",  i)
        for(j in 1:(length(groups)-1)){
          for(k in (j+1):length(groups)){
            g1 <- groups[j]
            g2 <- groups[k]
            # message(proteins[i])
            prots <- c(prots, proteins[i])
            # message(prots)
            comparison <- c(comparison, paste(g1, g2, sep="-"))
            
            ## test whether the protein abundance are signigicantly different between two geiven groups
            # Estimate fold change
            diff<- group.coeff[g1] - group.coeff[g2]
            # Estimate group variance
            variance<-MSE*sum(1/group.size[GROUP==g1,df] + 1/group.size[GROUP==g2, df])
            # Calculate t statistic
            t<-diff/sqrt(variance)
            # Calculate p-value
            pv <- 2*pt(-abs(t),df = df)
            
            pvalue <- c(pvalue, pv)
            log2FC <- c(log2FC, diff)
            SE <- c(SE, sqrt(variance))
            DF <- c(DF, df)
          }
        } # end for pairwise comparison
        
      } else { # if unable to fit the linear model
        for(j in 1:(length(groups)-1)){
          for(k in (j+1):length(groups)){
            # Assign NA to the testing results
            g1 <- groups[j]
            g2 <- groups[k]
            prots <- c(prots, proteins[i])
            comparison <- c(comparison, paste(g1, g2, sep="-"))
            pvalue <- c(pvalue, NA)
            log2FC <- c(log2FC, NA)
            SE <- c(SE, NA)
            DF <- c(DF, NA)
          }
        } # end for pairwise comparison
        
      } # end if model error for one protein
      
    } # end for all the proteins
  }
  
  # Save testing results
  statistic <- data.frame(Protein = prots, 
                          Compare = comparison, 
                          log2FC = log2FC, 
                          SE = SE, 
                          DF = DF, 
                          pvalue = pvalue)
  
  statistic$adjusted.pvalue <- p.adjust(statistic$pvalue, adj.method)
  
  #filename <- paste(paste("protein.linear", paste0(gsub("FG.Normalized", "", MS), collapse="."), sep="."), ".rda", sep = "")
  #filename <- gsub("PeakArea", "", filename)
  #save(statistic, file = filename)
  
  #comparisons <- unique(statistic$Compare)
  # Pairwise comparison
  #for(i in 1:length(comparisons)){
    #comp = comparisons[i]
    #statistic[statistic$Compare == comp, "adjusted.pvalue"] <- p.adjust(statistic[statistic$Compare == comp, "adjusted.pvalue"], adj.method)
  #}
  
  return(statistic)
}


## Check whether the input have all the requried columns
.check.input <- function(data){
  
  required.columns <- c("Run", "R.Condition", "R.Replicate", 
                        "PG.Qvalue", "EG.PrecursorId", "EG.Qvalue", 
                        "FG.MS2PeakArea", "FG.NormalizedMS2PeakArea", 
                        "FG.MS1PeakArea", "FG.NormalizedMS1PeakArea", 
                        "IDPicker.InferenceId")
  
  if (!all(required.columns %in% colnames(data))) {
    
    missedColumns <- which(!(required.columns %in% colnames(data)))
    stop(paste("Please check the required column in the input file. ** columns :",
               paste(required.columns[missedColumns], collapse = ", "), " are missed."))
    
  }
}