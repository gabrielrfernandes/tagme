#' Assign taxonomy to an input FASTA file
#'
#' @param file The FASTA file to be classified. It can be a group of ASV, OTU or Unique sequences.
#' @param db The directory containing the downloaded RDS and TXT files corresponding to your amplicon region. The directory name must end with "/".
#' @param specificity The Specificity that will determine the threshold score for assignment.
#' @param batch The maximum number of sequences to be assigned per batch. Bigger batches demand more memory.
#'
#' @return A table containing the taxonomic assignment and the score to each sequence.
#'
#' @export
tagmeFromFasta <- function(file, db = "./", specificity = 0.8, batch = 50000){

  require("randomForest")
  require("seqinr")

  kmer_size=4

  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  if (substrRight(db, 1) != "/") {
    stop("The path to the DATABASE directory must end with \"/\".\n", call. = FALSE)
  }

  if(!file.exists(file)){
    stop(paste(file, "not found."))
  }

  # creating a sequence object from fasta file
  seqobj <- seqinr::read.fasta(file = file, forceDNAtolower = T, set.att = FALSE)

  # Generating all kmer possibilities to format the header.
  # I am passing a toy sequence just to avoid heavy calculations.
  kmers<-seqinr::count("a", kmer_size)
  kmers<-names(kmers)

  # We consider the frequency of a kmer: kmer frequency plus the frequency of the
  # reverse complement. So we will split all the kmers into: [1] Main kmers;
  # [2] Kmer reverse complement.
  nonred_main<-NULL # will store the main kmers
  nonred_rev<-NULL # will store the reverse complement of the main kmers
  hash_nonred <- new.env(hash = TRUE) # hash to control the kmers and reverse complements
  index<-1


  for(i in 1:length(kmers)){
    inner_tetra<-kmers[i]
    compl<-c2s(rev(comp(s2c(inner_tetra))))
    # if the kmer is not in the array we put it and set the reverse complement
    if(is.null(hash_nonred[[inner_tetra]])){
      nonred_main[index]<-inner_tetra
      nonred_rev[index]<-compl
      # we set both the kmer and the reverse in the hash to avoid duplication
      hash_nonred[[inner_tetra]]<-1
      hash_nonred[[compl]]<-1
      index<-index+1
    }
  }

  # For matter of performance we create the whole matrix with zero because is expected to have a sparse
  # matrix
  final_matrix <- matrix(0, nrow = length(seqobj), ncol = length(nonred_main))

  # the colnames will be each kmer in upper case
  colnames(final_matrix)<-toupper(nonred_main)
  # the row names will be the ids of the sequence
  nomes = names(seqobj)
  # For each sequence we calculate the kmer frequency
  for(si in seq_along(seqobj)){
    sequence <- seqobj[[si]] # the sequence itself
    seqid <- names(seqobj)[si] # the sequence identifier
    # to calculate the frequencies we call the function count of the package seqinr
    frequencies <- seqinr::count(sequence, kmer_size)

    # for each kmer in the array nonred_main we sum its frequency with the reverse
    # complement
    for(i in 1:length(nonred_main)){
      kmer_freq <- frequencies[nonred_main[i]]
      # if the kmer and its reverse are identical means that its a palindrome
      if(!identical(nonred_main[i], nonred_rev[i])){
        kmer_freq <- kmer_freq + frequencies[nonred_rev[i]]
      }
      final_matrix[si,i]<-kmer_freq
    }
  }
  data.file = data.frame(final_matrix, seqID=nomes)

  fullsize = dim(data.file)[1]
  times = fullsize%/%batch
  remain = fullsize%%batch

  if (dim(data.file)[1] > batch) {

    for (r in 0:times) {
      if (r == times) {
        start = start+batch
        end = end+remain
      } else if (r==0) {
        start = 1
        end = batch
      } else {
        start = start+batch
        end = end+batch
      }

      cat (r,": ",start,"-",end,"\n\n")

      data = data.file[start:end,]

      cat("Starting Species... ")

      table = read.table(paste(db,"specie.txt",sep = ""), header = T, row.names = 1, sep = "\t")
      cutoff = rownames(table)[which.min(abs(table$Specificity/specificity -1))]

      model = readRDS(paste(db,"species2.rds", sep = ""))
      pred = predict(model, newdata=data)
      pred.prob = predict(model, newdata=data, type="prob")
      second = function (x) {
        x[order(x, decreasing=T)[2]]
      }
      result = data.frame(seqID=data$seqID, Taxonomy=pred, Best=apply(pred.prob, 1, max), Second=apply(pred.prob, 1, second), Ratio = log2(apply(pred.prob, 1, max)/apply(pred.prob, 1, second))*apply(pred.prob, 1, max))
      result$Ratio[result$Ratio=="Inf"] = 10
      result.final = result[which(result$Ratio>cutoff),]

      Rest = data.frame(seqID = result$seqID[which(result$Ratio<=cutoff)])
      data = unique(merge(Rest, data))
      model = readRDS(paste(db,"species.uniq.rds", sep = ""))
      pred = predict(model, newdata=data)
      pred.prob = predict(model, newdata=data, type="prob")
      result = data.frame(seqID=data$seqID, Taxonomy=pred, Best=apply(pred.prob, 1, max), Second=apply(pred.prob, 1, second), Ratio = log2(apply(pred.prob, 1, max)/apply(pred.prob, 1, second))*apply(pred.prob, 1, max))
      result$Ratio[result$Ratio=="Inf"] = 10
      result.final = rbind(result.final, result[which(result$Ratio>cutoff),])

      number = dim(result.final)[1]
      cat(number, " assigned\n\n")

      cat("Starting Genus... ")

      Rest = data.frame(seqID = result$seqID[which(result$Ratio<=cutoff)])

      table = read.table(paste(db,"genus.txt", sep=""), header = T, row.names = 1, sep = "\t")
      cutoff = rownames(table)[which.min(abs(table$Specificity/specificity -1))]

      data = unique(merge(Rest, data))
      model = readRDS(paste(db,"genus2.rds", sep=""))
      pred = predict(model, newdata=data)
      pred.prob = predict(model, newdata=data, type="prob")
      result = data.frame(seqID=data$seqID, Taxonomy=pred, Best=apply(pred.prob, 1, max), Second=apply(pred.prob, 1, second), Ratio = log2(apply(pred.prob, 1, max)/apply(pred.prob, 1, second))*apply(pred.prob, 1, max))
      result$Ratio[result$Ratio=="Inf"] = 10
      result.final = rbind(result.final, result[which(result$Ratio>cutoff),])

      number = dim(result.final)[1]
      cat(number, " assigned\n\n")

      cat("Starting Family... ")

      Rest = data.frame(seqID = result$seqID[which(result$Ratio<=(cutoff))])

      table = read.table(paste(db,"family.txt", sep=""), header = T, row.names = 1, sep = "\t")
      cutoff = rownames(table)[which.min(abs(table$Specificity/specificity -1))]

      data = unique(merge(Rest, data))
      model = readRDS(paste(db,"family.rds", sep=""))
      pred = predict(model, newdata=data)
      pred.prob = predict(model, newdata=data, type="prob")
      result = data.frame(seqID=data$seqID, Taxonomy=pred, Best=apply(pred.prob, 1, max), Second=apply(pred.prob, 1, second), Ratio = log2(apply(pred.prob, 1, max)/apply(pred.prob, 1, second))*apply(pred.prob, 1, max))
      result$Ratio[result$Ratio=="Inf"] = 10
      result.final = rbind(result.final, result[which(result$Ratio>(cutoff)),])

      number = dim(result.final)[1]
      cat(number, " assigned\n\n")

      cat("Starting Order... ")

      Rest = data.frame(seqID = result$seqID[which(result$Ratio<=(cutoff))])

      table = read.table(paste(db,"order.txt", sep = ""),header = T, row.names = 1, sep = "\t")
      cutoff = rownames(table)[which.min(abs(table$Specificity/specificity -1))]

      data = unique(merge(Rest, data))
      model = readRDS(paste(db,"order.rds", sep = ""))
      pred = predict(model, newdata=data)
      pred.prob = predict(model, newdata=data, type="prob")
      result = data.frame(seqID=data$seqID, Taxonomy=pred, Best=apply(pred.prob, 1, max), Second=apply(pred.prob, 1, second), Ratio = log2(apply(pred.prob, 1, max)/apply(pred.prob, 1, second))*apply(pred.prob, 1, max))
      result$Ratio[result$Ratio=="Inf"] = 10
      result.final = rbind(result.final, result[which(result$Ratio>(cutoff)),])

      number = dim(result.final)[1]
      cat(number, " assigned\n\n")

      cat("Starting Class... ")

      Rest = data.frame(seqID = result$seqID[which(result$Ratio<=(cutoff))])

      table = read.table(paste(db,"class.txt", sep=""), header = T, row.names = 1, sep = "\t")
      cutoff = rownames(table)[which.min(abs(table$Specificity/specificity -1))]

      data = unique(merge(Rest, data))
      model = readRDS(paste(db,"class.rds", sep = ""))
      pred = predict(model, newdata=data)
      pred.prob = predict(model, newdata=data, type="prob")
      result = data.frame(seqID=data$seqID, Taxonomy=pred, Best=apply(pred.prob, 1, max), Second=apply(pred.prob, 1, second), Ratio = log2(apply(pred.prob, 1, max)/apply(pred.prob, 1, second))*apply(pred.prob, 1, max))
      result$Ratio[result$Ratio=="Inf"] = 10
      result.final = rbind(result.final, result[which(result$Ratio>(cutoff)),])

      number = dim(result.final)[1]
      cat(number, " assigned\n\n")

      cat("Starting Phylum... ")

      Rest = data.frame(seqID = result$seqID[which(result$Ratio<=(cutoff))])

      table = read.table(paste(db,"phylum.txt",sep=""), header = T, row.names = 1, sep = "\t")
      cutoff = rownames(table)[which.min(abs(table$Specificity/specificity -1))]

      data = unique(merge(Rest, data))
      model = readRDS(paste(db,"phylum.rds", sep = ""))
      pred = predict(model, newdata=data)
      pred.prob = predict(model, newdata=data, type="prob")
      result = data.frame(seqID=data$seqID, Taxonomy=pred, Best=apply(pred.prob, 1, max), Second=apply(pred.prob, 1, second), Ratio = log2(apply(pred.prob, 1, max)/apply(pred.prob, 1, second))*apply(pred.prob, 1, max))
      result$Ratio[result$Ratio=="Inf"] = 10
      result.final = rbind(result.final, result[which(result$Ratio>(cutoff)),])

      number = dim(result.final)[1]
      cat(number, " assigned\n\n")

      cat("Starting Domain... ")

      Rest = data.frame(seqID = result$seqID[which(result$Ratio<=(cutoff))])

      data = unique(merge(Rest, data))
      model = readRDS(paste(db,"domain.rds", sep = ""))
      pred = predict(model, newdata=data)
      pred.prob = predict(model, newdata=data, type="prob")
      result = data.frame(seqID=data$seqID, Taxonomy=pred, Best=apply(pred.prob, 1, max), Second=apply(pred.prob, 1, second), Ratio = log2(apply(pred.prob, 1, max)/apply(pred.prob, 1, second))*apply(pred.prob, 1, max))
      result$Ratio[result$Ratio=="Inf"] = 10
      result.final = rbind(result.final, result[which(result$Ratio>(1)),])

      number = dim(result.final)[1]
      cat(number, " assigned\n\n")

      cat("Printing Unassigned... \n\n")

      Rest = data.frame(seqID = result$seqID[which(result$Ratio<=(1))])
      if (dim(Rest)[1]>0) {
        data = unique(merge(Rest, data))
        pred = predict(model, newdata=data)
        pred.prob = predict(model, newdata=data, type="prob")
        result = data.frame(seqID=data$seqID, Taxonomy="Unassigned", Best=apply(pred.prob, 1, max), Second=apply(pred.prob, 1, second), Ratio = log2(apply(pred.prob, 1, max)/apply(pred.prob, 1, second))*apply(pred.prob, 1, max))
        result.final = rbind(result.final, result)
      }
      if (r == 0) {
        export = result.final
      } else if (r == times) {
        export = rbind(export, result.final)
        return(export)
      } else {
        export = rbind(export, result.final)
      }
      r = r+1
    }

  } else {

    data = data.file

    cat("Starting Species... ")

    table = read.table(paste(db,"specie.txt",sep = ""), header = T, row.names = 1, sep = "\t")
    cutoff = rownames(table)[which.min(abs(table$Specificity/specificity -1))]

    model = readRDS(paste(db,"species2.rds", sep = ""))
    pred = predict(model, newdata=data)
    pred.prob = predict(model, newdata=data, type="prob")
    second = function (x) {
      x[order(x, decreasing=T)[2]]
    }
    result = data.frame(seqID=data$seqID, Taxonomy=pred, Best=apply(pred.prob, 1, max), Second=apply(pred.prob, 1, second), Ratio = log2(apply(pred.prob, 1, max)/apply(pred.prob, 1, second))*apply(pred.prob, 1, max))
    result$Ratio[result$Ratio=="Inf"] = 10
    result.final = result[which(result$Ratio>cutoff),]

    Rest = data.frame(seqID = result$seqID[which(result$Ratio<=cutoff)])
    data = unique(merge(Rest, data))
    model = readRDS(paste(db,"species.uniq.rds", sep = ""))
    pred = predict(model, newdata=data)
    pred.prob = predict(model, newdata=data, type="prob")
    result = data.frame(seqID=data$seqID, Taxonomy=pred, Best=apply(pred.prob, 1, max), Second=apply(pred.prob, 1, second), Ratio = log2(apply(pred.prob, 1, max)/apply(pred.prob, 1, second))*apply(pred.prob, 1, max))
    result$Ratio[result$Ratio=="Inf"] = 10
    result.final = rbind(result.final, result[which(result$Ratio>cutoff),])

    number = dim(result.final)[1]
    cat(number, " assigned\n\n")

    cat("Starting Genus... ")

    Rest = data.frame(seqID = result$seqID[which(result$Ratio<=cutoff)])

    table = read.table(paste(db,"genus.txt", sep=""), header = T, row.names = 1, sep = "\t")
    cutoff = rownames(table)[which.min(abs(table$Specificity/specificity -1))]

    data = unique(merge(Rest, data))
    model = readRDS(paste(db,"genus2.rds", sep=""))
    pred = predict(model, newdata=data)
    pred.prob = predict(model, newdata=data, type="prob")
    result = data.frame(seqID=data$seqID, Taxonomy=pred, Best=apply(pred.prob, 1, max), Second=apply(pred.prob, 1, second), Ratio = log2(apply(pred.prob, 1, max)/apply(pred.prob, 1, second))*apply(pred.prob, 1, max))
    result$Ratio[result$Ratio=="Inf"] = 10
    result.final = rbind(result.final, result[which(result$Ratio>cutoff),])

    if(file.exists("genus.uniq.rds")){
      Rest = data.frame(seqID = result$seqID[which(result$Ratio<=cutoff)])
      data = unique(merge(Rest, data))
      model = readRDS(paste(db,"genus.uniq.rds", sep=""))
      pred = predict(model, newdata=data)
      pred.prob = predict(model, newdata=data, type="prob")
      result = data.frame(seqID=data$seqID, Taxonomy=pred, Best=apply(pred.prob, 1, max), Second=apply(pred.prob, 1, second), Ratio = log2(apply(pred.prob, 1, max)/apply(pred.prob, 1, second))*apply(pred.prob, 1, max))
      result$Ratio[result$Ratio=="Inf"] = 10
      result.final = rbind(result.final,result[which(result$Ratio>(cutoff)),])
    }
    number = dim(result.final)[1]
    cat(number, " assigned\n\n")

    cat("Starting Family... ")

    Rest = data.frame(seqID = result$seqID[which(result$Ratio<=(cutoff))])

    table = read.table(paste(db,"family.txt", sep=""), header = T, row.names = 1, sep = "\t")
    cutoff = rownames(table)[which.min(abs(table$Specificity/specificity -1))]

    data = unique(merge(Rest, data))
    model = readRDS(paste(db,"family.rds", sep=""))
    pred = predict(model, newdata=data)
    pred.prob = predict(model, newdata=data, type="prob")
    result = data.frame(seqID=data$seqID, Taxonomy=pred, Best=apply(pred.prob, 1, max), Second=apply(pred.prob, 1, second), Ratio = log2(apply(pred.prob, 1, max)/apply(pred.prob, 1, second))*apply(pred.prob, 1, max))
    result$Ratio[result$Ratio=="Inf"] = 10
    result.final = rbind(result.final, result[which(result$Ratio>(cutoff)),])

    number = dim(result.final)[1]
    cat(number, " assigned\n\n")

    cat("Starting Order... ")

    Rest = data.frame(seqID = result$seqID[which(result$Ratio<=(cutoff))])

    table = read.table(paste(db,"order.txt", sep = ""),header = T, row.names = 1, sep = "\t")
    cutoff = rownames(table)[which.min(abs(table$Specificity/specificity -1))]

    data = unique(merge(Rest, data))
    model = readRDS(paste(db,"order.rds", sep = ""))
    pred = predict(model, newdata=data)
    pred.prob = predict(model, newdata=data, type="prob")
    result = data.frame(seqID=data$seqID, Taxonomy=pred, Best=apply(pred.prob, 1, max), Second=apply(pred.prob, 1, second), Ratio = log2(apply(pred.prob, 1, max)/apply(pred.prob, 1, second))*apply(pred.prob, 1, max))
    result$Ratio[result$Ratio=="Inf"] = 10
    result.final = rbind(result.final, result[which(result$Ratio>(cutoff)),])

    number = dim(result.final)[1]
    cat(number, " assigned\n\n")

    cat("Starting Class... ")

    Rest = data.frame(seqID = result$seqID[which(result$Ratio<=(cutoff))])

    table = read.table(paste(db,"class.txt", sep=""), header = T, row.names = 1, sep = "\t")
    cutoff = rownames(table)[which.min(abs(table$Specificity/specificity -1))]

    data = unique(merge(Rest, data))
    model = readRDS(paste(db,"class.rds", sep = ""))
    pred = predict(model, newdata=data)
    pred.prob = predict(model, newdata=data, type="prob")
    result = data.frame(seqID=data$seqID, Taxonomy=pred, Best=apply(pred.prob, 1, max), Second=apply(pred.prob, 1, second), Ratio = log2(apply(pred.prob, 1, max)/apply(pred.prob, 1, second))*apply(pred.prob, 1, max))
    result$Ratio[result$Ratio=="Inf"] = 10
    result.final = rbind(result.final, result[which(result$Ratio>(cutoff)),])

    number = dim(result.final)[1]
    cat(number, " assigned\n\n")

    cat("Starting Phylum... ")

    Rest = data.frame(seqID = result$seqID[which(result$Ratio<=(cutoff))])

    table = read.table(paste(db,"phylum.txt",sep=""), header = T, row.names = 1, sep = "\t")
    cutoff = rownames(table)[which.min(abs(table$Specificity/specificity -1))]

    data = unique(merge(Rest, data))
    model = readRDS(paste(db,"phylum.rds", sep = ""))
    pred = predict(model, newdata=data)
    pred.prob = predict(model, newdata=data, type="prob")
    result = data.frame(seqID=data$seqID, Taxonomy=pred, Best=apply(pred.prob, 1, max), Second=apply(pred.prob, 1, second), Ratio = log2(apply(pred.prob, 1, max)/apply(pred.prob, 1, second))*apply(pred.prob, 1, max))
    result$Ratio[result$Ratio=="Inf"] = 10
    result.final = rbind(result.final, result[which(result$Ratio>(cutoff)),])

    number = dim(result.final)[1]
    cat(number, " assigned\n\n")

    cat("Starting Domain... ")

    Rest = data.frame(seqID = result$seqID[which(result$Ratio<=(cutoff))])

    data = unique(merge(Rest, data))
    model = readRDS(paste(db,"domain.rds", sep = ""))
    pred = predict(model, newdata=data)
    pred.prob = predict(model, newdata=data, type="prob")
    result = data.frame(seqID=data$seqID, Taxonomy=pred, Best=apply(pred.prob, 1, max), Second=apply(pred.prob, 1, second), Ratio = log2(apply(pred.prob, 1, max)/apply(pred.prob, 1, second))*apply(pred.prob, 1, max))
    result$Ratio[result$Ratio=="Inf"] = 10
    result.final = rbind(result.final, result[which(result$Ratio>(1)),])

    number = dim(result.final)[1]
    cat(number, " assigned\n\n")

    cat("Printing Unassigned... \n\n")

    Rest = data.frame(seqID = result$seqID[which(result$Ratio<=(1))])
    if (dim(Rest)[1]>0) {
      data = unique(merge(Rest, data))
      pred = predict(model, newdata=data)
      pred.prob = predict(model, newdata=data, type="prob")
      result = data.frame(seqID=data$seqID, Taxonomy="Unassigned", Best=apply(pred.prob, 1, max), Second=apply(pred.prob, 1, second), Ratio = log2(apply(pred.prob, 1, max)/apply(pred.prob, 1, second))*apply(pred.prob, 1, max))
      result.final = rbind(result.final, result)
    }
    return(result.final)
  }
}

