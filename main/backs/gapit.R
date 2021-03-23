#!/usr/bin/Rscript


main <- function () {
	#Step 1: Set data directory and import files
	genoNumFile = "example-genotype-tetra-gwaspoly-ACGT-GAPIT-NUM.csv"
	genoMapFile = "example-genotype-tetra-gwaspoly-ACGT-GAPIT-MAP.csv"
	phenoFile   = "example-phenotype-single-trait-GAPIT.csv"

	genoNumFile = "x-genotype.csv"
	genoMapFile = "x-map.csv"
	phenoFile   = "x-phenotype.csv"

	myY  = read.table(phenoFile, head = TRUE)
	myGM = read.table(genoMapFile , head = TRUE)
	myGD = read.table(genoNumFile, head = TRUE)

	#Step 2: Run GAPIT
	myGAPIT <- GAPIT(Y=myY, GD=myGD, GM=myGM,model="GLM", file.output=F)
	write.csv (myGAPIT$GWAS, "gapit-gwas.csv", quote=F, row.names=F)
}
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
getOperatingSystem <- function()
{
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
  os <- sysinf['sysname']
  if (os == 'Darwin')
    os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(as.character (os))
}


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.Memory` <- function(Memory = NULL, Infor) 
{
    #Object: To report memory usage
    #Output: Memory
    #Authors: Zhiwu Zhang
    # Last update: June 6, 2011
    ##############################################################################################
    gc()
	if (getOperatingSystem() %in% c("linux", "osx"))
		size = strsplit(system('grep MemTotal /proc/meminfo', intern = TRUE),split="\\s+")[[1]][2]
	else
    	size <- memory.size()

    #print(paste("Memory usage: ",size," for", Infor))
    if (is.null(Memory)) {
      Increased = 0
      Memory = cbind(Infor, size , Increased)
    } else{
      Increased = 0
      Memory.current = cbind(Infor, size , Increased)
      Memory = rbind(Memory, Memory.current)
      Memory[nrow(Memory), 3] = as.numeric(as.matrix(Memory[nrow(Memory), 2])) -
        as.numeric(as.matrix(Memory[nrow(Memory) - 1, 2]))
    }
    
    return (Memory)
  }#end of GAPIT.Memory function


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.Timmer` <- function(Timmer = NULL, Infor) 
{
    #Object: To report current time
    #Output: Timmer
    #Authors: Zhiwu Zhang
    # Last update: may 8, 2011
    ##############################################################################################
    Time <- Sys.time()
    if (is.null(Timmer)) {
      Elapsed = 0
      Timmer = cbind(Infor, Time, Elapsed)
    } else{
      Elapsed = 0
      Timmer.current = cbind(Infor, Time, Elapsed)
      Timmer = rbind(Timmer, Timmer.current)
      Timmer[nrow(Timmer), 3] = as.numeric(as.matrix(Timmer[nrow(Timmer), 2])) -
        as.numeric(as.matrix(Timmer[nrow(Timmer) - 1, 2]))
    }
    
    #print(paste('Time used: ', Timmer[nrow(Timmer),3], ' seconds for ',Infor,sep="" ))
    return (Timmer)
  }#end of GAPIT.Timmer function

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.Genotype` <-
  function(G = NULL,
           GD = NULL,
           GM = NULL,
           KI = NULL,
           kinship.algorithm = "Zhang",
           SNP.effect = "Add",
           SNP.impute = "Middle",
           PCA.total = 0,
           PCA.col = NULL,
           PCA.3d = PCA.3d,
           seed = 123,
           SNP.fraction = 1,
           file.path = NULL,
           file.from = NULL,
           file.to = NULL,
           file.total = NULL,
           file.fragment = 1000,
           SNP.test = TRUE,
           file.G = NULL,
           file.Ext.G = NULL,
           file.GD = NULL,
           file.Ext.GD = NULL,
           file.GM = NULL,
           file.Ext.GM = NULL,
           SNP.MAF = 0.05,
           FDR.Rate = 0.05,
           SNP.FDR = 1,
           Timmer = NULL,
           Memory = NULL,
           LD.chromosome = NULL,
           LD.location = NULL,
           LD.range = NULL,
           SNP.CV = NULL,
           GP = NULL,
           GK = NULL,
           GTindex = NULL,
           bin.size = 1000,
           inclosure.size = 100,
           sangwich.top = NULL,
           sangwich.bottom = NULL,
           file.output = TRUE,
           kinship.cluster = "average",
           NJtree.group = NULL,
           NJtree.type = c("fan", "unrooted"),
           Create.indicator = FALSE,
           Major.allele.zero = FALSE,
           Geno.View.output = TRUE) {
    #Object: To unify genotype and calculate kinship and PC if required:
    #       1.For G data, convert it to GD and GI
    #       2.For GD and GM data, nothing change
    #       3.Samling GD and create KI and PC
    #       4.Go through multiple files
    #       5.In any case, GD must be returned (for QC)
    #Output: GD, GI, GT, KI and PC
    #Authors: Zhiwu Zhang
    #Last update: August 11, 2011
    ##############################################################################################
    
    #print("Genotyping: numericalization, sampling kinship, PCs and much more...")
    
    
    
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Genotype start")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "Genotype start")
    compress_z = NULL
    type_col = NULL
    #Create logical variables
    byData = !is.null(G) | !is.null(GD)
    byFile = !is.null(file.G) | !is.null(file.GD)
    hasGenotype = (byData | byFile)
    needKinPC = (is.null(KI) |
                   PCA.total > 0 | kinship.algorithm == "Separation")
    
    if (!is.null(KI) &
        !byData & !byFile & !SNP.test & kinship.algorithm != "SUPER")
    {
      print("It return unexpected")
      return (
        list(
          GD = NULL,
          GI = NULL,
          GT = NULL,
          hasGenotype = FALSE,
          genoFormat = NULL,
          KI = KI,
          PC = NULL,
          byFile = FALSE,
          fullGD = TRUE,
          Timmer = Timmer,
          Memory = Memory
        )
      )
    }
    
    
    #Set indicator for full GD
    fullGD = FALSE
    if (byData)
      fullGD = TRUE
    if (byFile & SNP.fraction == 1 & needKinPC)
      fullGD = TRUE
    
    #SET GT to NULL in case of no genotype
    if (!byData & !byFile & is.null(GK) & kinship.algorithm != "SUPER")
    {
      if (is.null(KI) &
          is.null(GP) &
          is.null(GK))
        stop("GAPIT says: Kinship has to be provided or estimated from genotype!!!")
      return (
        list(
          GD = NULL,
          GI = NULL,
          GT = NULL,
          hasGenotype = FALSE,
          genoFormat = NULL,
          KI = KI,
          PC = NULL,
          byFile = FALSE,
          fullGD = TRUE,
          Timmer = Timmer,
          Memory = Memory
        )
      )
    }
    
    genoFormat = "hapmap"
    if (is.null(G) & is.null(file.G))
      genoFormat = "EMMA"
    
    #Multiple genotype files
    #In one of the 3 situations, calculate KI with the algorithm specified, otherwise skip cit by setting algorithm to "SUPER"
    kinship.algorithm.save = kinship.algorithm
    kinship.algorithm = "SUPER"
    #Normal
    if (is.null(sangwich.top) &
        is.null(sangwich.bottom))
      kinship.algorithm = kinship.algorithm.save
    #TOP or Bottom is MLM
    pass.top = FALSE
    if (!is.null(sangwich.top))
      pass.top = !(sangwich.top == "FaST" |
                     sangwich.top == "SUPER" | sangwich.top == "DC")
    pass.bottom = FALSE
    if (!is.null(sangwich.bottom))
      pass.bottom = !(sangwich.bottom == "FaST" |
                        sangwich.bottom == "SUPER" | sangwich.bottom == "DC")
    if (pass.top | pass.bottom)
      kinship.algorithm = kinship.algorithm.save
    #Compatibility of input
    
    #agreement among file from, to and total
    if (!is.null(file.from) & !is.null(file.to) & !is.null(file.total))
    {
      if (file.total != (file.to - file.from + 1))
        stop("GAPIT says: Conflict among file (from, to and total)")
    }
    if (!is.null(file.from) & !is.null(file.to))
    {
      if (file.to < file.from)
        stop("GAPIT says: file.from should smaller than file.to")
    }
    #file.from and file.to must be in pair
    if (is.null(file.from) &
        !is.null(file.to))
      stop("GAPIT says: file.from and file.to must be in pair)")
    if (!is.null(file.from) &
        is.null(file.to))
      stop("GAPIT says: file.from and file.to must be in pair)")
    
    #assign file.total
    if (!is.null(file.from) &
        !is.null(file.to))
      file.total = file.to - file.from + 1
    if (byFile &
        is.null(file.total))
      stop("GAPIT says: file.from and file.to must be provided!)")
    
    if (!is.null(GP) &
        !is.null(GK))
      stop("GAPIT Says: You can not provide GP and GK at same time")
    if (!is.null(GP) &
        !is.null(KI))
      stop("GAPIT Says: You can not provide GP and KI at same time")
    if (!is.null(GK) &
        !is.null(KI))
      stop("GAPIT says: You can not specify GK and KI at same time!!!")
    
    #GP does not allow TOP
    if (!is.null(GP) &
        !is.null(sangwich.top))
      stop("GAPIT Says: You provided GP. You can not spycify sangwich.top")
    
    #Top require a bottom
    if (!is.null(sangwich.top) &
        is.null(sangwich.bottom))
      stop("GAPIT Says: Top require its Bottom")
    
    #naked bottom require GP or GK
    if (is.null(sangwich.top) &
        !is.null(sangwich.bottom) &
        (is.null(GP) &
         is.null(GK)))
      stop("GAPIT Says: Uncovered Bottom (without TOP) requires GP or GK")
    
    #Pseudo top (GK or GP) requires a bottom
    if (is.null(sangwich.top) &
        is.null(sangwich.bottom) &
        (!is.null(GP) |
         !is.null(GK)))
      stop("GAPIT Says: You have provide GP or GK, you need to provide Bottom")
    
    #if(!is.null(KI) &!is.null(kinship.algorithm))  stop("GAPIT says: You can not specify kinship.algorithm and provide kinship at same time!!!")
    
    
    
    if (!needKinPC &
        SNP.fraction < 1)
      stop(
        "GAPIT says: You did not require calculate kinship or PCs. SNP.fraction should not be specified!!!"
      )
    if (!SNP.test &
        is.null(KI) &
        !byData &
        !byFile)
      stop("GAPIT says: For SNP.test optioin, please input either use KI or use genotype")
    
    #if(is.null(file.path) & !byData & byFile) stop("GAPIT Ssays: A path for genotype data should be provided!")
    if (is.null(file.total) &
        !byData &
        byFile)
      stop("GAPIT Ssays: Number of file should be provided: >=1")
    if (!is.null(G) &
        !is.null(GD))
      stop("GAPIT Ssays: Both hapmap and EMMA format exist, choose one only.")
    
    if (!is.null(file.GD) &
        is.null(file.GM) &
        (!is.null(GP) |
         !is.null(GK)))
      stop("GAPIT Ssays: Genotype data and map files should be in pair")
    if (is.null(file.GD) &
        !is.null(file.GM) &
        (!is.null(GP) |
         !is.null(GK)))
      stop("GAPIT Ssays: Genotype data and map files should be in pair")
    
    if (!is.null(GD) &
        is.null(GM) &
        (is.null(GP) &
         is.null(GK)) &
        kinship.algorithm != "SUPER")
      stop("GAPIT Says: Genotype data and map files should be in pair")
    if (is.null(GD) &
        !is.null(GM) &
        (is.null(GP) &
         is.null(GK)) &
        kinship.algorithm != "SUPER")
      stop("GAPIT Says: Genotype data and map files should be in pair")
    
    
    #if(!byData & !byFile) stop("APIT Ssays: Either genotype data or files should be given!")
    #if(byData&(!is.null(file.path))) stop ("APIT Ssays: You have provided geotype data. file.path should not be provided!")
    
    #print("Pass compatibility of input")
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Genotype loaded")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "Genotype loaded")
    
    #Inital GLD
    GLD = NULL
    SNP.QTN = NULL #Intitial
    GT = NULL
    
    #Handler of read data in numeric format (EMMA)
    #Rename GM as GI
    if (!is.null(GM))
      GI = GM
    rm(GM)
    gc()
    #Extract GD and GT from read data GD
    if (!is.null(GD))
    {
      GT = as.matrix(GD[, 1])  #get taxa
      GD = as.matrix(GD[, -1]) #remove taxa column
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "GT created from GD)")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "GT created from GD")
    }
    
    #Hapmap format
    if (!is.null(G))
    {
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Before HapMap")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "Before HapMap")
      #Convert HapMap to numerical
      print(paste("Converting genotype...", sep = ""))
      hm = GAPIT.HapMap(
        G,
        SNP.effect = SNP.effect,
        SNP.impute = SNP.impute,
        Create.indicator = Create.indicator,
        Major.allele.zero = Major.allele.zero
      )
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "after HapMap")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "after HapMap")
      #Extracting SNP for LD plot
      if (!is.null(LD.chromosome))
      {
        #print("Extracting SNP for LD plot...")
        chromosome = (G[, 3] == LD.chromosome[1])
        bp = as.numeric(as.vector(G[, 4]))
        deviation = abs(bp - as.numeric(as.vector(LD.location[1])))
        location = deviation < as.numeric(as.vector(LD.range[1]))
        index = chromosome & location
        GLD = G[index, ]
      } else{
        #print("No data in GLD")
        GLD = NULL
      }
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "HapMap")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "HapMap")
      print(paste("Converting genotype done.", sep = ""))
      #rm(G)
      #gc()
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "G removed")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "G removed")
      GT = hm$GT
      GD = hm$GD
      GI = hm$GI
      #
      #print(unique(GI[,2]))
      rm(hm)
      gc()
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "hm removed")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "hm removed")
    }
    
    #From files
    if (!byData & byFile)
    {
      #print("Loading genotype from files...")
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "byFile")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "byFile")
      numFileUsed = file.to
      if (!needKinPC)
        numFileUsed = file.from
      #Initial GI as storage
      GD = NULL
      GT = NULL
      GI = NULL
      GLD = NULL
      #multiple fragments or files
      for (file in file.from:numFileUsed)
      {
        frag = 1
        numSNP = file.fragment
        myFRG = NULL
        #print(paste("numSNP  before while is ",numSNP))
        while (numSNP == file.fragment)
        {
          #this is problematic if the read end at the last line
          print(paste("Reading file: ", file, "Fragment: ", frag))
          Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Before Fragment")
          Memory = GAPIT.Memory(Memory = Memory, Infor = "Before Fragment")
          myFRG = GAPIT.Fragment(
            file.path = file.path,
            file.from = file.from,
            file.to = file.to,
            file.total = file.total,
            file.G = file.G,
            file.Ext.G = file.Ext.G,
            seed = seed,
            SNP.fraction = SNP.fraction,
            SNP.effect = SNP.effect,
            SNP.impute = SNP.impute,
            genoFormat = genoFormat,
            file.GD = file.GD,
            file.Ext.GD = file.Ext.GD,
            file.GM = file.GM,
            file.Ext.GM = file.Ext.GM,
            file.fragment = file.fragment,
            file = file,
            frag = frag,
            LD.chromosome = LD.chromosome,
            LD.location = LD.location,
            LD.range = LD.range,
            Create.indicator = Create.indicator,
            Major.allele.zero = Major.allele.zero
          )
          Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "After Fragment")
          Memory = GAPIT.Memory(Memory = Memory, Infor = "After Fragment")
          
          if (is.null(GT) &
              !is.null(myFRG$GT))
            GT = as.matrix(myFRG$GT)
          
          if (is.null(GD))
          {
            GD = myFRG$GD
          } else{
            if (!is.null(myFRG$GD))
            {
              GD = cbind(GD, myFRG$GD)
            }
          }
          if (is.null(GI))
          {
            GI = myFRG$GI
          } else{
            if (!is.null(myFRG$GI))
            {
              colnames(myFRG$GI) = c("SNP", "Chromosome", "Position")
              GI = as.data.frame(rbind(as.matrix(GI), as.matrix(myFRG$GI)))
            }
          }
          
          if (is.null(G))
          {
            G = myFRG$G
          } else{
            if (!is.null(myFRG$G))
            {
              G = as.data.frame(rbind(as.matrix(G), as.matrix(myFRG$G[-1, ])))
            }
          }
          
          if (is.null(GLD))
          {
            GLD = myFRG$GLD
          } else{
            if (!is.null(myFRG$GLD))
            {
              if (myFRG$heading)
              {
                GLD = as.data.frame(rbind(as.matrix(GLD), as.matrix(myFRG$GLD[-1, ])))
              } else{
                GLD = as.data.frame(rbind(as.matrix(GLD), as.matrix(myFRG$GLD)))
              }
            }
          }
          
          if (file == file.from & frag == 1)
            GT = as.matrix(myFRG$GT)
          frag = frag + 1
          if (!is.null(myFRG$GI))
          {
            numSNP = myFRG$linesRead[1]
          } else{
            numSNP = 0
          }
          
          if (!needKinPC)
            numSNP = 0  #force to end the while loop
          if (is.null(myFRG))
            numSNP = 0  #force to end the while loop
          
          Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "END this Fragment")
          Memory = GAPIT.Memory(Memory = Memory, Infor = "END this Fragment")
          
          
          
        } #end whileof repeat on fragment
        # print("This file is OK")
      } #end of file loop
      print("All files loaded")
    } #end of if(!byData&byFile)
    
    #GM=as.matrix(GI)
    #GI=GM
    GM = GI
    
    # modified by Jiabo in 20190927. sorted number of chrom by numeric and charicter
    
    chor_taxa = as.character(unique(GM[, 2]))
    
    chor_taxa[order(gsub("([A-Z]+)([0-9]+)", "\\1", chor_taxa),
                    as.numeric(gsub("([A-Z]+)([0-9]+)", "\\2", chor_taxa)))]
    chr_letter = grep("[A-Z]|[a-z]", chor_taxa)
    if (!setequal(integer(0), chr_letter))
    {
      GI = as.matrix(GI)
      for (i in 1:(length(chor_taxa)))
      {
        index = GM[, 2] == chor_taxa[i]
        GI[index, 2] = i
      }
    }
    
    #print(chor_taxa)
    #print(head(GI))
    #print("@@@@@@@@@@@")
    #print(GD[1:5,1:5])
    #print(dim(GI))
    #Follow the MAF to filter markers
    if (!is.null(GD))
    {
      #maf=apply(as.matrix(GD),2,function(one) abs(1-sum(one)/(2*nrow(GD))))
      #maf[maf>0.5]=1-maf[maf>0.5]
      ss = apply(GD, 2, sum)
      maf = apply(cbind(.5 * ss / (nrow(GD)), 1 - .5 * ss / (nrow(GD))), 1, min)
      #print(max(maf))
      #print(min(maf))
      maf_index = maf >= SNP.MAF
      print(paste("GAPIT will filter marker with MAF setting !!"))
      print(paste("The markers will be filtered by SNP.MAF: ", SNP.MAF, sep =
                    ""))
      print(table(maf_index))
      
      #print(head(maf[!maf_index]))
      
      GD = GD[, maf_index]
      GI = as.data.frame(GI[maf_index, ])
      GM = as.data.frame(GM[maf_index, ])
      #GI=GM
    }
    #print("file loaded")
    
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Sampling genotype")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "Sampling genotype")
    #print(KI)
    #Plot third part kinship
    if (!is.null(KI) & file.output)
    {
      if (KI != 1)
      {
        if (nrow(KI) < 2000)
        {
          print("Plotting Kinship")
          #print(dim(KI))
          theKin = as.matrix(KI[, -1])
          line.names <- KI[, 1]
          colnames(theKin) = KI[, 1]
          rownames(theKin) = KI[, 1]
          distance.matrix = dist(theKin, upper = TRUE)
          hc = hclust(distance.matrix, method = kinship.cluster)
          hcd = as.dendrogram(hc)
          ##plot NJtree
          if (!is.null(NJtree.group))
          {
            clusMember <- cutree(hc, k = NJtree.group)
            compress_z = table(clusMember, paste(line.names))
            type_col = rainbow(NJtree.group)
            Optimum = c(nrow(theKin), kinship.cluster, NJtree.group)
          }
          Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "set kinship")
          Memory = GAPIT.Memory(Memory = Memory, Infor = "set kinship")
          if (file.output)
          {
            print("Creating heat map for kinship...")
            pdf(
              paste("GAPIT.Kin.thirdPart.pdf", sep = ""),
              width = 12,
              height = 12
            )
            par(mar = c(25, 25, 25, 25))
            Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "prepare heatmap")
            Memory = GAPIT.Memory(Memory = Memory, Infor = "prepare heatmap")
            heatmap.2(
              theKin,
              cexRow = .2,
              cexCol = 0.2,
              col = rev(heat.colors(256)),
              scale = "none",
              symkey = FALSE,
              trace = "none"
            )
            dev.off()
            print("Kinship heat map PDF created!")
            Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "plot heatmap")
            Memory = GAPIT.Memory(Memory = Memory, Infor = "plot heatmap")
          }
          ## Jiabo Wang add NJ Tree of kinship at 4.5.2017
          if (!is.null(NJtree.group) & file.output)
          {
            for (tr in 1:length(NJtree.type))
            {
              print("Creating NJ Tree for kinship...")
              pdf(
                paste("GAPIT.Kin.NJtree.", NJtree.type[tr], ".pdf", sep = ""),
                width = 12,
                height = 12
              )
              par(mar = c(5, 5, 5, 5))
              Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "prepare NJ TREE")
              Memory = GAPIT.Memory(Memory = Memory, Infor = "prepare NJ TREE")
              plot(
                as.phylo(hc),
                type = NJtree.type[tr],
                tip.color = type_col[clusMember],
                use.edge.length = TRUE,
                col = "gray80",
                cex = 0.8
              )
              legend(
                "topright",
                legend = paste(
                  c(
                    "Tatal individuals is: ",
                    "Cluster method: ",
                    "Group number: "
                  ),
                  Optimum[c(1:3)],
                  sep = ""
                ),
                lty = 0,
                cex = 1.3,
                bty = "n",
                bg = par("bg")
              )
              dev.off()
            }
          }
          if (!is.null(compress_z))
            write.table(compress_z,
                        paste("GAPIT.Kin.NJtree.compress_z.txt", sep = ""),
                        quote = F)
          print("Kinship NJ TREE PDF created!")
          
          Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "plot NJ TREE")
          Memory = GAPIT.Memory(Memory = Memory, Infor = "plot NJ TREE")
          #rm(hc,clusMember)
        }#end
        ## NJ Tree end    } #end of if(nrow(KI)<1000)
      } #end of if(KI!=1)
    } #end of if(!is.null(KI))
    
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Before SUPER")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "Before SUPER")
    
    #SUPER
    # if(!is.null(GP) & kinship.algorithm=="SUPER" & !is.null(bin.size) & !is.null(inclosure.size))
    # {
    #   mySpecify=GAPIT.Specify(GI=GI,GP=GP,bin.size=bin.size,inclosure.size=inclosure.size)
    #   SNP.QTN=mySpecify$index
    #   if(!is.null(GD))
    #   {
    # 	  GK=GD[,SNP.QTN]
    #     SNPVar=apply(as.matrix(GK),2,var)
    #     GK=GK[,SNPVar>0]
    #     GK=cbind(as.data.frame(GT),as.data.frame(GK)) #add taxa
    #   }
    # }
    # Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Before creating kinship")
    # Memory=GAPIT.Memory(Memory=Memory,Infor="Before creating kinship")
    
    #Create kinship from genotype if not provide
    if (is.null(KI) &
        (!is.null(GD) |
         !is.null(GK)) & !kinship.algorithm %in% c("FarmCPU", "Blink", "MLMM"))
    {
      print("Calculating kinship...")
      if (!is.null(GK))
      {
        thisGD = GK[, -1]
        myGT = as.matrix(GK[, 1])
        print("GK is used to create KI")
      } else{
        thisGD = GD
        myGT = GT
      }
      print(paste(
        "Number of individuals and SNPs are ",
        nrow(thisGD),
        " and ",
        ncol(thisGD)
      ))
      theKin = NULL
      #if(is.null(PCA.col)&!is.null(NJtree.group))PCA.col=rainbow(NJtree.group)[clusMember]
      if (kinship.algorithm == "EMMA")
      {
        half.thisGD = as.matrix(.5 * thisGD)
        if (length(which(is.na(half.thisGD))) > 0)
        {
          print(
            "Substituting missing values with heterozygote for kinship matrrix calculation...."
          )
          half.thisGD[which(is.na(half.thisGD))] = 1
        }
        theKin = emma.kinship(snps = t(as.matrix(.5 * thisGD)),
                              method = "additive",
                              use = "all")
      }
      if (kinship.algorithm == "Loiselle")
        theKin = GAPIT.kinship.loiselle(snps = t(as.matrix(.5 * thisGD)),
                                        method = "additive",
                                        use = "all")
      if (kinship.algorithm == "VanRaden")
        theKin = GAPIT.kinship.VanRaden(snps = as.matrix(thisGD))
      if (kinship.algorithm == "Zhang")
        theKin = GAPIT.kinship.Zhang(snps = as.matrix(thisGD))
      if (kinship.algorithm == "Separation")
      {
        thePCA = GAPIT.PCA(
          X = GD,
          taxa = GT,
          PC.number = PCA.total,
          file.output = F,
          PCA.total = PCA.total,
          PCA.col = NULL,
          PCA.3d = F
        )
        PC = thePCA$PCs[, 1:(1 + PCA.total)]
        theKin = GAPIT.kinship.separation(PCs = thePCA$PCs,
                                          EV = thePCA$EV,
                                          nPCs = PCA.total)
      }
      if (!is.null(theKin))
      {
        colnames(theKin) = myGT
        rownames(theKin) = myGT
        line.names <- myGT
        if (!is.null(NJtree.group))
        {
          distance.matrix = dist(theKin, upper = TRUE)
          hc = hclust(distance.matrix, method = kinship.cluster)
          hcd = as.dendrogram(hc)
          clusMember <- cutree(hc, k = NJtree.group)
          compress_z = table(clusMember, paste(line.names))
          type_col = rainbow(NJtree.group)
          Optimum = c(nrow(theKin), kinship.cluster, NJtree.group)
        }
        print("kinship calculated")
        if (length(GT) < 2000 & file.output)
        {
          #Create heat map for kinship
          print("Creating heat map for kinship...")
          pdf(
            paste("GAPIT.Kin.", kinship.algorithm, ".pdf", sep = ""),
            width = 12,
            height = 12
          )
          par(mar = c(25, 25, 25, 25))
          heatmap.2(
            theKin,
            cexRow = .2,
            cexCol = 0.2,
            col = rev(heat.colors(256)),
            scale = "none",
            symkey = FALSE,
            trace = "none"
          )
          dev.off()
          print("Kinship heat map created")
          ## Jiabo Wang add NJ Tree of kinship at 4.5.2017
          if (!is.null(NJtree.group))
          {
            print("Creating NJ Tree for kinship...")
            for (tr in 1:length(NJtree.type))
            {
              pdf(
                paste("GAPIT.Kin.NJtree.", NJtree.type[tr], ".pdf", sep = ""),
                width = 12,
                height = 12
              )
              par(mar = c(0, 0, 0, 0))
              Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "prepare NJ TREE")
              Memory = GAPIT.Memory(Memory = Memory, Infor = "prepare NJ TREE")
              plot(
                as.phylo(hc),
                type = NJtree.type[tr],
                tip.color = type_col[clusMember],
                use.edge.length = TRUE,
                col = "gray80",
                cex = 0.6
              )
              #legend("topright",legend=c(paste("Tatal numerber of individuals is ",),lty=0,cex=1.3,bty="n",bg=par("bg"))
              legend(
                "topright",
                legend = paste(
                  c(
                    "Tatal individuals is: ",
                    "Group method: ",
                    "Group number: "
                  ),
                  Optimum[c(1:3)],
                  sep = ""
                ),
                lty = 0,
                cex = 1.3,
                bty = "n",
                bg = par("bg")
              )
              dev.off()
            }
            # print(Optimum)
            write.table(compress_z,
                        paste("GAPIT.Kin.NJtree.compress_z.txt", sep = ""),
                        quote = F)
            print("Kinship NJ TREE PDF created!")
            Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "plot NJ TREE")
            Memory = GAPIT.Memory(Memory = Memory, Infor = "plot NJ TREE")
            rm(hc)
          }#end NJtree
        }
        print("Adding IDs to kinship...")
        #Write the kinship into a text file
        KI = cbind(myGT, as.data.frame(theKin)) #This require big memory. Need a way to solve it.
        print("Writing kinship to file...")
        if (file.output)
          write.table(
            KI,
            paste("GAPIT.Kin.", kinship.algorithm, ".csv", sep = ""),
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = FALSE
          )
        print("Kinship save as file")
        rm(theKin)
        gc()
      }
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Estimating kinship")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "Estimating kinship")
      print("Kinship created!")
    }  #end of if(is.null(KI)&!is.null(GD))
    
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "after creating kinship")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "after creating kinship")
    
    
    PC = NULL
    thePCA = NULL
    
    if (PCA.total > 0)
    {
      if (is.null(PCA.col) & !is.null(type_col))
        PCA.col = type_col[clusMember]
      thePCA = GAPIT.PCA(
        X = GD,
        taxa = GT,
        PC.number = PCA.total,
        file.output = file.output,
        PCA.total = PCA.total,
        PCA.col = PCA.col,
        PCA.3d = PCA.3d
      )
      PC = thePCA$PCs[, 1:(1 + PCA.total)]
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "PCA")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "PCA")
      print("PC created")
    }
    
    #LD plot
    #print("LD section")
    if (!is.null(GLD) & file.output)
    {
      if (nrow(GLD) > 500)
      {
        GLD = GLD[1, ]
        print("WARNING: The number of SNPs requested is beyond limitation. No LD plot created.")
      }
      if (nrow(GLD) > 1)
      {
        print("Plot LD...")
        hapmapgeno = data.frame(as.matrix(t(GLD[, -c(1:11)])))
        hapmapgeno[hapmapgeno == "NN"] = NA
        hapmapgeno[hapmapgeno == "XX"] = NA
        hapmapgeno[hapmapgeno == "--"] = NA
        hapmapgeno[hapmapgeno == "++"] = NA
        hapmapgeno[hapmapgeno == "//"] = NA
        LDdist = as.numeric(as.vector(GLD[, 4]))
        LDsnpName = GLD[, 1]
        colnames(hapmapgeno) = LDsnpName
        #Prune SNM names
        #LDsnpName=LDsnpName[GAPIT.Pruning(LDdist,DPP=7)]
        LDsnpName = LDsnpName[c(1, length(LDsnpName))] #keep the first and last snp names only
        #print(hapmapgeno)
        print("Getting genotype object")
        LDsnp = makeGenotypes(hapmapgeno, sep = "", method = as.genotype)   #This need to be converted to genotype object
        print("Caling LDheatmap...")
        pdf(
          paste(
            "GAPIT.LD.chromosom",
            LD.chromosome,
            "(",
            round(max(0, LD.location - LD.range) / 1000000),
            "_",
            round((LD.location + LD.range) / 1000000),
            "Mb)",
            ".pdf",
            sep = ""
          ),
          width = 12,
          height = 12
        )
        #pdf(paste("GAPIT.LD.pdf",sep=""), width = 12, height = 12)
        par(mar = c(25, 25, 25, 25))
        MyHeatmap <-
          try(LDheatmap(
            LDsnp,
            LDdist,
            LDmeasure = "r",
            add.map = TRUE,
            SNP.name = LDsnpName,
            color = rev(cm.colors(20)),
            name = "myLDgrob",
            add.key = TRUE,
            geneMapLabelY = 0.1
          ))
        if (!inherits(MyHeatmap, "try-error"))
        {
          #Modify the plot
          grid.edit(gPath("myLDgrob", "Key", "title"),
                    gp = gpar(cex = .5, col = "blue"))  #edit key title size and color
          grid.edit(gPath("myLDgrob", "geneMap", "title"),
                    gp = gpar(
                      just = c("center", "bottom"),
                      cex = 0.8,
                      col = "black"
                    )) #Edit gene map title
          grid.edit(gPath("myLDgrob", "geneMap", "SNPnames"),
                    gp = gpar(cex = 0.3, col = "black")) #Edit SNP name
        } else{
          print("Warning: error in converting genotype. No LD plot!")
        }
        dev.off()
        print("LD heatmap crated")
      } else{
        # alternative of if(nrow(GLD)>1)
        print("Warning: There are less than two SNPs on the region you sepcified. No LD plot!")
      } #end of #if(nrow(GLD)>1)
    }#end of if(!is.null(GLD))
    
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "after LD plot")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "after LD plot")
    
    
    ###output Marker density and decade of linkage disequilibrium over distance
    if (!is.null(GI) & !is.null(GD) & file.output & Geno.View.output)
    {
      ViewGenotype <- GAPIT.Genotype.View(myGI = GI,
                                          myGD = GD,)
    }
    
    #print("Genotype successfully acomplished")
    return (
      list(
        G = G,
        GD = GD,
        GI = GI,
        GT = GT,
        hasGenotype = hasGenotype,
        genoFormat = genoFormat,
        KI = KI,
        PC = PC,
        byFile = byFile,
        fullGD = fullGD,
        Timmer = Timmer,
        Memory = Memory,
        SNP.QTN = SNP.QTN,
        chor_taxa = chor_taxa
      )
    )
  }

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.0000` <-
  function() {
    ##############################################################################################
    #GAPIT: Genome Association and Prediction Integrated Tool
    #Objective 1: State of art methods for high  power, accuracy and speed;
    #Objective 2: User friendly by design, help documents, and web forum;
    #Objective 3: Comprehensive output to interpret data and results;
    #Objective 4: Informative tables and high quality figures for reports and publication;
    
    #Methods implimented:
    # 1. GLM (Structure or Q method for GWAS, Pritchard et. al. Genetics, 2000)
    # 2. MLM (Q+K, Yu et. al. Nature Genetics, 2006)
    # 3. gBLUP (Marker based kinship, Zhang et. al. Journal of Animal Science, 2007)
    # 4. PCA (Zhao et. al. Plos Genetics, 2007)
    # 5. EMMA (Kang et. al. Genetics, 2008)
    # 6. CMLM (Zhang et. al. Nature Genetics, 2010)
    # 7. EMMAx (Kang et. al. Nature Genetics, 2010)
    # 8. P3D (Zhang et. al. Nature Genetics, 2010)
    # 9. FaST-LMM (Lippert et. al. Nature Methods, 2011)
    # 10. ECMLM (Li et. al. BMC Bioogy, 2014)
    # 11. SUPER (Wang et. al. PLoS One, 2014)
    
    #Designed by Zhiwu Zhang
    #Authors of paper on Bioinformatics (2012, 28:2397-2399): Alex Lipka, Feng Tian, Qishan Wang, Xiaolei Liu, Meng Li,You Tang and Zhiwu Zhang
    #Authors of paper on Plant Genome (2016, Vol 9, No. 2): You Tang, Xiaolei Liu, Jiabo Wang, Meng Li, Qishan Wang, Feng Tian, Zhongbin Su, Yuchun Pan, Di Liu, Alexander E. Lipka, Edward S. Buckler, and Zhiwu Zhang
    library (multtest)
    #if(!require(multtest))
    #{
    #	if (!requireNamespace("BiocManager", quietly = TRUE))
    #    install.packages("BiocManager")
    #    BiocManager::install("multtest")
    #	#source("http://www.bioconductor.org/biocLite.R")
    #    #biocLite("multtest")
    #}
    
    #if(!require(gplots)) install.packages("gplots")
    #if(!require(LDheatmap)) install.packages("LDheatmap")
    #if(!require(genetics)) install.packages("genetics")
    #if(!require(ape)) install.packages("ape")
    #if(!require(compiler)) install.packages("compiler")
    
    #if(!require(EMMREML)) install.packages("EMMREML")
    #if(!require(scatterplot3d)) install.packages("scatterplot3d")
    
    #if(!'multtest'%in% installed.packages()[,"Package"]){
    #	if (!requireNamespace("BiocManager", quietly = TRUE))
    #    install.packages("BiocManager")
    #    BiocManager::install("multtest")
    #    BiocManager::install("snpStats")
    #}
    
    GAPIT.Version = "2020.10.24, GAPIT 3.0"
    print(paste(
      "All packages are loaded already !  ",
      "GAPIT.Version is ",
      GAPIT.Version,
      sep = ""
    ))
    return(GAPIT.Version)
  }

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT` <- function(Y = NULL,
           G = NULL,
           GD = NULL,
           GM = NULL,
           KI = NULL,
           Z = NULL,
           CV = NULL,
           CV.Inheritance = NULL,
           GP = NULL,
           GK = NULL,
           testY = NULL,
           group.from = 1000000 ,
           group.to = 1000000,
           group.by = 20,
           DPP = 100000,
           kinship.cluster = "average",
           kinship.group = 'Mean',
           kinship.algorithm = "VanRaden",
           buspred = FALSE,
           lmpred = FALSE,
           bin.from = 10000,
           bin.to = 10000,
           bin.by = 10000,
           inclosure.from = 10,
           inclosure.to = 10,
           inclosure.by = 10,
           SNP.P3D = TRUE,
           SNP.effect = "Add",
           SNP.impute = "Middle",
           PCA.total = 0,
           SNP.fraction = 1,
           seed = NULL,
           BINS = 20,
           SNP.test = TRUE,
           SNP.MAF = 0,
           FDR.Rate = 1,
           SNP.FDR = 1,
           SNP.permutation = FALSE,
           SNP.CV = NULL,
           SNP.robust = "GLM",
           file.from = 1,
           file.to = 1,
           file.total = NULL,
           file.fragment = 99999,
           file.path = NULL,
           file.G = NULL,
           file.Ext.G = NULL,
           file.GD = NULL,
           file.GM = NULL,
           file.Ext.GD = NULL,
           file.Ext.GM = NULL,
           ngrid = 100,
           llim = -10,
           ulim = 10,
           esp = 1e-10,
           LD.chromosome = NULL,
           LD.location = NULL,
           LD.range = NULL,
           PCA.col = NULL,
           PCA.3d = FALSE,
           NJtree.group = NULL,
           NJtree.type = c("fan", "unrooted"),
           sangwich.top = NULL,
           sangwich.bottom = NULL,
           QC = TRUE,
           GTindex = NULL,
           LD = 0.1,
           plot.bin = 10 ^ 5,
           file.output = TRUE,
           cutOff = 0.05,
           Model.selection = FALSE,
           output.numerical = FALSE,
           output.hapmap = FALSE,
           Create.indicator = FALSE,
           Multi_iter = FALSE,
           num_regwas = 10,
           opt = "extBIC",
           QTN = NULL,
           QTN.round = 1,
           QTN.limit = 0,
           QTN.update = TRUE,
           QTN.method = "Penalty",
           Major.allele.zero = FALSE,
           Random.model = FALSE,
           method.GLM = "FarmCPU.LM",
           method.sub = "reward",
           method.sub.final = "reward",
           method.bin = "static",
           bin.size = c(1000000),
           bin.selection = c(10, 20, 50, 100, 200, 500, 1000),
           memo = NULL,
           Prior = NULL,
           ncpus = 1,
           maxLoop = 3,
           threshold.output = .01,
           Inter.Plot = FALSE,
           Inter.type = c("m", "q"),
           WS = c(1e0, 1e3, 1e4, 1e5, 1e6, 1e7),
           alpha = c(.01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1),
           maxOut = 100,
           QTN.position = NULL,
           CG = NULL,
           converge = 1,
           iteration.output = FALSE,
           acceleration = 0,
           iteration.method = "accum",
           PCA.View.output = TRUE,
           Geno.View.output = TRUE,
           plot.style = "Oceanic",
           SUPER_GD = NULL,
           SUPER_GS = FALSE,
           h2 = NULL,
           NQTN = NULL,
           QTNDist = "normal",
           effectunit = 1,
           category = 1,
           r = 0.25,
           cveff = NULL,
           a2 = 0,
           adim = 2,
           Multiple_analysis = FALSE,
           model = "MLM",
           Para = NULL) 
{
    #Object: To perform GWAS and GPS (Genomic Prediction/Selection)
    #Designed by Zhiwu Zhang
    #Writen by Jiabo Wang
    #Last update: Novenber 3, 2016
    ##############################################################################################
    print("--------------------- Welcome to GAPIT ----------------------------")
    echo = TRUE
    all.memo = NULL
    
    GAPIT.Version = GAPIT.0000()
    #if(!is.null(model))if(!match(model,c("MLM","CMLM","SUPER","GLM","FarmCPU","Blink","BlinkC","MLMM","gBLUP","cBLUP","sBLUP"))) stop(paste("PLease choose one model from ","MLM","CMLM","SUPER","GLM","FarmCPU","Blink","gBLUP","cBLUP","sBLUP",sep=""))
    #Allow either KI or K, but not both
    if (model %in% c("gBLUP", "cBLUP", "sBLUP"))
    {
      SNP.test = FALSE
      SUPER_GS = TRUE
    }
    if (!is.null(KI) &
        is.null(GD) &
        is.null(G) & is.null(file.G) & is.null(file.GD))
      SNP.test = FALSE
    model_store = model
    KI0 = KI
    
    
    print(model_store)
    if (!is.null(Y))
    {
      for (m in 1:length(model_store))
      {
        model = model_store[m]
        if (toupper(model) == "BLINK")
          model = "Blink"
        if (toupper(model) == "FARMCPU")
          model = "FarmCPU"
        if (toupper(model) == "BLINKC")
          model = "BlinkC"
        if (toupper(model) == "GBLUP")
          model = "gBLUP"
        if (toupper(model) == "CBLUP")
          model = "cBLUP"
        if (toupper(model) == "SBLUP")
          model = "sBLUP"
        
        if (group.from < nrow(Y))
          model = "CMLM"
        # }
        if (group.to != group.from)
          model = "CMLM"
        if (group.to == 1 & group.from == 1)
          model = "GLM"
        if (!is.null(sangwich.bottom) &
            !is.null(sangwich.bottom))
          model = "SUPER"
        if (model == "gBLUP")
          model = "MLM"
        if (model == "cBLUP")
          model = "CMLM"
        if (model == "sBLUP")
        {
          model = "SUPER"
          Para$group.from = 1000000
          Para$group.to = 1000000
          Para$group.by = nrow(Y) / 10
        }
        #CMLM
        if (model == "GLM")
        {
          Para$group.from = 1
          Para$group.to = 1
          Para$group.by = group.by
        }
        if (model == "MLM")
        {
          Para$group.from = 1000000
          Para$group.to = 1000000
          Para$group.by = group.by
        }
        if (model == "CMLM")
        {
          if (group.from >= group.to)
            Para$group.from = 1
          Para$group.to = group.to
          Para$group.by = group.by
          #if(Para$group.from==Para$group.to)Para$group.from=10
          print(group.from)
          print(group.to)
        }
        if (model == "SUPER")
        {
          if (!is.null(inclosure.from) &
              is.null(Para$inclosure.from))
            Para$inclosure.from = inclosure.from
          if (is.null(Para$inclosure.from))
            Para$inclosure.from = 10
          if (!is.null(inclosure.to) &
              is.null(Para$inclosure.to))
            Para$inclosure.to = inclosure.to
          if (is.null(Para$inclosure.to))
            Para$inclosure.to = 100
          if (!is.null(inclosure.by) &
              is.null(Para$inclosure.by))
            Para$inclosure.by = inclosure.by
          if (is.null(Para$inclosure.by))
            Para$inclosure.by = 10
          if (!is.null(bin.from) &
              is.null(Para$bin.from))
            Para$bin.from = bin.from
          if (is.null(Para$bin.from))
            Para$bin.from = 10000
          if (!is.null(bin.to) &
              is.null(Para$bin.to))
            Para$bin.to = bin.to
          if (is.null(Para$bin.to))
            Para$bin.to = 10000
          if (!is.null(bin.by) &
              is.null(Para$bin.by))
            Para$bin.by = bin.by
          if (is.null(Para$bin.by))
            Para$bin.by = 10000
          if (!is.null(sangwich.top) &
              is.null(Para$sangwich.top))
            Para$sangwich.top = sangwich.top
          if (is.null(Para$sangwich.top))
            Para$sangwich.top = "MLM"
          if (!is.null(sangwich.bottom) &
              is.null(Para$sangwich.bottom))
            Para$sangwich.bottom = sangwich.bottom
          if (is.null(Para$sangwich.bottom))
            Para$sangwich.bottom = "SUPER"
        }
        if (model == "FarmCPU")
          Para$kinship.algorithm = "FarmCPU"
        if (model == "MLMM")
          Para$kinship.algorithm = "MLMM"
        if (model == "Blink")
          Para$kinship.algorithm = "Blink"
        if (model == "BlinkC")
          Para$kinship.algorithm = "BlinkC"
        if (is.null(memo))
        {
          Para$memo = model
        } else{
          # print(memo)
          # print(model)
          Para$memo = paste(memo, ".", model, sep = "")
        }
        all.memo = c(all.memo, Para$memo)
        # print(Para$memo)
        GAPIT_list = list(
          group.from = group.from ,
          group.to = group.to,
          group.by = group.by,
          DPP = DPP,
          kinship.cluster = kinship.cluster,
          kinship.group = kinship.group,
          kinship.algorithm = kinship.algorithm,
          bin.from = bin.from,
          bin.to = bin.to,
          bin.by = bin.by,
          inclosure.from = inclosure.from,
          inclosure.to = inclosure.to,
          inclosure.by = inclosure.by,
          SNP.P3D = SNP.P3D,
          SNP.effect = SNP.effect,
          SNP.impute = SNP.impute,
          PCA.total = PCA.total,
          SNP.fraction = SNP.fraction,
          seed = seed,
          BINS = 20,
          SNP.test = SNP.test,
          SNP.MAF = SNP.MAF,
          FDR.Rate = FDR.Rate,
          SNP.FDR = SNP.FDR,
          SNP.permutation = SNP.permutation,
          SNP.CV = NULL,
          SNP.robust = "GLM",
          file.from = file.from,
          file.to = file.to,
          file.total = file.total,
          file.fragment = file.fragment,
          file.path = file.path,
          file.G = file.G,
          file.Ext.G = file.Ext.G,
          file.GD = file.GD,
          file.GM = file.GM,
          file.Ext.GD = file.Ext.GD,
          file.Ext.GM = file.Ext.GM,
          ngrid = 100,
          llim = -10,
          ulim = 10,
          esp = 1e-10,
          Inter.Plot = Inter.Plot,
          Inter.type = Inter.type,
          LD.chromosome = LD.chromosome,
          LD.location = LD.location,
          LD.range = LD.range,
          PCA.col = PCA.col,
          PCA.3d = PCA.3d,
          NJtree.group = NJtree.group,
          NJtree.type = NJtree.type,
          opt = opt,
          sangwich.top = sangwich.top,
          sangwich.bottom = sangwich.bottom,
          QC = QC,
          GTindex = GTindex,
          LD = LD,
          plot.bin = plot.bin,
          file.output = file.output,
          cutOff = cutOff,
          Model.selection = Model.selection,
          output.numerical = output.numerical,
          output.hapmap = output.hapmap,
          Create.indicator = Create.indicator,
          QTN = QTN,
          QTN.round = 1,
          QTN.limit = 0,
          QTN.update = TRUE,
          QTN.method = "Penalty",
          Major.allele.zero = Major.allele.zero,
          method.GLM = method.GLM,
          method.sub = method.sub,
          method.sub.final = "reward",
          method.bin = "static",
          bin.size = bin.size,
          bin.selection = bin.selection,
          model = model,
          Random.model = Random.model,
          h2 = h2,
          NQTN = NQTN,
          QTNDist = "normal",
          effectunit = effectunit,
          category = category,
          r = r,
          cveff = NULL,
          a2 = 0,
          adim = 2,
          Multi_iter = Multi_iter,
          num_regwas = num_regwas,
          memo = "",
          Prior = NULL,
          ncpus = 1,
          maxLoop = maxLoop,
          threshold.output = threshold.output,
          WS = c(1e0, 1e3, 1e4, 1e5, 1e6, 1e7),
          alpha = alpha,
          maxOut = 100,
          QTN.position = QTN.position,
          CG = CG,
          converge = converge,
          iteration.output = iteration.output,
          acceleration = 0,
          iteration.method = "accum",
          PCA.View.output = PCA.View.output,
          Geno.View.output = Geno.View.output,
          plot.style = "Oceanic",
          SUPER_GD = NULL,
          SUPER_GS = SUPER_GS,
          Multiple_analysis = Multiple_analysis
        )
        
        G_list_M = rownames(as.matrix(GAPIT_list))
        P_list_M = rownames(as.matrix(Para))
        
        Para = c(GAPIT_list[!G_list_M %in% P_list_M], Para)
        #print(Para$kinship.algorithm)
        if (SUPER_GS == TRUE)
          Para$SNP.test = FALSE
        IC = NULL
        #GAPIT.Version=GAPIT.0000()
        print("--------------------Processing traits----------------------------------")
        # if(!is.null(Y)){
        print("Phenotype provided!")
        if (ncol(Y) < 2)
          stop (
            "Phenotype should have taxa name and one trait at least. Please correct phenotype file!"
          )
        print(paste("The ", m, " model in all.", sep = ""))
        print(model)
        if (m == 1)
        {
          DP = GAPIT.DP(
            G = G,
            GD = GD,
            GM = GM,
            KI = KI0,
            Z = Z,
            CV = CV,
            CV.Inheritance = Para$CV.Inheritance,
            GP = GP,
            GK = GK,
            group.from = Para$group.from ,
            group.to = Para$group.to,
            group.by = Para$group.by,
            DPP = Para$DPP,
            kinship.cluster = Para$kinship.cluster,
            kinship.group = Para$kinship.group,
            kinship.algorithm = Para$kinship.algorithm,
            NJtree.group = Para$NJtree.group,
            NJtree.type = Para$NJtree.type,
            plot.bin = Para$plot.bin,
            PCA.col = Para$PCA.col,
            PCA.3d = Para$PCA.3d,
            sangwich.top = Para$sangwich.top,
            sangwich.bottom = Para$sangwich.bottom,
            LD = Para$LD,
            bin.from = Para$bin.from,
            bin.to = Para$bin.to,
            bin.by = Para$bin.by,
            inclosure.from = Para$inclosure.from,
            inclosure.to = Para$inclosure.to,
            inclosure.by = Para$inclosure.by,
            SNP.P3D = Para$SNP.P3D,
            SNP.effect = Para$SNP.effect,
            SNP.impute = Para$SNP.impute,
            PCA.total = Para$PCA.total,
            SNP.fraction = Para$SNP.fraction,
            seed = Para$seed,
            BINS = Para$BINS,
            SNP.test = Para$SNP.test,
            SNP.MAF = Para$SNP.MAF,
            FDR.Rate = Para$FDR.Rate,
            SNP.FDR = Para$SNP.FDR,
            SNP.permutation = Para$SNP.permutation,
            opt = Para$opt,
            SNP.CV = Para$SNP.CV,
            SNP.robust = Para$SNP.robust,
            Inter.Plot = Para$Inter.Plot,
            Inter.type = Para$Inter.type,
            file.from = Para$file.from,
            file.to = Para$file.to,
            file.total = Para$file.total,
            file.fragment = Para$file.fragment,
            file.path = Para$file.path,
            file.G = Para$file.G,
            file.Ext.G = Para$file.Ext.G,
            file.GD = Para$file.GD,
            file.GM = Para$file.GM,
            file.Ext.GD = Para$file.Ext.GD,
            file.Ext.GM = Para$file.Ext.GM,
            ngrid = Para$ngrid,
            llim = Para$llim,
            ulim = Para$ulim,
            esp = Para$esp,
            Multi_iter = Para$Multi_iter,
            num_regwas = Para$num_regwas,
            LD.chromosome = Para$LD.chromosome,
            LD.location = Para$LD.location,
            LD.range = Para$LD.range,
            QC = Para$QC,
            GTindex = Para$GTindex,
            cutOff = Para$cutOff,
            Model.selection = Para$Model.selection,
            output.numerical = Para$output.numerical,
            Random.model = Para$Random.model,
            Create.indicator = Para$Create.indicator,
            QTN = Para$QTN,
            QTN.round = Para$QTN.round,
            QTN.limit = Para$QTN.limit,
            QTN.update = Para$QTN.update,
            QTN.method = Para$QTN.method,
            Major.allele.zero = Para$Major.allele.zero,
            method.GLM = Para$method.GLM,
            method.sub = Para$method.sub,
            method.sub.final = Para$method.sub.final,
            method.bin = Para$method.bin,
            bin.size = Para$bin.size,
            bin.selection = Para$bin.selection,
            memo = Para$memo,
            Prior = Para$Prior,
            ncpus = Para$ncpus,
            maxLoop = Para$maxLoop,
            threshold.output = Para$threshold.output,
            WS = Para$WS,
            alpha = Para$alpha,
            maxOut = Para$maxOut,
            QTN.position = Para$QTN.position,
            converge = Para$converge,
            iteration.output = Para$iteration.output,
            acceleration = Para$acceleration,
            iteration.method = Para$iteration.method,
            PCA.View.output = Para$PCA.View.output,
            output.hapmap = Para$output.hapmap,
            file.output = Para$file.output,
            Geno.View.output = Para$Geno.View.output,
            plot.style = Para$plot.style,
            SUPER_GD = Para$SUPER_GD,
            SUPER_GS = Para$SUPER_GS,
            CG = Para$CG,
            model = model
          )
        } else{
          DP$kinship.algorithm = Para$kinship.algorithm
          DP$group.from = Para$group.from
          DP$group.to = Para$group.to
          DP$group.by = Para$group.by
          DP$sangwich.top = Para$sangwich.top
          DP$sangwich.bottom = Para$sangwich.bottom
          DP$bin.from = Para$bin.from
          DP$bin.to = Para$bin.to
          DP$bin.by = Para$bin.by
          DP$inclosure.from = Para$inclosure.from
          DP$inclosure.to = Para$inclosure.toDP$inclosure.by = Para$inclosure.by
        }
        
        for (trait in 2:ncol(Y))
        {
          traitname = colnames(Y)[trait]
          ###Statistical distributions of phenotype
          ###Correlation between phenotype and principal components
          print(paste("Processing trait: ", traitname, sep = ""))
          if (!is.null(Para$memo))
            traitname = paste(Para$memo, ".", traitname, sep = "")
          if (!is.null(Y) &
              Para$file.output)
            ViewPhenotype <-
            GAPIT.Phenotype.View(myY = Y[, c(1, trait)],
                                 traitname = traitname,
                                 memo = Para$memo)
          # print(dim(KI0))
          if (!is.null(KI0))
            DP$KI = KI0
          Judge = GAPIT.Judge(
            Y = Y[, c(1, trait)],
            G = DP$G,
            GD = DP$GD,
            KI = DP$KI,
            GM = DP$GM,
            group.to = DP$group.to,
            group.from = DP$group.from,
            sangwich.top = DP$sangwich.top,
            sangwich.bottom = DP$sangwich.bottom,
            kinship.algorithm = DP$kinship.algorithm,
            PCA.total = DP$PCA.total,
            model = DP$model,
            SNP.test = DP$SNP.test
          )
          DP$group.from = Judge$group.from
          DP$group.to = Judge$group.to
          DP$name.of.trait = traitname
          DP$Y = Y[!is.na(Y[, trait]), c(1, trait)]
          DP$model = model
          # print(Para$SNP.test)
          IC = GAPIT.IC(DP = DP)
          SS = GAPIT.SS(
            DP = DP,
            IC = IC,
            buspred = buspred,
            lmpred = lmpred
          )
          if (Para$SNP.test &
              Para$file.output)
            ID = GAPIT.ID(DP = DP,
                          IC = IC,
                          SS = SS)
        }#for loop trait
        #print(SNP.test)
        print("GAPIT accomplished successfully for multiple traits. Result are saved")
        print(
          "It is OK to see this: 'There were 50 or more warnings (use warnings() to see the first 50)'"
        )
        out <- list()
        out$QTN <- QTN.position
        out$GWAS <- SS$GWAS
        out$Pred <- SS$Pred
        out$QTN <- IC$QTN
        out$Power <- SS$Power
        out$FDR <- SS$FDR
        out$Power.Alpha <- SS$Power.Alpha
        out$alpha <- SS$alpha
        out$mc = SS$mc
        out$bc = SS$bc
        out$mp = SS$mp
        out$h2 = SS$h2
        out$PCA = IC$myallCV
        out$GD = DP$GD
        out$GM = DP$GM
        out$KI = IC$K
        out$GM = DP$GM
        out$Compression = SS$Compression
        if (Para$SNP.test == TRUE)
          names(out$GWAS$P.value) = "mp"
        if (kinship.algorithm == "FarmCPU")
          names(out$Pred) = c("Taxa", traitname, "Prediction")
        #return (out)
      }#end of model loop
    } else{
      # is.null(Y)
      #print(Para$SNP.MAF)
      out <- list()
      GAPIT_list = list(
        group.from = group.from ,
        group.to = group.to,
        group.by = group.by,
        DPP = DPP,
        kinship.cluster = kinship.cluster,
        kinship.group = kinship.group,
        kinship.algorithm = kinship.algorithm,
        bin.from = bin.from,
        bin.to = bin.to,
        bin.by = bin.by,
        inclosure.from = inclosure.from,
        inclosure.to = inclosure.to,
        inclosure.by = inclosure.by,
        SNP.P3D = SNP.P3D,
        SNP.effect = SNP.effect,
        SNP.impute = SNP.impute,
        PCA.total = PCA.total,
        SNP.fraction = SNP.fraction,
        seed = seed,
        BINS = 20,
        SNP.test = SNP.test,
        SNP.MAF = SNP.MAF,
        FDR.Rate = FDR.Rate,
        SNP.FDR = SNP.FDR,
        SNP.permutation = SNP.permutation,
        SNP.CV = NULL,
        SNP.robust = "GLM",
        file.from = file.from,
        file.to = file.to,
        file.total = file.total,
        file.fragment = file.fragment,
        file.path = file.path,
        file.G = file.G,
        file.Ext.G = file.Ext.G,
        file.GD = file.GD,
        file.GM = file.GM,
        file.Ext.GD = file.Ext.GD,
        file.Ext.GM = file.Ext.GM,
        ngrid = 100,
        llim = -10,
        ulim = 10,
        esp = 1e-10,
        Inter.Plot = Inter.Plot,
        Inter.type = Inter.type,
        LD.chromosome = LD.chromosome,
        LD.location = LD.location,
        LD.range = LD.range,
        PCA.col = PCA.col,
        PCA.3d = PCA.3d,
        NJtree.group = NJtree.group,
        NJtree.type = NJtree.type,
        opt = opt,
        sangwich.top = sangwich.top,
        sangwich.bottom = sangwich.bottom,
        QC = QC,
        GTindex = GTindex,
        LD = LD,
        plot.bin = plot.bin,
        file.output = file.output,
        cutOff = cutOff,
        Model.selection = Model.selection,
        output.numerical = output.numerical,
        output.hapmap = output.hapmap,
        Create.indicator = Create.indicator,
        QTN = QTN,
        QTN.round = 1,
        QTN.limit = 0,
        QTN.update = TRUE,
        QTN.method = "Penalty",
        Major.allele.zero = Major.allele.zero,
        method.GLM = method.GLM,
        method.sub = method.sub,
        method.sub.final = "reward",
        method.bin = "static",
        bin.size = bin.size,
        bin.selection = bin.selection,
        model = model,
        Random.model = Random.model,
        h2 = h2,
        NQTN = NQTN,
        QTNDist = "normal",
        effectunit = effectunit,
        category = category,
        r = r,
        cveff = NULL,
        a2 = 0,
        adim = 2,
        Multi_iter = Multi_iter,
        num_regwas = num_regwas,
        memo = "",
        Prior = NULL,
        ncpus = 1,
        maxLoop = maxLoop,
        threshold.output = threshold.output,
        WS = c(1e0, 1e3, 1e4, 1e5, 1e6, 1e7),
        alpha = alpha,
        maxOut = 100,
        QTN.position = QTN.position,
        CG = CG,
        converge = converge,
        iteration.output = iteration.output,
        acceleration = 0,
        iteration.method = "accum",
        PCA.View.output = PCA.View.output,
        Geno.View.output = Geno.View.output,
        plot.style = "Oceanic",
        SUPER_GD = NULL,
        SUPER_GS = SUPER_GS,
        Multiple_analysis = Multiple_analysis
      )
      if (model == "MLM")
      {
        Para$group.from = 1000000
        Para$group.to = 1000000
        Para$group.by = group.by
      }
      G_list_M = rownames(as.matrix(GAPIT_list))
      P_list_M = rownames(as.matrix(Para))
      if (is.null(memo))
      {
        Para$memo = model
      } else{
        Para$memo = paste(memo, ".", mode, sep = "")
      }
      all.memo = c(all.memo, Para$memo)
      Para = c(GAPIT_list[!G_list_M %in% P_list_M], Para)
      myGenotype <-
        GAPIT.Genotype(
          G = G,
          GD = GD,
          GM = GM,
          KI = KI,
          kinship.algorithm = kinship.algorithm,
          PCA.total = PCA.total,
          SNP.fraction = SNP.fraction,
          SNP.test = SNP.test,
          file.path = file.path,
          file.from = file.from,
          file.to = file.to,
          file.total = file.total,
          file.fragment = file.fragment,
          file.G = file.G,
          file.Ext.G = file.Ext.G,
          file.GD = file.GD,
          file.GM = file.GM,
          file.Ext.GD = file.Ext.GD,
          file.Ext.GM = file.Ext.GM,
          SNP.MAF = SNP.MAF,
          FDR.Rate = FDR.Rate,
          SNP.FDR = SNP.FDR,
          SNP.effect = SNP.effect,
          SNP.impute = SNP.impute,
          NJtree.group = NJtree.group,
          NJtree.type = NJtree.type,
          LD.chromosome = LD.chromosome,
          LD.location = LD.location,
          LD.range = LD.range,
          GP = GP,
          GK = GK,
          bin.size = NULL,
          inclosure.size = NULL,
          sangwich.top = NULL,
          sangwich.bottom = sangwich.bottom,
          GTindex = NULL,
          file.output = file.output,
          Create.indicator = Create.indicator,
          Major.allele.zero = Major.allele.zero,
          Geno.View.output = Geno.View.output,
          PCA.col = PCA.col,
          PCA.3d = PCA.3d
        )
      GD = myGenotype$GD
      GI = myGenotype$GI
      GT = myGenotype$GT
      #G=myGenotype$G
      chor_taxa = myGenotype$chor_taxa
      rownames(GD) = GT
      colnames(GD) = GI[, 1]
      taxa = GT
      if (!is.null(chor_taxa))
      {
        chro = as.numeric(as.matrix(GI[, 2]))
        for (i in 1:length(chro))
        {
          chro[chro == i] = chor_taxa[i]
        }
        GI[, 2] = chro
      }
      #print(GD[1:5,1:5])
      if (output.numerical)
      {
        write.table(
          cbind(taxa, GD),
          "GAPIT.Genotype.Numerical.txt",
          quote = FALSE,
          sep = "\t",
          row.names = F,
          col.names = T
        )
        write.table(
          GI,
          "GAPIT.Genotype.map.txt",
          quote = FALSE,
          sep = "\t",
          row.names = F,
          col.names = T
        )
      }
      if (output.hapmap)
        write.table(
          myGenotype$G,
          "GAPIT.Genotype.hmp.txt",
          quote = FALSE,
          sep = "\t",
          row.names = FALSE,
          col.names = FALSE
        )
      #GD=cbind(as.data.frame(GT),GD)
      if (!is.null(seed))
        set.seed(seed)
      #print(Para$NQTN)
      if (!is.null(Para$NQTN) & !is.null(Para$h2))
      {
        myG_simulation <-
          GAPIT.Phenotype.Simulation(
            GD = cbind(as.data.frame(myGenotype$GT), myGenotype$GD),
            GM = myGenotype$GI,
            h2 = Para$h2,
            NQTN = Para$NQTN,
            QTNDist = Para$QTNDist,
            effectunit = Para$effectunit,
            category = Para$category,
            r = Para$r,
            cveff = Para$cveff,
            a2 = Para$a2,
            adim = Para$adim
          )
        out = c(out, myG_simulation)
      }
      out$GD = data.frame(cbind(as.data.frame(GT), as.data.frame(GD)))
      out$GM = GI
      out$G = myGenotype$G
      out$kinship = myGenotype$KI
      out$PCA = myGenotype$PC
      out$chor_taxa = chor_taxa
    }# is.null(Y)
    
    #print(tail(IC$GM))
    model_store = all.memo
    if (!is.null(Y) &
        SNP.test)
      if(Multiple_analysis &
         Para$file.output &
         length(model_store) * (ncol(Y) - 1) > 1 &
         length(model_store) * (ncol(Y) - 1) < 9)
      {
        #print(DP$QTN.position)
        GMM = GAPIT.Multiple.Manhattan(
          model_store = model_store,
          Y = Y,
          GM = IC$GM,
          seqQTN = DP$QTN.position,
          cutOff = DP$cutOff
        )
        #print(str(GMM$multip_mapP))
        GAPIT.Circle.Manhatton.Plot(
          band = 1,
          r = 3,
          GMM$multip_mapP,
          plot.type = c("c", "q"),
          signal.line = 1,
          xz = GMM$xz,
          threshold = DP$cutOff
        )
      }# end of mutiple manhantton plot
    
    if (file.output &
        !SNP.test & model_store %in% c("gBLUP", "cBLUP", "sBLUP") &
        Inter.Plot)
    {
      print("here will start interactive for GS !!!")
      GAPIT.Interactive.GS(model_store = model_store, Y = Y)
      if (!is.null(testY))
        GAPIT.Interactive.GS(model_store = model_store,
                             Y = Y,
                             testY = testY)
      
      #print(str(GMM$multip_mapP))
      #GAPIT.Circle.Manhatton.Plot(band=1,r=3,GMM$multip_mapP,plot.type=c("c","q"),signal.line=1,xz=GMM$xz,threshold=DP$cutOff)
    }# end of mutiple manhantton plot
    
    return (out)
  }  #end of GAPIT function

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.DP` <-
  function(G = NULL,
           GD = NULL,
           GM = NULL,
           KI = NULL,
           Z = NULL,
           CV = NULL,
           CV.Inheritance = NULL,
           GP = NULL,
           GK = NULL,
           group.from = 30 ,
           group.to = 1000000,
           group.by = 10,
           DPP = 100000,
           kinship.cluster = "average",
           kinship.group = 'Mean',
           kinship.algorithm = "VanRaden",
           bin.from = 10000,
           bin.to = 10000,
           bin.by = 10000,
           inclosure.from = 10,
           inclosure.to = 10,
           inclosure.by = 10,
           SNP.P3D = TRUE,
           SNP.effect = "Add",
           SNP.impute = "Middle",
           PCA.total = 0,
           SNP.fraction = 1,
           seed = 123,
           BINS = 20,
           SNP.test = TRUE,
           SNP.MAF = 0,
           FDR.Rate = 1,
           SNP.FDR = 1,
           SNP.permutation = FALSE,
           SNP.CV = NULL,
           SNP.robust = "GLM",
           NJtree.group = NULL,
           NJtree.type = c("fan", "unrooted"),
           plot.bin = 10 ^ 6,
           PCA.col = NULL,
           PCA.3d = FALSE,
           file.from = 1,
           file.to = 1,
           file.total = NULL,
           file.fragment = 99999,
           file.path = NULL,
           Inter.Plot = FALSE,
           Inter.type = c("m", "q"),
           file.G = NULL,
           file.Ext.G = NULL,
           file.GD = NULL,
           file.GM = NULL,
           file.Ext.GD = NULL,
           file.Ext.GM = NULL,
           ngrid = 100,
           llim = -10,
           ulim = 10,
           esp = 1e-10,
           Multi_iter = FALSE,
           num_regwas = 10,
           LD.chromosome = NULL,
           LD.location = NULL,
           LD.range = NULL,
           p.threshold = NA,
           QTN.threshold = 0.01,
           maf.threshold = 0.03,
           sangwich.top = NULL,
           sangwich.bottom = NULL,
           QC = TRUE,
           GTindex = NULL,
           LD = 0.1,
           opt = "extBIC",
           file.output = FALSE,
           cutOff = 0.01,
           Model.selection = FALSE,
           output.numerical = FALSE,
           Random.model = FALSE,
           output.hapmap = FALSE,
           Create.indicator = FALSE,
           QTN = NULL,
           QTN.round = 1,
           QTN.limit = 0,
           QTN.update = TRUE,
           QTN.method = "Penalty",
           Major.allele.zero = FALSE,
           method.GLM = "fast.lm",
           method.sub = "reward",
           method.sub.final = "reward",
           method.bin = "static",
           bin.size = c(1000000),
           bin.selection = c(10, 20, 50, 100, 200, 500, 1000),
           memo = "",
           Prior = NULL,
           ncpus = 1,
           maxLoop = 3,
           threshold.output = .01,
           WS = c(1e0, 1e3, 1e4, 1e5, 1e6, 1e7),
           alpha = c(.01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1),
           maxOut = 100,
           QTN.position = NULL,
           converge = 1,
           iteration.output = FALSE,
           acceleration = 0,
           iteration.method = "accum",
           PCA.View.output = TRUE,
           Geno.View.output = TRUE,
           plot.style = "Oceanic",
           SUPER_GD = NULL,
           SUPER_GS = FALSE,
           CG = NULL,
           model = "MLM") {
    #Object: To Data and Parameters
    #Designed by Zhiwu Zhang
    #Writen by Jiabo Wang
    #Last update: Novenber 3, 2016
    ##############################################################################################
    print("GAPIT.DP in process...")
    #Judge phenotype  genotype and GAPIT logical
    #print(file.from)
    #print(kinship.algorithm)
    #print(NJtree.group)
    myGenotype <-
      GAPIT.Genotype(
        G = G,
        GD = GD,
        GM = GM,
        KI = KI,
        PCA.total = PCA.total,
        kinship.algorithm = kinship.algorithm,
        SNP.fraction = SNP.fraction,
        SNP.test = FALSE,
        file.path = file.path,
        file.from = file.from,
        file.to = file.to,
        file.total = file.total,
        file.fragment = file.fragment,
        file.G = file.G,
        file.Ext.G = file.Ext.G,
        file.GD = file.GD,
        file.GM = file.GM,
        file.Ext.GD = file.Ext.GD,
        file.Ext.GM = file.Ext.GM,
        SNP.MAF = SNP.MAF,
        FDR.Rate = FDR.Rate,
        SNP.FDR = SNP.FDR,
        SNP.effect = SNP.effect,
        SNP.impute = SNP.impute,
        NJtree.group = NJtree.group,
        NJtree.type = NJtree.type,
        LD.chromosome = LD.chromosome,
        LD.location = LD.location,
        LD.range = LD.range,
        GP = GP,
        GK = GK,
        bin.size = NULL,
        inclosure.size = NULL,
        sangwich.top = sangwich.top,
        sangwich.bottom = sangwich.bottom,
        GTindex = NULL,
        file.output = file.output,
        Create.indicator = Create.indicator,
        Major.allele.zero = Major.allele.zero,
        Geno.View.output = Geno.View.output,
        PCA.col = PCA.col,
        PCA.3d = PCA.3d
      )
    
    # }
    
    KI = myGenotype$KI
    PC = myGenotype$PC
    print(dim(PC))
    
    genoFormat = myGenotype$genoFormat
    hasGenotype = myGenotype$hasGenotype
    byFile = myGenotype$byFile
    fullGD = myGenotype$fullGD
    GD = myGenotype$GD
    GI = myGenotype$GI
    
    GT = myGenotype$GT
    G = myGenotype$G
    chor_taxa = myGenotype$chor_taxa
    
    #if G exist turn to GD and GM
    
    if (output.numerical)
      write.table(
        GD,
        "GAPIT.Genotype.Numerical.txt",
        quote = FALSE,
        sep = "\t",
        row.names = TRUE,
        col.names = NA
      )
    if (output.hapmap)
      write.table(
        myGenotype$G,
        "GAPIT.Genotype.hmp.txt",
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE
      )
    
    
    rownames(GD) = GT
    colnames(GD) = GI[, 1]
    GD = cbind(as.data.frame(GT), GD)
    
    print("GAPIT.DP accomplished successfully for multiple traits. Results are saved")
    return (
      list(
        Y = NULL,
        G = G,
        GD = GD,
        GM = GI,
        KI = KI,
        Z = Z,
        CV = CV,
        CV.Inheritance = CV.Inheritance,
        GP = GP,
        GK = GK,
        PC = PC,
        GI = GI,
        group.from = group.from ,
        group.to = group.to,
        group.by = group.by,
        DPP = DPP,
        name.of.trait = NULL,
        kinship.cluster = kinship.cluster,
        kinship.group = kinship.group,
        kinship.algorithm = kinship.algorithm,
        NJtree.group = NJtree.group,
        NJtree.type = NJtree.type,
        PCA.col = PCA.col,
        PCA.3d = PCA.3d,
        bin.from = bin.from,
        bin.to = bin.to,
        bin.by = bin.by,
        inclosure.from = inclosure.from,
        inclosure.to = inclosure.to,
        inclosure.by = inclosure.by,
        opt = opt,
        SNP.P3D = SNP.P3D,
        SNP.effect = SNP.effect,
        SNP.impute = SNP.impute,
        PCA.total = PCA.total,
        SNP.fraction = SNP.fraction,
        seed = seed,
        BINS = BINS,
        SNP.test = SNP.test,
        SNP.MAF = SNP.MAF,
        FDR.Rate = FDR.Rate,
        SNP.FDR = SNP.FDR,
        SNP.permutation = SNP.permutation,
        SNP.CV = SNP.CV,
        SNP.robust = SNP.robust,
        file.from = file.from,
        file.to = file.to,
        file.total = file.total,
        file.fragment = file.fragment,
        file.path = file.path,
        file.G = file.G,
        file.Ext.G = file.Ext.G,
        file.GD = file.GD,
        file.GM = file.GM,
        file.Ext.GD = file.Ext.GD,
        file.Ext.GM = file.Ext.GM,
        ngrid = ngrid,
        llim = llim,
        ulim = ulim,
        esp = esp,
        Inter.Plot = Inter.Plot,
        Inter.type = Inter.type,
        LD.chromosome = LD.chromosome,
        LD.location = LD.location,
        LD.range = LD.range,
        Multi_iter = Multi_iter,
        sangwich.top = sangwich.top,
        sangwich.bottom = sangwich.bottom,
        QC = QC,
        GTindex = GTindex,
        LD = LD,
        GT = GT,
        file.output = file.output,
        cutOff = cutOff,
        Model.selection = Model.selection,
        output.numerical = output.numerical,
        output.hapmap = output.hapmap,
        Create.indicator = Create.indicator,
        Random.model = Random.model,
        QTN = QTN,
        QTN.round = QTN.round,
        QTN.limit = QTN.limit,
        QTN.update = QTN.update,
        QTN.method = QTN.method,
        Major.allele.zero = Major.allele.zero,
        method.GLM = method.GLM,
        method.sub = method.sub,
        method.sub.final = method.sub.final,
        method.bin = method.bin,
        bin.size = bin.size,
        bin.selection = bin.selection,
        memo = memo,
        Prior = Prior,
        ncpus = 1,
        maxLoop = maxLoop,
        threshold.output = threshold.output,
        WS = WS,
        alpha = alpha,
        maxOut = maxOut,
        QTN.position = QTN.position,
        converge = 1,
        iteration.output = iteration.output,
        acceleration = 0,
        iteration.method = iteration.method,
        PCA.View.output = PCA.View.output,
        p.threshold = p.threshold,
        QTN.threshold = QTN.threshold,
        maf.threshold = maf.threshold,
        chor_taxa = chor_taxa,
        num_regwas = num_regwas,
        Geno.View.output = Geno.View.output,
        plot.style = plot.style,
        SUPER_GD = SUPER_GD,
        SUPER_GS = SUPER_GS,
        CG = CG,
        plot.bin = plot.bin
      )
    )
  }  #end of GAPIT DP function


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.kinship.VanRaden` <-
  function(snps, hasInbred = TRUE) {
    # Object: To calculate the kinship matrix using the method of VanRaden (2009, J. Dairy Sci. 91:4414???C4423)
    # Input: snps is n individual rows by m snps columns
    # Output: n by n relationship matrix
    # Authors: Zhwiu Zhang
    # Last update: March 2, 2016
    ##############################################################################################
    print("Calculating kinship with VanRaden method...")
    #Remove invariants
    fa = colSums(snps) / (2 * nrow(snps))
    index.non = fa >= 1 | fa <= 0
    snps = snps[, !index.non]
    
    nSNP = ncol(snps)
    nInd = nrow(snps)
    n = nInd
    
    ##allele frequency of second allele
    p = colSums(snps) / (2 * nInd)
    P = 2 * (p - .5) #Difference from .5, multiple by 2
    snps = snps - 1 #Change from 0/1/2 coding to -1/0/1 coding
    
    print("substracting P...")
    Z = t(snps) - P#operation on matrix and vector goes in direction of column
    print("Getting X'X...")
    #K=tcrossprod((snps), (snps))
    K = crossprod((Z), (Z)) #Thanks to Peng Zheng, Meng Huang and Jiafa Chen for finding the problem
    
    print("Adjusting...")
    adj = 2 * sum(p * (1 - p))
    K = K / adj
    
    print("Calculating kinship with VanRaden method: done")
    
    return(K)
  }

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.Judge` <-
  function(Y = Y,
           G = NULL,
           GD = NULL,
           KI = NULL,
           GM = NULL,
           group.to = group.to,
           group.from = group.from,
           sangwich.top = sangwich.top,
           sangwich.bottom = sangwich.bottom,
           kinship.algorithm = kinship.algorithm,
           PCA.total = PCA.total,
           model = "MLM",
           SNP.test = TRUE) {
    #Object: To judge Pheno and Geno data practicability
    #Designed by Zhiwu Zhang
    #Writen by Jiabo Wang
    #Last update: Novenber 3, 2016
    ##############################################################################################
    print("--------------------Phenotype and Genotype ----------------------------------")
    if (ncol(Y) < 2)
      stop ("Phenotype should have taxa name and one trait at least. Please correct phenotype file!")
    print(kinship.algorithm)
    if (is.null(KI) &
        is.null(GD) &
        kinship.algorithm != "SUPER" &
        is.null(G))
      stop (
        "GAPIT says: Kinship is required. As genotype is not provided, kinship can not be created."
      )
    if (kinship.algorithm == "FarmCPU" &
        SNP.test == FALSE)
      stop("FarmCPU is only for GWAS, plase set: SNP.test= TRUE")
    #if((!is.null(GD))&(!is.null(G))) stop("GAPIT Says:Please put in only one type of geno data.")
    if (is.null(GD) &
        is.null(G) & is.null(KI))
      stop ("GAPIT Says:GAPIT need genotype!!!")
    if (!is.null(GD) &
        is.null(GM) &
        (is.null(G)) &
        SNP.test)
      stop("GAPIT Says: Genotype data and map files should be in pair")
    if (is.null(GD) &
        !is.null(GM) &
        (is.null(G)) &
        SNP.test)
      stop("GAPIT Says: Genotype data and map files should be in pair")
    
    if (!is.null(GD) & !is.null(Y))
    {
      if (is.null(GD[, 1] %in% Y[, 1]))
        stop("GAPIT Says: There are no common taxa between genotype and phenotype")
    }
    if (!is.null(G) & !is.null(Y))
    {
      if (is.null(colnames(G)[-c(1:11)] %in% Y[, 1]))
        stop("GAPIT Says: There are no common taxa between genotype and phenotype")
    }
    
    if (!is.null(Y))
      nY = nrow(Y)
    if (!is.null(Y))
      ntrait = ncol(Y) - 1
    print(paste("There are ", ntrait, " traits in phenotype data."))
    print(paste("There are ", nY, " individuals in phenotype data."))
    if (!is.null(G))
      nG = nrow(G) - 11
    if (!is.null(GD))
    {
      nG = ncol(GD) - 1
      print(paste("There are ", nG, " markers in genotype data."))
    }
    print("Phenotype and Genotype are test OK !!")
    
    print("--------------------GAPIT Logical ----------------------------------")
    #if (group.to>nY&is.null(KI))group.to=nY
    #if (group.from>group.to&is.null(KI)) group.from=group.to
    if (!is.null(sangwich.top) &
        is.null(sangwich.bottom))
      stop("GAPIT Says: SUPER method need sangwich.top and bottom")
    if (is.null(sangwich.top) &
        !is.null(sangwich.bottom))
      stop("GAPIT Says: SUPER method need sangwich.top and bottom")
    if (kinship.algorithm == "Separation" &
        PCA.total == 0)
      stop ("GAPIT Says: Separation kinship need PCA.total>0")
    
    
    
    
    
    
    
    return (list(group.to = group.to, group.from = group.from))
  }#end of GAPIT.Pheno.Geno.judge function

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.IC` <-
  function(DP = NULL) {
    #Object: To Intermediate Components
    #Designed by Zhiwu Zhang
    #Writen by Jiabo Wang
    #Last update: Novenber 3, 2016
    ##############################################################################################
    print("GAPIT.IC in process...")
    
    Y = DP$Y
    PC = DP$PC
    CV = DP$CV
    GD = DP$GD
    # print(dim(GD))
    noCV = FALSE
    if (is.null(CV)) {
      noCV = TRUE
      CV = GD[, 1:2]
      CV[, 2] = 1
      colnames(CV) = c("taxa", "overall")
      print(paste("There is 0 Covarinces.", sep = ""))
      
    }
    
    Y = Y[!is.na(Y[, 2]), ]
    taxa_Y = as.character(Y[, 1])
    # print(head(Y))
    if (DP$PCA.total > 0 &
        !is.null(DP$CV))
      CV = GAPIT.CVMergePC(DP$CV, PC)
    if (DP$PCA.total > 0 & is.null(DP$CV))
      CV = PC
    # print(all.equal(as.character(GD[,1]),as.character(CV[,1])))
    # print(all.equal(as.character(DP$CV[,1]),as.character(CV[,1])))
    # print(all.equal(as.character(PC[,1]),as.character(CV[,1])))
    # print(dim(GD))
    # print(dim(DP$KI))
    if (ncol(GD) == 0 & !is.null(DP$KI))
    {
      taxa_KI = as.character(DP$KI[, 1])
      
      taxa_CV = as.character(CV[, 1])
      taxa_comall = intersect(intersect(taxa_KI, taxa_Y), taxa_CV)
      # print(length(taxa_comall))
      comCV = CV[taxa_CV %in% taxa_comall, ]
      comCV <- comCV[match(taxa_comall, as.character(comCV[, 1])), ]
      comY = Y[taxa_Y %in% taxa_comall, ]
      comY <- comY[match(taxa_comall, as.character(comY[, 1])), ]
      # print(head(comY))
      # print(dim(comY))
      comGD = NULL
      
    } else{
      taxa_GD = as.character(GD[, 1])
      taxa_comGD = as.character(GD[, 1])
      taxa_CV = as.character(CV[, 1])
      # print(length(taxa_GD))
      # print(length(taxa_CV))
      # print(length(taxa_Y))
      taxa_comall = intersect(intersect(taxa_GD, taxa_Y), taxa_CV)
      comCV = CV[taxa_CV %in% taxa_comall, ]
      comCV <- comCV[match(taxa_comall, as.character(comCV[, 1])), ]
      
      comGD = GD[taxa_GD %in% taxa_comall, ]
      comGD <- comGD[match(taxa_comall, as.character(comGD[, 1])), ]
      
      comY = Y[taxa_Y %in% taxa_comall, ]
      comY <- comY[match(taxa_comall, as.character(comY[, 1])), ]
      # print("@@@@@")
      # print(dim(comY))
      # print(all.equal(as.character(Y[,1]),as.character(comCV[,1])))
      # print(all.equal(as.character(GD[,1]),as.character(comCV[,1])))
      # print(all.equal(as.character(comY[,1]),as.character(comCV[,1])))
    }
    
    
    GT = as.matrix(as.character(taxa_comall))
    print(
      paste(
        "There are ",
        length(GT),
        " common individuals in genotype , phenotype and CV files.",
        sep = ""
      )
    )
    
    if (nrow(comCV) != length(GT))
      stop (
        "GAPIT says: The number of individuals in CV does not match to the number of individuals in genotype files."
      )
    
    print("The dimension of total CV is ")
    print(dim(comCV))
    
    print("GAPIT.IC accomplished successfully for multiple traits. Results are saved")
    if (DP$kinship.algorithm %in% c("FarmCPU", "Blink", "MLMM")) {
      return (
        list(
          Y = comY,
          GT = GT,
          PCA = comCV,
          K = DP$KI,
          GD = comGD,
          GM = DP$GM,
          myallCV = CV,
          myallGD = GD
        )
      )
    } else{
      return (
        list(
          Y = comY,
          GT = GT,
          PCA = comCV,
          K = DP$KI,
          GD = comGD,
          GM = DP$GM,
          myallCV = CV,
          myallGD = GD,
          myallY = Y
        )
      )
    }
  }  #end of GAPIT IC function

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.SS` <-
  function(DP = NULL,
           IC = NULL,
           buspred = FALSE,
           lmpred = TRUE) {
    #Object: To Sufficient Statistics (SS) for GWAS and GS
    #Designed by Zhiwu Zhang
    #Writen by Jiabo Wang
    #Last update: Novenber 3, 2016
    ##############################################################################################
    print("GAPIT.SS in process...")
    #Define the funcitno here
    Timmer = GAPIT.Timmer(Infor = "GAPIT.SS")
    Memory = GAPIT.Memory(Infor = "GAPIT.SS")
    GR = list()
    GR$GVs = NULL
    
    if (DP$SNP.test)
    {
      ic_GD = IC$GD
      ic_GM = IC$GM
      ic_Y = IC$Y
      ic_KI = IC$K
      ic_PCA = IC$PCA
      Z = DP$Z
      
      taxa_Y = as.character(ic_Y[, 1])
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "GAPIT.QC")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "GAPIT.QC")
      
      if (DP$kinship.algorithm != "None" &
          DP$kinship.algorithm != "SUPER" & is.null(Z))
      {
        Z = as.data.frame(diag(1, nrow(ic_Y)))
        Z = rbind(taxa_Y, Z)
        taxa = c('Taxa', as.character(taxa_Y))
        Z = cbind(taxa, Z)
      }
      # print(head(ic_PCA))
      # print(dim(DP$CV))
      # print(head(DP$PC))
      if (max(ic_PCA[, 2]) == min(ic_PCA[, 2]))
        ic_PCA = NULL
      #print(head(ic_PCA))
      # print("@@@@@")
      # print(DP$kinship.algorithm)
      if (DP$SNP.test &
          DP$kinship.algorithm %in% c("FarmCPU", "Blink", "MLMM", "BlinkC"))
      {
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "GAPIT.FarmCPU")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "GAPIT.FarmCPU")
        myBus = GAPIT.Bus(
          Y = ic_Y,
          CV = ic_PCA,
          Z = NULL,
          GK = NULL,
          KI = ic_KI,
          GD = ic_GD,
          GM = ic_GM,
          GT = IC$GT,
          method = DP$kinship.algorithm,
          GTindex = DP$GTindex,
          LD = DP$LD,
          opt = DP$opt,
          bin.size = DP$bin.size,
          bin.selection = DP$bin.selection,
          alpha = DP$alpha,
          WS = DP$WS,
          cutOff = DP$cutOff,
          p.threshold = DP$p.threshold,
          QTN.threshold = DP$QTN.threshold,
          maf.threshold = DP$maf.threshold,
          method.GLM = DP$method.GLM,
          method.sub = DP$method.sub,
          method.sub.final = DP$method.sub.final,
          method.bin = DP$method.bin,
          Random.model = DP$Random.model,
          DPP = DP$DPP,
          file.output = DP$file.output,
          Multi_iter = DP$Multi_iter,
          num_regwas = DP$num_regwas
        )
        GWAS = myBus$GWAS
        Pred = myBus$Pred
        
        # BUS Prediction with gBLUP
        # lmpred=TRUE
        if (!is.null(Pred))
          buspred = FALSE
        print(myBus$seqQTN)
        if (buspred)
        {
          X = DP$GD[, -1]
          # print(dim(X))
          # print(dim(IC$myallCV))
          # print(dim(ic_PCA))
          if (lmpred)
          {
            print("Linear Regression to Predict phenotype !!")
            # colnames(busCV)[1]=c("Taxa")
            # print(length(IC$GT))
            index = as.character(DP$GD[, 1]) %in% as.character(ic_Y[, 1])
            # print(cbind(ic_Y,IC$PCA))
            if (!is.null(myBus$seqQTN))
            {
              # busCV=cbind(IC$myallCV,X[,myBus$seqQTN])
              GD1 = as.matrix(X[index, myBus$seqQTN])
              GD2 = as.matrix(X[, myBus$seqQTN])
            } else{
              numMarker = nrow(GWAS)
              sp = sort(GWAS$P.value)
              spd = abs(DP$cutOff - sp)
              index_fdr = grep(min(spd), spd)[1]
              FDRcutoff = DP$cutOff * index_fdr / numMarker
              seqQTN = as.numeric(rownames(GWAS[GWAS$P.value < FDRcutoff, ]))
              # busCV=cbind(IC$myallCV,X[,seqQTN])
              GD1 = as.matrix(X[index, seqQTN])
              GD2 = as.matrix(X[, seqQTN])
            }
            if (!is.null(IC$myallCV))
            {
              CV1 = as.matrix(IC$PCA[, -1])
              Group = 1:nrow(DP$GD)
              RefInf = rep(2, nrow(DP$GD))
              print(table(index))
              RefInf[index] = 1
              ID = 1:nrow(IC$myallCV)
              BLUP = rep(NA, nrow(DP$GD))
              PEV = rep(NA, nrow(DP$GD))
              BLUE = rep(NA, nrow(DP$GD))
              print("The dimension of CV in lm model :")
              print(dim(CV1))
              print(dim(GD1))
              # print(ic_Y[!is.na(ic_Y[,2]),2])
              mylm = lm(ic_Y[, 2] ~ cbind(CV1, GD1))
              print(cor(ic_Y[, 2], as.numeric(
                predict(mylm, as.data.frame(cbind(
                  CV1, GD1
                )))
              )))
              # Pred = cbind(as.data.frame(DP$GD[index,1]),as.data.frame(predict(mylm,as.data.frame(cbind(CV1,GD1)))))
              # colnames(Pred)=c("Taxa","Prediction")
              # print(mylm$coefficients)
              # print(head(cbind(IC$myallCV,GD2))
              if (var(IC$myallCV[, 2]) == 0)
              {
                kk = 1:2
              } else{
                kk = 1
              }
              aa = as.numeric(mylm$coefficients[-kk] %*% t(as.matrix(cbind(
                IC$myallCV[, -kk], GD2
              ))))
              # print(aa)
              pred0 = cbind(Group,
                            RefInf,
                            ID,
                            BLUP,
                            PEV,
                            BLUE,
                            as.data.frame(aa))
              Pred = cbind(as.data.frame(DP$GD[, 1]), as.matrix(pred0))
              colnames(Pred) = c("Taxa",
                                 "Group",
                                 "RefInf",
                                 "ID",
                                 "BLUP",
                                 "PEV",
                                 "BLUE",
                                 "Prediction")
              
            } else{
              busCV = cbind(as.data.frame(DP$GD[, 1]), X[, myBus$seqQTN])
              CV1 = NULL
              Group = 1:nrow(IC$myallCV)
              RefInf = rep(2, nrow(IC$myallCV))
              RefInf[index] = 1
              ID = 1:nrow(IC$myallCV)
              BLUP = NA
              PEV = NA
              BLUE = NA
              print("The dimension of CV in lm model :")
              print(dim(CV1))
              print(dim(GD1))
              # print(dim(GD1))
              # print(ic_Y[!is.na(ic_Y[,2]),2])
              mylm = lm(ic_Y[!is.na(ic_Y[, 2]), 2] ~ GD1)
              # print("!!")
              print(predict(mylm, as.data.frame(cbind(
                IC$myallCV[, -1], GD2
              ))))
              Pred = cbind(
                as.character(DP$GD[, 1]),
                Group,
                RefInf,
                ID,
                BLUP,
                PEV,
                BLUE,
                predict(mylm, as.data.frame(cbind(
                  IC$myallCV[, -1], GD2
                )))
              )
              colnames(Pred) = c("Taxa",
                                 "Group",
                                 "RefInf",
                                 "ID",
                                 "BLUP",
                                 "PEV",
                                 "BLUE",
                                 "Prediction")
              
            }
            
            # print(dim(CV1))
            # print(table(index))
            print("Linear Regression to Predict phenotype Done !!")
            
          } else{
            print("aBLUP to Predict phenotype !!")
            
            if (!is.null(IC$myallCV))
            {
              if (!is.null(myBus$seqQTN))
              {
                busCV = cbind(IC$myallCV, X[, myBus$seqQTN])
              } else{
                numMarker = nrow(GWAS)
                sp = sort(GWAS$P.value)
                spd = abs(DP$cutOff - sp)
                index_fdr = grep(min(spd), spd)[1]
                FDRcutoff = DP$cutOff * index_fdr / numMarker
                seqQTN = as.numeric(rownames(GWAS[GWAS$P.value < FDRcutoff, ]))
                busCV = cbind(IC$myallCV, X[, seqQTN])
              }
              
            } else{
              busCV = cbind(as.data.frame(DP$GD[, 1]), X[, myBus$seqQTN])
            }
            pv = GWAS$P.value
            noneff = as.numeric(rownames(GWAS[GWAS$P.value > DP$cutOff, ]))
            
            if (is.null(DP$KI))
            {
              KI = GAPIT.kinship.VanRaden(snps = as.matrix(X))
              colnames(KI) = as.character(DP$GD[, 1])
              busKI = cbind(as.data.frame(DP$GD[, 1]), KI)
              colnames(busKI)[1] = c("Taxa")
            } else{
              busKI = DP$KI
            }
            print("The dimension of CV in lm model :")
            print(dim(busCV))
            # print(dim(busKI))
            # print(busKI[1:10,1:10])
            # print(cor(busCV[,-1]))
            busGAPIT = GAPIT(
              Y = ic_Y,
              K = busKI,
              CV = busCV,
              model = "gBLUP",
              file.output = F
            )
            Pred = busGAPIT$Pred
            print("aBLUP to Predict phenotype Done!!")
            
          }#lmpred
        }#buspred
        
        if (DP$file.output)
          write.csv(
            Pred,
            paste(
              "GAPIT.",
              DP$kinship.algorithm,
              ".Pred.result.csv",
              sep = ""
            ),
            row.names = FALSE,
            col.names = TRUE
          )
        
        
        va = myBus$vg
        ve = myBus$ve
        h2 = va / (va + ve)
        mc = NULL
        #mc=(exp(1)^(1/GWAS$P.value))/10000
        bc = NULL
        mp = NULL
        #myP=1/(exp(10000*fm$tau2)
        #print(str(GWAS))
        TV = NULL
        Compression = NULL
        GVs = myBus$GVs
      }
      #print(ic_GD[1:10,1:10])
      
      
      
      if (!DP$kinship.algorithm %in% c("FarmCPU", "MLMM", "Blink", "BlinkC"))
      {
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "GAPIT.Main")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "GAPIT.Main")
        
        GT = as.matrix(ic_GD[, 1])
        #print("!!!!!!!")
        #print(DP$sangwich.top)
        if (DP$PCA.total == 0)
          ic_PCA = NULL
        # print(ic_Y)
        #print(dim(ic_PCA))
        gapitMain <-
          GAPIT.Main(
            Y = ic_Y,
            GD = DP$GD[, -1],
            GM = DP$GM,
            KI = ic_KI,
            CV = DP$CV,
            CV.Inheritance = DP$CV.Inheritance,
            GP = DP$GP,
            GK = DP$GK,
            SNP.P3D = DP$SNP.P3D,
            kinship.algorithm = DP$kinship.algorithm,
            bin.from = DP$bin.from,
            bin.to = DP$bin.to,
            bin.by = DP$bin.by,
            inclosure.from = DP$inclosure.from,
            inclosure.to = DP$inclosure.to,
            inclosure.by = DP$inclosure.by,
            group.from = DP$group.from,
            group.to = DP$group.to,
            group.by = DP$group.by,
            kinship.cluster = DP$kinship.cluster,
            kinship.group = DP$kinship.group,
            name.of.trait = DP$name.of.trait,
            file.path = DP$file.path,
            file.from = DP$file.from,
            file.to = DP$file.to,
            file.total = DP$file.total,
            file.fragment = DP$file.fragment,
            file.G = DP$file.G,
            file.Ext.G = DP$file.Ext.G,
            file.GD = DP$file.GD,
            file.GM = DP$file.GM,
            file.Ext.GD = DP$file.Ext.GD,
            file.Ext.GM = DP$file.Ext.GM,
            SNP.MAF = DP$SNP.MAF,
            FDR.Rate = DP$FDR.Rate,
            SNP.FDR = DP$SNP.FDR,
            SNP.effect = DP$SNP.effect,
            SNP.impute = DP$SNP.impute,
            PCA.total = DP$PCA.total,
            GAPIT.Version = GAPIT.Version,
            GT = DP$GT,
            SNP.fraction = DP$SNP.fraction,
            seed = DP$seed,
            BINS = DP$BINS,
            SNP.test = DP$SNP.test,
            DPP = DP$DPP,
            SNP.permutation = DP$SNP.permutation,
            LD.chromosome = DP$LD.chromosome,
            LD.location = LD.location,
            LD.range = LD.range,
            SNP.CV = SNP.CV,
            SNP.robust = DP$SNP.robust,
            model = DP$model,
            genoFormat = "EMMA",
            hasGenotype = TRUE,
            byFile = FALSE,
            fullGD = TRUE,
            PC = DP$PC,
            GI = ic_GM,
            Timmer = DP$Timmer,
            Memory = DP$Memory,
            sangwich.top = DP$sangwich.top,
            sangwich.bottom = DP$sangwich.bottom,
            QC = DP$QC,
            GTindex = DP$GTindex,
            LD = DP$LD,
            file.output = FALSE,
            cutOff = DP$cutOff,
            GAPIT3.output = DP$file.output,
            Model.selection = DP$Model.selection,
            Create.indicator = DP$Create.indicator,
            QTN = DP$QTN,
            QTN.round = DP$QTN.round,
            QTN.limit = DP$QTN.limit,
            QTN.update = QTN.update,
            QTN.method = DP$QTN.method,
            Major.allele.zero = DP$Major.allele.zero,
            NJtree.group = DP$NJtree.group,
            NJtree.type = DP$NJtree.type,
            plot.bin = DP$plot.bin,
            QTN.position = DP$QTN.position,
            plot.style = DP$plot.style,
            SUPER_GS = DP$SUPER_GS
          )
        #print(str(gapitMain))
        GWAS = gapitMain$GWAS
        if (DP$Random.model)
          GR = GAPIT.RandomModel(
            Y = ic_Y,
            X = DP$GD[, -1],
            GWAS = GWAS,
            CV = gapitMain$PC,
            cutOff = DP$cutOff,
            GT = IC$GT
          )
        Pred = gapitMain$Pred
        #print(head(Pred))
        va = NA#gapitMain$vg
        ve = NA#gapitMain$ve
        h2 = gapitMain$h2
        mc = gapitMain$effect.snp
        bc = gapitMain$effect.cv
        mp = gapitMain$P
        TV = gapitMain$TV
        Compression = gapitMain$Compression
        GVs = GR$GVs
        
      }
      myPower = NULL
      #print(head(GWAS))
      #print(DP$QTN.position)
      if (!is.null(GWAS))
        myPower = GAPIT.Power(
          WS = DP$WS,
          alpha = DP$alpha,
          maxOut = DP$maxOut,
          seqQTN = DP$QTN.position,
          GM = DP$GM,
          GWAS = GWAS
        )
      
      #print(str(myPower))
      #print("GAPIT.III accomplished successfully for multiple traits. Results are saved")
      return (
        list(
          GWAS = GWAS,
          Pred = Pred,
          FDR = myPower$FDR,
          Power = myPower$Power,
          Power.Alpha = myPower$Power.Alpha,
          alpha = myPower$alpha,
          h2 = h2,
          va = va,
          ve = ve,
          mc = mc,
          bc = bc,
          mp = mp,
          TV = TV,
          Compression = Compression,
          Timmer = Timmer,
          Memory = Memory,
          GVs = GVs
        )
      )
    } else{
      # Here is Genomic Prediction function
      
      gapitMain <-
        GAPIT.Main(
          Y = IC$Y,
          GD = DP$GD[, -1],
          GM = DP$GM,
          KI = DP$KI,
          Z = DP$Z,
          CV = DP$CV,
          CV.Inheritance = DP$CV.Inheritance,
          GP = DP$GP,
          GK = DP$GK,
          SNP.P3D = DP$SNP.P3D,
          kinship.algorithm = DP$kinship.algorithm,
          bin.from = DP$bin.from,
          bin.to = DP$bin.to,
          bin.by = DP$bin.by,
          inclosure.from = DP$inclosure.from,
          inclosure.to = DP$inclosure.to,
          inclosure.by = DP$inclosure.by,
          group.from = DP$group.from,
          group.to = DP$group.to,
          group.by = DP$group.by,
          kinship.cluster = DP$kinship.cluster,
          kinship.group = DP$kinship.group,
          name.of.trait = DP$name.of.trait,
          file.path = DP$file.path,
          file.from = DP$file.from,
          file.to = DP$file.to,
          file.total = DP$file.total,
          file.fragment = DP$file.fragment,
          file.G = DP$file.G,
          file.Ext.G = DP$file.Ext.G,
          file.GD = DP$file.GD,
          file.GM = DP$file.GM,
          file.Ext.GD = DP$file.Ext.GD,
          file.Ext.GM = DP$file.Ext.GM,
          SNP.MAF = DP$SNP.MAF,
          FDR.Rate = DP$FDR.Rate,
          SNP.FDR = DP$SNP.FDR,
          SNP.effect = DP$SNP.effect,
          SNP.impute = DP$SNP.impute,
          PCA.total = DP$PCA.total,
          GAPIT.Version = GAPIT.Version,
          GT = DP$GT,
          SNP.fraction = DP$SNP.fraction,
          seed = DP$seed,
          BINS = DP$BINS,
          SNP.test = DP$SNP.test,
          DPP = DP$DPP,
          SNP.permutation = DP$SNP.permutation,
          LD.chromosome = DP$LD.chromosome,
          LD.location = LD.location,
          LD.range = LD.range,
          SNP.CV = SNP.CV,
          SNP.robust = DP$SNP.robust,
          model = DP$model,
          genoFormat = "EMMA",
          hasGenotype = TRUE,
          byFile = FALSE,
          fullGD = TRUE,
          PC = DP$PC,
          GI = DP$GI,
          Timmer = DP$Timmer,
          Memory = DP$Memory,
          GAPIT3.output = DP$file.output,
          sangwich.top = DP$sangwich.top,
          sangwich.bottom = DP$sangwich.bottom,
          QC = DP$QC,
          GTindex = DP$GTindex,
          LD = DP$LD,
          file.output = FALSE,
          cutOff = DP$cutOff,
          Model.selection = DP$Model.selection,
          Create.indicator = DP$Create.indicator,
          QTN = DP$QTN,
          QTN.round = DP$QTN.round,
          QTN.limit = DP$QTN.limit,
          QTN.update = QTN.update,
          QTN.method = DP$QTN.method,
          Major.allele.zero = DP$Major.allele.zero,
          NJtree.group = DP$NJtree.group,
          NJtree.type = DP$NJtree.type,
          plot.bin = DP$plot.bin,
          QTN.position = DP$QTN.position,
          plot.style = DP$plot.style,
          SUPER_GS = DP$SUPER_GS
        )
      #print(str(gapitMain))
      GWAS = gapitMain$GWAS
      Pred = gapitMain$Pred
      #print(head(Pred))
      va = NA#gapitMain$vg
      ve = NA#gapitMain$ve
      h2 = gapitMain$h2
      mc = gapitMain$effect.snp
      bc = gapitMain$effect.cv
      mp = gapitMain$P
      Compression = gapitMain$Compression
      GAPIT.Compression.Visualization(Compression = Compression,
                                      name.of.trait = DP$name.of.trait)
      # # print(list(GWAS=GWAS,Pred=Pred,FDR=NULL,Power=NULL,
      #   Power.Alpha=NULL,alpha=NULL,h2=h2,va=va,ve=ve,Compression=Compression,
      #   mc=mc,bc=bc,mp=mp,TV=gapitMain$TV,
      #   Timmer=Timmer,Memory=Memory))
      return (
        list(
          GWAS = GWAS,
          Pred = Pred,
          FDR = NULL,
          Power = NULL,
          Power.Alpha = NULL,
          alpha = NULL,
          h2 = h2,
          va = va,
          ve = ve,
          Compression = Compression,
          mc = mc,
          bc = bc,
          mp = mp,
          TV = gapitMain$TV,
          Timmer = Timmer,
          Memory = Memory
        )
      )
    }#end of SNP.TEST
    
  }  #end of GAPIT.SS function

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.Main` <-
  function(Y,
           G = NULL,
           GD = NULL,
           GM = NULL,
           KI = NULL,
           Z = NULL,
           CV = NULL,
           CV.Inheritance = NULL,
           SNP.P3D = TRUE,
           GP = NULL,
           GK = NULL,
           group.from = 1000000 ,
           group.to = 1,
           group.by = 10,
           kinship.cluster = "average",
           kinship.group = 'Mean',
           kinship.algorithm = NULL,
           DPP = 50000,
           ngrid = 100,
           llin = -10,
           ulim = 10,
           esp = 1e-10,
           GAPIT3.output = TRUE,
           file.path = NULL,
           file.from = NULL,
           file.to = NULL,
           file.total = NULL,
           file.fragment = 512,
           file.G = NULL,
           file.Ext.G = NULL,
           file.GD = NULL,
           file.GM = NULL,
           file.Ext.GD = NULL,
           file.Ext.GM = NULL,
           SNP.MAF = 0,
           FDR.Rate = 1,
           SNP.FDR = 1,
           SNP.effect = "Add",
           SNP.impute = "Middle",
           PCA.total = 0,
           GAPIT.Version = GAPIT.Version,
           name.of.trait,
           GT = NULL,
           SNP.fraction = 1,
           seed = 123,
           BINS = 20,
           SNP.test = TRUE,
           SNP.robust = "FaST",
           LD.chromosome = NULL,
           LD.location = NULL,
           LD.range = NULL,
           model = model,
           bin.from = 10000,
           bin.to = 5000000,
           bin.by = 1000,
           inclosure.from = 10,
           inclosure.to = 1000,
           inclosure.by = 10,
           SNP.permutation = FALSE,
           SNP.CV = NULL,
           NJtree.group = NJtree.group,
           NJtree.type = NJtree.type,
           plot.bin = plot.bin,
           genoFormat = NULL,
           hasGenotype = NULL,
           byFile = NULL,
           fullGD = NULL,
           PC = NULL,
           GI = NULL,
           Timmer = NULL,
           Memory = NULL,
           sangwich.top = NULL,
           sangwich.bottom = NULL,
           QC = TRUE,
           GTindex = NULL,
           LD = 0.05,
           file.output = TRUE,
           cutOff = 0.05,
           Model.selection = FALSE,
           Create.indicator = FALSE,
           QTN = NULL,
           QTN.round = 1,
           QTN.limit = 0,
           QTN.update = TRUE,
           QTN.method = "Penalty",
           Major.allele.zero = FALSE,
           QTN.position = NULL,
           SUPER_GD = NULL,
           SUPER_GS = SUPER_GS,
           plot.style = "Beach",
           CG = CG,
           chor_taxa = chor_taxa) {
    #Object: To perform GWAS and GPS (Genomic Prediction or Selection)
    #Output: GWAS table (text file), QQ plot (PDF), Manhattan plot (PDF), genomic prediction (text file), and
    #        genetic and residual variance components
    #Authors: Zhiwu Zhang
    # Last update: Oct 23, 2015  by Jiabo Wang add REML threshold and SUPER GD KI
    ##############################################################################################
    
    #Initial p3d and h2.opt temporaryly
    h2.opt = NULL
    p3d = list(
      ps = NULL,
      REMLs = NULL,
      stats = NULL,
      effect.est = NULL,
      rsquare_base = NULL,
      rsquare = NULL,
      dfs = NULL,
      df = NULL,
      tvalue = NULL,
      stderr = NULL,
      maf = NULL,
      nobs = NULL,
      Timmer = NULL,
      Memory = NULL,
      vgs = NULL,
      ves = NULL,
      BLUP = NULL,
      BLUP_Plus_Mean = NULL,
      PEV = NULL,
      BLUE = NULL,
      logLM = NULL,
      effect.snp = NULL,
      effect.cv = NULL
    )
    
    
    if (SUPER_GS)
    {
      Compression = NULL
      kinship.optimum = NULL
      kinship = NULL
      PC = PC
      REMLs = NULL
      GWAS = NULL
      QTN = NULL
      Timmer = GAPIT.Timmer(Infor = "GAPIT.SUPER.GS")
      Memory = GAPIT.Memory(Infor = "GAPIT.SUPER.GS")
      #print(model)
      SUPER_GS_GAPIT = GAPIT.SUPER.GS(
        Y = Y,
        GD = GD,
        GM = GM,
        KI = KI,
        Z = Z,
        CV = CV,
        GK = GK,
        kinship.algorithm = kinship.algorithm,
        bin.from = bin.from,
        bin.to = bin.to,
        bin.by = bin.by,
        inclosure.from = inclosure.from,
        inclosure.to = inclosure.to,
        inclosure.by = inclosure.by,
        group.from = group.from,
        group.to = group.to,
        group.by = group.by,
        kinship.cluster = kinship.cluster,
        kinship.group = kinship.group,
        PCA.total = PCA.total,
        GT = GT,
        PC = PC,
        GI = GI,
        Timmer = Timmer,
        Memory = Memory,
        model = model,
        sangwich.top = sangwich.top,
        sangwich.bottom = sangwich.bottom,
        QC = QC,
        GTindex = GTindex,
        LD = LD,
        file.output = GAPIT3.output,
        cutOff = cutOff
      )
      # Compression=as.matrix(SUPER_GS_GAPIT$Compression)
      # opt=
      print("SUPER_GS_GAPIT FUNCTION DONE")
      return (
        list(
          Compression = SUPER_GS_GAPIT$Compression,
          kinship.optimum = SUPER_GS_GAPIT$SUPER_kinship,
          kinship = SUPER_GS_GAPIT$kinship,
          PC = SUPER_GS_GAPIT$PC,
          GWAS = GWAS,
          GPS = SUPER_GS_GAPIT$GPS,
          Pred = SUPER_GS_GAPIT$Pred,
          Timmer = Timmer,
          Memory = Memory,
          h2 = SUPER_GS_GAPIT$h2,
          SUPER_GD = SUPER_GS_GAPIT$SUPER_GD,
          GWAS = NULL,
          QTN = NULL
        )
      )
      
    } else{
      #print("@@@@@@@")
      #print(group.from)
      
      #Handler of SNP.test=F
      #Iniciate with two by seven NA matrix
      #The seventh is for p values of SNP
      DTS = rbind(rep(NA, 7), rep(NA, 7))
      
      
      #End imediatly in one of these situtiona
      shortcut = FALSE
      LL.save = 1e10
      #In case of null Y and null GP, sent back genotype only
      thisY = Y[, 2]
      thisY = thisY[!is.na(thisY)]
      if (length(thisY) < 3) {
        shortcut = TRUE
      } else{
        if (var(thisY) == 0)
          shortcut = TRUE
      }
      
      if (shortcut) {
        print(paste("Y is empty. No GWAS/GS performed for ", name.of.trait, sep =
                      ""))
        return (
          list(
            compression = NULL,
            kinship.optimum = NULL,
            kinship = KI,
            PC = PC,
            GWAS = NULL,
            GPS = NULL,
            Pred = NULL,
            REMLs = NULL,
            Timmer = Timmer,
            Memory = Memory,
            h2 = NULL
          )
        )
      }
      
      #QC
      print("------------Examining data (QC)------------------------------------------")
      if (is.null(Y))
        stop ("GAPIT says: Phenotypes must exist.")
      if (is.null(KI) &
          missing(GD) &
          kinship.algorithm != "SUPER")
        stop (
          "GAPIT says: Kinship is required. As genotype is not provided, kinship can not be created."
        )
      
      #When GT and GD are missing, force to have fake ones (creating them from Y),GI is not required in this case
      if (is.null(GD) & is.null(GT)) {
        GT = as.matrix(Y[, 1])
        GD = matrix(1, nrow(Y), 1)
        GI = as.data.frame(matrix(0, 1, 3))
        colnames(GI) = c("SNP", "Chromosome", "Position")
      }
      
      if (is.null(GT)) {
        GT = as.character(Y[, 1])
      }
      #print("@@@@@@@@")
      #print(GD)
      #merge CV with PC: Put CV infront of PC
      if (PCA.total > 0 & !is.null(CV))
        CV = GAPIT.CVMergePC(CV, PC)
      if (PCA.total > 0 & is.null(CV))
        CV = PC
      #for GS merge CV with GD name
      if (is.null(CV))
      {
        my_allCV = CV
      } else{
        taxa_GD = rownames(GD)
        
        my_allCV = CV[order(CV[, 1]), ]
        my_allCV = my_allCV[my_allCV[, 1] %in% taxa_GD, ]
        #print(dim(my_allCV))
      }
      
      #Handler of CV.Inheritance
      if (is.null(CV) & !is.null(CV.Inheritance)) {
        stop ("GAPIT says: CV.Inheritance is more than avaiable.")
      }
      
      if (!is.null(CV) & !is.null(CV.Inheritance)) {
        if (CV.Inheritance > (ncol(CV) - 1))
          stop ("GAPIT says: CV.Inheritance is more than avaiable.")
      }
      
      #Create Z as identity matrix from Y if it is not provided
      if (kinship.algorithm != "None" &
          kinship.algorithm != "SUPER" & is.null(Z)) {
        taxa = as.character(Y[, 1]) #this part will make GS without CV not present all prediction
        Z = as.data.frame(diag(1, nrow(Y)))
        #taxa=as.character(KI[,1])
        #Z=as.data.frame(diag(1,nrow(KI)))
        Z = rbind(taxa, Z)
        taxa = c('Taxa', as.character(taxa))
        Z = cbind(taxa, Z)
      }
      ZI = Z
      
      #Add the part of non proportion in Z matrix
      if (kinship.algorithm != "None" &
          kinship.algorithm != "SUPER" & !is.null(Z))
      {
        if (nrow(Z) - 1 < nrow(Y))
          Z = GAPIT.ZmatrixFormation(Z = Z, Y = Y)
      }
      
      #Create CV with all 1's if it is not provided
      noCV = FALSE
      if (is.null(CV)) {
        noCV = TRUE
        CV = Y[, 1:2]
        CV[, 2] = 1
        colnames(CV) = c("taxa", "overall")
      }
      
      #Remove duplicat and integragation of data
      print("QC is in process...")
      
      CVI <- CV
      
      # print(dim(Z))
      if (QC)
      {
        qc <- GAPIT.QC(
          Y = Y,
          KI = KI,
          GT = GT,
          CV = CV,
          Z = Z,
          GK = GK
        )
        GTindex = qc$GTindex
        Y = qc$Y
        KI = qc$KI
        CV = qc$CV
        Z = qc$Z
        GK = qc$GK
        if (noCV)
          CVI = qc$CV #this part will make GS without CV not present all prediction
        my_taxa = as.character(KI[, 1])
      }
      #print(GTindex)
      
      #print(dim(KI))
      #Output phenotype
      colnames(Y) = c("Taxa", name.of.trait)
      if (file.output)
      {
        try(write.table(
          Y,
          paste("GAPIT.", name.of.trait, ".phenotype.csv" , sep = ""),
          quote = FALSE,
          sep = ",",
          row.names = FALSE,
          col.names = TRUE
        ))
      }
      #TDP
      if (kinship.algorithm == "None")
      {
        if (min(CV[, 2]) == max(CV[, 2]))
          CV = NULL
        
        theTDP = GAPIT.TDP(
          Y = Y,
          CV = CV,
          SNP = as.data.frame(cbind(GT[GTindex], as.matrix(
            as.data.frame(GD[GTindex, ])
          ))),
          QTN = QTN,
          Round = QTN.round,
          QTN.limit = QTN.limit,
          QTN.update = QTN.update,
          Method = QTN.method
        )
        #print(dim(GM))
        #print(length(theTDP$p))
        
        theGWAS = cbind(GM, theTDP$p, NA, NA, NA)
        
        return (
          list(
            Compression = NULL,
            kinship.optimum = NULL,
            kinship = NULL,
            PC = NULL,
            GWAS = theGWAS,
            GPS = NULL,
            Pred = NULL,
            REMLs = NULL,
            QTN = theTDP$QTN,
            Timmer = Timmer,
            Memory = Memory,
            h2 = NULL
          )
        )
      }
      
      rm(qc)
      gc()
      
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "QC")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "QC")
      
      #Get indicator of sangwich top and bottom
      byPass.top = FALSE
      byPass = FALSE
      NOBLUP = FALSE
      if (group.from < 2 & group.to < 2)
        NOBLUP = TRUE
      #if(!is.null(sangwich.bottom)) byPass=((sangwich.bottom=="FaST" | sangwich.bottom=="SUPER" | sangwich.bottom=="DC" )& is.null(GP)   )
      if (!is.null(sangwich.top))
        byPass.top = ((
          sangwich.top == "FaST" |
            sangwich.top == "SUPER" | sangwich.top == "DC"
        ))
      if (!is.null(sangwich.bottom))
        byPass = ((
          sangwich.bottom == "FaST" |
            sangwich.bottom == "SUPER" |
            sangwich.bottom == "DC"
        ))
      
      print("Try to group from and to were set to 1")
      
      if (byPass) {
        print("group from and to were set to 1")
        group.from = 1
        group.to = 1
      }
      
      print("------------Examining data (QC) done-------------------------------------")
      
      #Sagnwich top bun: To gep GP if it is not provided
      if (!is.null(sangwich.top) & is.null(GP))
      {
        print("-------------------Sandwich top bun-----------------------------------")
        #print(dim(GD))
        #print(GD[1:5,1:5])
        
        #Create GK if not provided
        if (is.null(GK)) {
          #    set.seed(1)
          nY = floor(nrow(Y) * .9)
          nG = ncol(GD)
          if (nG > nY) {
            snpsam = sample(1:nG, nY)
          } else{
            snpsam = 1:nG
          }
          GK = GD[GTindex, snpsam]
          SNPVar = apply(as.matrix(GK), 2, var)
          GK = GK[, SNPVar > 0]
          GK = cbind(as.data.frame(GT[GTindex]), as.data.frame(GK)) #add taxa
          
        }
        
        #myGD=cbind(as.data.frame(GT),as.data.frame(GD))
        file.output.temp = file.output
        file.output = FALSE
        #print(sangwich.top)
        GP = GAPIT.Bread(
          Y = Y,
          CV = CV,
          Z = Z,
          KI = KI,
          GK = GK,
          GD = cbind(as.data.frame(GT), as.data.frame(GD)),
          GM = GI,
          method = sangwich.top,
          GTindex = GTindex,
          LD = LD,
          file.output = file.output
        )$GWAS
        file.output = file.output.temp
        
        
        
        GK = NULL
        
        print("-------------------Sagnwich top bun: done-----------------------------")
        
      }
      
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "SagnwichTop")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "SagnwichTop")
      
      #Sandwich burger and dressing
      print("-------------------Sandwich burger and dressing------------------------")
      
      #Handler of group boundry
      if (group.from > group.to)
        stop("GAPIT says: group.to should  be larger than group.from. Please correct them!")
      
      if (is.null(CV) | (!is.null(CV) & group.to < (ncol(CV) + 1))) {
        #The minimum of group is 1 + number of columns in CV
        group.from = 1
        group.to = 1
        message(
          "The upper bound of groups (group.to) is not sufficient. both boundries were set to a and GLM is performed!"
        )
      }
      
      if (!is.null(CV) & group.from < 1) {
        group.from = 1 #minimum of group is number of columns in CV
        warning("The lower bound of groups should be 1 at least. It was set to 1!")
      }
      
      nk = 1000000000
      if (!is.null(KI))
        nk = min(nk, nrow(KI))
      if (!is.null(GK))
        nk = min(nk, nrow(GK))
      
      if (!is.null(KI))
      {
        if (group.to > nk) {
          #group.to=min(nrow(KI),length(GTindex)) #maximum of group is number of rows in KI
          group.to = nk #maximum of group is number of rows in KI
          warning("The upper bound of groups is too high. It was set to the size of kinship!")
        }
        if (group.from > nk) {
          group.from = nk
          warning("The lower bound of groups is too high. It was set to the size of kinship!")
        }
      }
      
      if (!is.null(CV)) {
        if (group.to <= ncol(CV) + 1) {
          #The minimum of group is number of columns in CV
          #group.from=ncol(CV)+2
          #group.to=ncol(CV)+2
          message(
            "The upper bound of groups (group.to) is not sufficient. both boundries were set to their minimum and GLM is performed!"
          )
        }
      }
      
      #bin.fold=ceiling(log2(bin.to/bin.from))
      #bin.seq=0:bin.fold
      #bin.level=bin.from*2^bin.seq
      
      #Set upper bound for inclosure.to
      if (inclosure.to > nrow(Y))
        inclosure.to = nrow(Y) - 1
      
      #set inclosure loop levels
      bin.level = seq(bin.from, bin.to, by = bin.by)
      inclosure = seq(inclosure.from, inclosure.to, by = inclosure.by)
      
      #Optimization for group number, cluster algorithm and kinship type
      GROUP = seq(group.to, group.from, by = -group.by)#The reverse order is to make sure to include full model
      if (missing("kinship.cluster"))
        kinship.cluster = c("ward",
                            "single",
                            "complete",
                            "average",
                            "mcquitty",
                            "median",
                            "centroid")
      if (missing("kinship.group"))
        kinship.group = c("Mean", "Max", "Min", "Median")
      numSetting = length(GROUP) * length(kinship.cluster) * length(kinship.group) *
        length(bin.level) * length(inclosure)
      
      #Reform Y, GD and CV into EMMA format
      ys = as.matrix(Y[, 2])
      X0 = as.matrix(CV[, -1])
      CV.taxa = CVI[, 1]
      #print(length(ys))
      #Initial
      count = 0
      Compression = matrix(, numSetting, 6)
      colnames(Compression) = c("Type", "Cluster", "Group", "REML", "VA", "VE")
      
      #add indicator of overall mean
      if (min(X0[, 1]) != max(X0[, 1]))
        X0 <-
        cbind(1, X0) #do not add overall mean if X0 has it already at first column
      
      
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "DataProcessing")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "DataProcessing")
      
      print("-------------------------Iteration in process--------------------------")
      print(paste("Total iterations: ", numSetting, sep = ""))
      
      #Loop to optimize cluster algorithm, group number and kinship type
      for (bin in bin.level) {
        for (inc in inclosure) {
          #Grill: update KI if GK or GP is provided
          if (!byPass & (!is.null(GK) | !is.null(GP)))
          {
            print("Grilling KI...")
            
            myGenotype <-
              GAPIT.Genotype(
                G = NULL,
                GD = cbind(as.data.frame(GT), as.data.frame(GD)),
                GM = GI,
                KI = NULL,
                kinship.algorithm = kinship.algorithm,
                PCA.total = 0,
                SNP.fraction = SNP.fraction,
                SNP.test = SNP.test,
                file.path = file.path,
                file.from = file.from,
                file.to = file.to,
                file.total = file.total,
                file.fragment = file.fragment,
                file.G = file.G,
                file.Ext.G = file.Ext.G,
                file.GD = file.GD,
                file.GM = file.GM,
                file.Ext.GD = file.Ext.GD,
                file.Ext.GM = file.Ext.GM,
                SNP.MAF = SNP.MAF,
                FDR.Rate = FDR.Rate,
                SNP.FDR = SNP.FDR,
                SNP.effect = SNP.effect,
                SNP.impute = SNP.impute,
                kinship.cluster = kinship.cluster,
                NJtree.group = NJtree.group,
                NJtree.type = NJtree.type,
                LD.chromosome = LD.chromosome,
                LD.location = LD.location,
                LD.range = LD.range,
                GP = GP,
                GK = GK,
                bin.size = bin,
                inclosure.size = inc,
                SNP.CV = SNP.CV,
                Timmer = Timmer,
                Memory = Memory,
                GTindex = GTindex,
                sangwich.top = NULL,
                sangwich.bottom = sangwich.bottom,
                file.output = file.output,
                Create.indicator = Create.indicator,
                Major.allele.zero = Major.allele.zero
              )
            
            Timmer = myGenotype$Timmer
            Memory = myGenotype$Memory
            
            KI = myGenotype$KI
            #update group set by new KI
            nk = nrow(KI)
            GROUP = GROUP[GROUP <= nk]
          }
          for (ca in kinship.cluster) {
            for (group in GROUP) {
              for (kt in kinship.group) {
                #Do not screen SNP unless existing genotype and one combination
                if (numSetting == 1 & hasGenotype) {
                  optOnly = FALSE
                } else{
                  optOnly = TRUE
                }
                if (!SNP.test)
                  optOnly = TRUE
                
                if (optOnly | Model.selection) {
                  colInclude = 1
                  optOnly = TRUE
                } else{
                  colInclude = c(1:ncol(GD))
                }
                
                if (!optOnly) {
                  print("Compressing and Genome screening...")
                }
                count = count + 1
                
                #Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PreP3D 1")
                #Memory=GAPIT.Memory(Memory=Memory,Infor="PreP3D 1")
                
                if (!byPass)
                {
                  if (count == 1)
                    print("-------Mixed model with Kinship-----------------------------")
                  if (group < ncol(X0) + 1)
                    group = 1 # the emma function (emma.delta.REML.dLL.w.Z) does not allow K has dim less then CV. turn to GLM (group=1)
                  
                  cp <-
                    GAPIT.Compress(
                      KI = KI,
                      kinship.cluster = ca,
                      kinship.group = kt,
                      GN = group,
                      Timmer = Timmer,
                      Memory = Memory
                    )
                  Timmer = cp$Timmer
                  Memory = cp$Memory
                  
                  Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "PreP3D 2_cp")
                  Memory = GAPIT.Memory(Memory = Memory, Infor = "PreP3D 2_cp")
                  
                  #print("BK...")
                  
                  bk <- GAPIT.Block(Z = Z,
                                    GA = cp$GA,
                                    KG = cp$KG)
                  
                  Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "PreP3D 2_bk")
                  Memory = GAPIT.Memory(Memory = Memory, Infor = "PreP3D 2 bk")
                  
                  #print("ZC...")
                  zc <- GAPIT.ZmatrixCompress(Z = Z, GAU = bk$GA)
                  
                  Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "PreP3D 2_zc")
                  Memory = GAPIT.Memory(Memory = Memory, Infor = "PreP3D 2 zc")
                  
                  #print("wraping...")
                  #Reform KW and Z into EMMA format
                  
                  zrow = nrow(zc$Z)
                  zcol = ncol(zc$Z) - 1
                  #Z1=matrix(as.numeric(as.matrix(zc$Z[,-1])),nrow=zrow,ncol=zcol)
                  
                  Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Prio PreP3D")
                  Memory = GAPIT.Memory(Memory = Memory, Infor = "Prio PreP3D")
                  
                  #Evaluating maximum likelohood
                  #print("Calling EMMAxP3D...")
                  
                  #print("It made it to here")
                  #print("The dimension of xs is:")
                  #print("The value of SNP.impute is")
                  #print(SNP.impute)
                  
                  #write.table(zc$Z, "Z.csv", quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
                  
                  # print(head(ys))
                  # print(as.matrix(as.data.frame(GD[GTindex,colInclude]))[1:5,1:5])
                  # print(as.matrix(bk$KW)[1:5,1:5])
                  # print(ys)
                  # # print(as.matrix(as.data.frame(GD[GTindex,colInclude])))
                  # print(dim(as.matrix(bk$KW)))
                  # print(as.matrix(bk$KW)[,1])
                  # print(dim(as.matrix(as.data.frame(GD[GTindex,colInclude]))))
                  
                  p3d <-
                    GAPIT.EMMAxP3D(
                      ys = ys,
                      xs = as.matrix(as.data.frame(GD[GTindex, colInclude])),
                      K = as.matrix(bk$KW) ,
                      Z = matrix(
                        as.numeric(as.matrix(zc$Z[, -1])),
                        nrow = zrow,
                        ncol = zcol
                      ),
                      X0 = X0,
                      CVI = CVI,
                      CV.Inheritance = CV.Inheritance,
                      GI = GI,
                      SNP.P3D = SNP.P3D,
                      Timmer = Timmer,
                      Memory = Memory,
                      fullGD = fullGD,
                      SNP.permutation = SNP.permutation,
                      GP = GP,
                      SNP.fraction = SNP.fraction,
                      file.path = file.path,
                      file.from = file.from,
                      file.to = file.to,
                      file.total = file.total,
                      file.fragment = file.fragment,
                      byFile = byFile,
                      file.G = file.G,
                      file.Ext.G = file.Ext.G,
                      file.GD = file.GD,
                      file.GM = file.GM,
                      file.Ext.GD = file.Ext.GD,
                      file.Ext.GM = file.Ext.GM,
                      GTindex = GTindex,
                      genoFormat = genoFormat,
                      optOnly = optOnly,
                      SNP.effect = SNP.effect,
                      SNP.impute = SNP.impute,
                      name.of.trait = name.of.trait,
                      Create.indicator = Create.indicator,
                      Major.allele.zero = Major.allele.zero
                    )
                  
                  Timmer = p3d$Timmer
                  Memory = p3d$Memory
                  
                  Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Post PreP3D")
                  Memory = GAPIT.Memory(Memory = Memory, Infor = "Post PreP3D")
                  
                  #print("Cluster algorithm, kinship type, groups, VG, Ve and REML:")
                  print(
                    paste(
                      count,
                      "of",
                      numSetting,
                      "--",
                      "Vg=",
                      round(p3d$vgs, 4),
                      "VE=",
                      round(p3d$ves, 4),
                      "-2LL=",
                      round(p3d$REMLs, 2),
                      "  Clustering=",
                      ca,
                      "  Group number=",
                      group ,
                      "  Group kinship=",
                      kt,
                      sep = " "
                    )
                  )
                  #print(table(GTindex))
                  
                  #Recoding the optimum KI
                  if (count == 1) {
                    KI.save = KI
                    LL.save = p3d$REMLs
                  } else{
                    if (p3d$REMLs < LL.save) {
                      KI.save = KI
                      LL.save = p3d$REMLs
                    }
                  }
                  
                  #print(paste("CA is ",ca))
                  #print(paste("group is ",group))
                  #print(paste("kt is ",kt))
                  
                  #recording Compression profile on array
                  Compression[count, 1] = kt
                  Compression[count, 2] = ca
                  Compression[count, 3] = group
                  Compression[count, 4] = p3d$REMLs
                  Compression[count, 5] = p3d$vgs
                  Compression[count, 6] = p3d$ves
                  #print("result saved")
                  
                } else{
                  # end of if(!byPass)
                  
                  #Set QTNs
                  if (count == 1)
                    print("-------The burger is SNP-----------------------------------")
                  #bin.size=bin
                  #inclosure.size=inc
                  
                  
                  #@@@This section is not useful
                  if (!is.null(GP))
                  {
                    #print("Being specific...")
                    
                    myGenotype <-
                      GAPIT.Genotype(
                        G = NULL,
                        GD = NULL,
                        GM = GI,
                        KI = NULL,
                        kinship.algorithm = "SUPER",
                        PCA.total = 0,
                        SNP.fraction = SNP.fraction,
                        SNP.test = SNP.test,
                        file.path = file.path,
                        file.from = file.from,
                        file.to = file.to,
                        file.total = file.total,
                        file.fragment = file.fragment,
                        file.G = file.G,
                        file.Ext.G = file.Ext.G,
                        file.GD = file.GD,
                        file.GM = file.GM,
                        file.Ext.GD = file.Ext.GD,
                        file.Ext.GM = file.Ext.GM,
                        SNP.MAF = SNP.MAF,
                        FDR.Rate = FDR.Rate,
                        SNP.FDR = SNP.FDR,
                        SNP.effect = SNP.effect,
                        SNP.impute = SNP.impute,
                        LD.chromosome = LD.chromosome,
                        LD.location = LD.location,
                        LD.range = LD.range,
                        kinship.cluster = kinship.cluster,
                        #NJtree.group=NJtree.group,NJtree.type=NJtree.type,
                        GP = GP,
                        GK = NULL,
                        bin.size = bin,
                        inclosure.size = inc,
                        SNP.CV = SNP.CV,
                        GTindex = GTindex,
                        sangwich.top = NULL,
                        sangwich.bottom = sangwich.bottom,
                        Timmer = Timmer,
                        Memory = Memory,
                        file.output = file.output,
                        Create.indicator = Create.indicator,
                        Major.allele.zero = Major.allele.zero
                      )
                    
                    Timmer = myGenotype$Timmer
                    Memory = myGenotype$Memory
                    
                    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Genotype for burger")
                    Memory = GAPIT.Memory(Memory = Memory, Infor = "Genotype for burger")
                    
                    print(paste("bin---", bin, "---inc---", inc, sep = ""))
                    GK = GD[GTindex, myGenotype$SNP.QTN]
                    SUPER_GD = GD[, myGenotype$SNP.QTN]
                    SNPVar = apply(as.matrix(GK), 2, var)
                    
                    GK = GK[, SNPVar > 0]
                    SUPER_GD = SUPER_GD[, SNPVar > 0]
                    GK = cbind(as.data.frame(GT[GTindex]), as.data.frame(GK)) #add taxa
                    SUPER_GD = cbind(as.data.frame(GT),
                                     as.data.frame(SUPER_GD)) #add taxa
                    
                    #GP=NULL
                  }# end of if(is.null(GK))
                  
                  
                  if (!is.null(GK) & numSetting > 1)
                  {
                    print(
                      "-------Calculating likelihood-----------------------------------"
                    )
                    # myBurger=GAPIT.Burger(Y=Y,CV=CV,GK=GK)
                    myBurger = GAPIT.Burger(Y = Y,
                                            CV = NULL,
                                            GK = GK)   #########modified by Jiabo Wang
                    
                    myREML = myBurger$REMLs
                    myVG = myBurger$vg
                    myVE = myBurger$ve
                  } else{
                    myREML = NA
                    myVG = NA
                    myVE = NA
                  }
                  
                  #Recoding the optimum GK
                  if (count == 1) {
                    GK.save = GK
                    LL.save = myREML
                    SUPER_optimum_GD = SUPER_GD     ########### get SUPER GD
                    
                  } else{
                    if (myREML < LL.save) {
                      GK.save = GK
                      LL.save = myREML
                      SUPER_optimum_GD = SUPER_GD     ########### get SUPER GD
                    }
                  }
                  
                  #Put to storage
                  Compression[count, 1] = 1
                  Compression[count, 2] = bin
                  Compression[count, 3] = inc
                  Compression[count, 4] = myREML
                  Compression[count, 5] = myVG
                  Compression[count, 6] = myVG
                  print(Compression[count, ])
                  
                  #print("---------------SUPER 2nd stage: calculating LL ------------------------")
                  
                  
                }   # end of if(byPass)
                
              }#end of for (ca in kinship.cluster)
              
              #Skip the rest group in case group 1 is finished
              if (group == 1)
                break #To skip the rest group interations
              
            }#end of for (group in GROUP)
          }#end of for (kt in kinship.group)
          
          
        }#end of for (inc in inclosure)
      }#end of for (bin in bin.level)
      
      
      if (Model.selection == TRUE) {
        print(
          "------------------------Model selection for optimal number of PCs and Covariates-------------------------------------------------"
        )
        #update KI with the best likelihood
        KI = KI.save
        if (numSetting > 1) {
          Compression = Compression[order(as.numeric(Compression[, 4]), decreasing = FALSE), ]  #sort on REML
          kt = Compression[1, 1]
          ca = Compression[1, 2]
          group = Compression[1, 3]
        }
        
        cp <-
          GAPIT.Compress(
            KI = KI,
            kinship.cluster = ca,
            kinship.group = kt,
            GN = group,
            Timmer = Timmer,
            Memory = Memory
          )
        Timmer = cp$Timmer
        Memory = cp$Memory
        
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "PreP3D 2_cp")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "PreP3D 2_cp")
        
        bk <- GAPIT.Block(Z = Z,
                          GA = cp$GA,
                          KG = cp$KG)
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "PreP3D 2_bk")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "PreP3D 2 bk")
        
        zc <- GAPIT.ZmatrixCompress(Z = Z, GAU = bk$GA)
        
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "PreP3D 2_zc")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "PreP3D 2 zc")
        
        z0 = as.matrix(zc$Z[, -1])
        Z1 = matrix(as.numeric(z0),
                    nrow = nrow(z0),
                    ncol = ncol(z0))
        
        
        
        BIC <- rep(NA, ncol(X0))
        LogLike <- rep(NA, ncol(X0))
        for (i in 1:ncol(X0)) {
          #1 because the first column of X0 is the intercept
          
          X0.test <- as.matrix(X0[, 1:i])
          
          #print("The dim of bk$KW is ")
          #print(dim(bk$KW))
          #print(dim(X0.test))
          #print(dim(CVI))
          
          p3d <-
            GAPIT.EMMAxP3D(
              ys = ys,
              xs = as.matrix(as.data.frame(GD[, 1])),
              K = as.matrix(bk$KW) ,
              Z = Z1,
              X0 = X0.test,
              CVI = CVI,
              CV.Inheritance = CV.Inheritance,
              GI = GI,
              SNP.P3D = SNP.P3D,
              Timmer = Timmer,
              Memory = Memory,
              fullGD = fullGD,
              SNP.permutation = SNP.permutation,
              GP = GP,
              file.path = file.path,
              file.from = file.from,
              file.to = file.to,
              file.total = file.total,
              file.fragment = file.fragment,
              byFile = byFile,
              file.G = file.G,
              file.Ext.G = file.Ext.G,
              file.GD = file.GD,
              file.GM = file.GM,
              file.Ext.GD = file.Ext.GD,
              file.Ext.GM = file.Ext.GM,
              GTindex = GTindex,
              genoFormat = genoFormat,
              optOnly = TRUE,
              SNP.effect = SNP.effect,
              SNP.impute = SNP.impute,
              name.of.trait = name.of.trait,
              Create.indicator = Create.indicator,
              Major.allele.zero = Major.allele.zero
            )
          
          
          
          k.num.param <- 2 + i
          #k is (i-1) because we have the following parameters in the likelihood function:
          #  intercept
          #  (i-1) covariates
          #  sigma_g
          #  delta
          
          #print(paste("The value of round(p3d$REMLs,5) is ", round(p3d$REMLs,5), sep = ""))
          #print(paste("The value of log(GTindex) is ", log(GTindex), sep = ""))
          #print(paste("The value of 0.5*k.num.param*log(GTindex) is ", 0.5*k.num.param*log(nrow(Z1)), sep = ""))
          
          LogLike[i] <- p3d$logLM
          BIC[i] <- p3d$logLM - (0.5 * k.num.param * log(nrow(Z1)))
          
          #print("The value of k.num.param  is: ")
          #print(k.num.param)
          
          #print(paste("The value of nrow(Z1) is ", nrow(Z1), sep = ""))
          
        }
        Optimum.from.BIC <- which(BIC == max(BIC))
        
        print(
          paste(
            "-----------------------The optimal number of PCs/covariates is ",
            (Optimum.from.BIC - 1),
            " -------------------------",
            sep = ""
          )
        )
        
        BIC.Vector <-
          cbind(as.matrix(rep(0:(ncol(
            X0
          ) - 1))), as.matrix(BIC), as.matrix(LogLike))
        
        
        #print(seq(0:ncol(X0)))
        
        #print(BIC.Vector)
        
        colnames(BIC.Vector) <-
          c(
            "Number of PCs/Covariates",
            "BIC (larger is better) - Schwarz 1978",
            "log Likelihood Function Value"
          )
        
        write.table(
          BIC.Vector,
          paste(
            "GAPIT.",
            name.of.trait,
            ".BIC.Model.Selection.Results.csv",
            sep = ""
          ),
          quote = FALSE,
          sep = ",",
          row.names = FALSE,
          col.names = TRUE
        )
        
        #print(BIC.Vector)
        
        X0 <- X0[, 1:(Optimum.from.BIC)]
        
        if (Optimum.from.BIC == 1) {
          X0 <- as.matrix(X0)
        }
        print("The dimension of X0 after model selection is:")
        print(dim(X0))
        
        print("The head of X0 after model selection is")
        print(head(X0))
        
        
      } # where does it start: 522
      
      print("---------------------Sandwich bottom bun-------------------------------")
      print("Compression")
      print(Compression)
      
      #Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Compression")
      #Memory=GAPIT.Memory(Memory=Memory,Infor="Copmression")
      
      if (numSetting == 1)
      {
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "GWAS")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "GWAS")
      }
      
      #Perform GWAS with the optimum setting
      #This section is omited if there is only one setting
      if ((numSetting > 1) |
          (!is.null(sangwich.bottom) & !byPass) | Model.selection) {
        print("Genomic screening...")
        
        optOnly = FALSE  #set default to false and change it to TRUE in these situations:
        if (!hasGenotype)
          optOnly = TRUE
        if (!SNP.test)
          optOnly = TRUE
        
        if (optOnly) {
          colInclude = 1
        } else{
          colInclude = c(1:ncol(GD))
        }
        
        if (numSetting > 1) {
          #Find the best ca,kt and group
          #print(paste(as.numeric(Compression[1,4]))) ###added by Jiabo Wang 2015.7.20
          #print(paste(min(as.numeric(Compression[,4]),rm.na=TRUE)))
          adjust_value = as.numeric(Compression[1, 4]) - min(as.numeric(Compression[, 4]), rm.na =
                                                               TRUE)
          # nocompress_value=as.numeric(Compression[1,4])
          # REML_storage=as.numeric(Compression[,4])
          
          adjust_sq = sqrt(var(as.numeric(Compression[, 4])))
          # threshold=adjust_mean*0.1
          if (which.min(as.numeric(Compression[, 4])) != 1)
            ###added by Jiabo Wang 2015.7.20
          {
            if (which.min(as.numeric(Compression[, 4])) == which.max(as.numeric(Compression[, 5])))
            {
              kt = Compression[which.min(as.numeric(Compression[, 4])), 1]
              ca = Compression[which.min(as.numeric(Compression[, 4])), 2]
              group = Compression[which.min(as.numeric(Compression[, 4])), 3]
              va = Compression[which.min(as.numeric(Compression[, 4])), 5]
              ve = Compression[which.min(as.numeric(Compression[, 4])), 6]
            } else{
              # Compression0=Compression
              cnn = which.min(as.numeric(Compression[, 4]))
              if (cnn - which.min(as.numeric(Compression[-cnn, 4])) < 2)
              {
                kt = Compression[which.min(as.numeric(Compression[, 4])), 1]
                ca = Compression[which.min(as.numeric(Compression[, 4])), 2]
                group = Compression[which.min(as.numeric(Compression[, 4])), 3]
                va = Compression[which.min(as.numeric(Compression[, 4])), 5]
                ve = Compression[which.min(as.numeric(Compression[, 4])), 6]
              } else{
                kt = Compression[1, 1]
                ca = Compression[1, 2]
                group = Compression[1, 3]
                va = Compression[1, 5]
                ve = Compression[1, 6]
                print("The difference of compression is not enough!!")
              }
              
            }
            
            
            
            # Compression=Compression0
            
            print(paste("Compress Optimum: ", ca, kt, group, va, va, ve, sep = " "))
          } else{
            Compression = Compression[order(as.numeric(Compression[, 4]), decreasing = FALSE), ]  #sort on REML
            
            kt = Compression[1, 1]
            ca = Compression[1, 2]
            group = Compression[1, 3]
            print(
              paste(
                "Optimum: ",
                Compression[1, 2],
                Compression[1, 1],
                Compression[1, 3],
                Compression[1, 5],
                Compression[1, 6],
                Compression[1, 4],
                sep = " "
              )
            )
          }
        }#end  if(numSetting>1)
        Compression = Compression[order(as.numeric(Compression[, 4]), decreasing = FALSE), ]
        print(Compression)
        
        
        print("--------------  Sandwich bottom ------------------------")
        
        if (!byPass)
        {
          print("--------------  Sandwich bottom with raw burger------------------------")
          
          if (Model.selection == FALSE) {
            #update KI with the best likelihood
            if (is.null(sangwich.bottom))
              KI = KI.save
            
            cp <-
              GAPIT.Compress(
                KI = KI,
                kinship.cluster = ca,
                kinship.group = kt,
                GN = group,
                Timmer = Timmer,
                Memory = Memory
              )
            Timmer = cp$Timmer
            Memory = cp$Memory
            
            Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "PreP3D 2_cp")
            Memory = GAPIT.Memory(Memory = Memory, Infor = "PreP3D 2_cp")
            
            bk <- GAPIT.Block(Z = Z,
                              GA = cp$GA,
                              KG = cp$KG)
            Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "PreP3D 2_bk")
            Memory = GAPIT.Memory(Memory = Memory, Infor = "PreP3D 2 bk")
            
            zc <- GAPIT.ZmatrixCompress(Z = Z, GAU = bk$GA)
            
            Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "PreP3D 2_zc")
            Memory = GAPIT.Memory(Memory = Memory, Infor = "PreP3D 2 zc")
            
            #Reform KW and Z into EMMA format
            
            z0 = as.matrix(zc$Z[, -1])
            Z1 = matrix(as.numeric(z0),
                        nrow = nrow(z0),
                        ncol = ncol(z0))
          }
          
          print("--------------EMMAxP3D with the optimum setting-----------------------")
          #print(dim(ys))
          #print(dim(as.matrix(as.data.frame(GD[GTindex,colInclude]))))
          p3d <-
            GAPIT.EMMAxP3D(
              ys = ys,
              xs = as.matrix(as.data.frame(GD[GTindex, colInclude]))   ,
              K = as.matrix(bk$KW) ,
              Z = Z1,
              X0 = as.matrix(X0),
              CVI = CVI,
              CV.Inheritance = CV.Inheritance,
              GI = GI,
              SNP.P3D = SNP.P3D,
              Timmer = Timmer,
              Memory = Memory,
              fullGD = fullGD,
              SNP.permutation = SNP.permutation,
              GP = GP,
              file.path = file.path,
              file.from = file.from,
              file.to = file.to,
              file.total = file.total,
              file.fragment = file.fragment,
              byFile = byFile,
              file.G = file.G,
              file.Ext.G = file.Ext.G,
              file.GD = file.GD,
              file.GM = file.GM,
              file.Ext.GD = file.Ext.GD,
              file.Ext.GM = file.Ext.GM,
              GTindex = GTindex,
              genoFormat = genoFormat,
              optOnly = optOnly,
              SNP.effect = SNP.effect,
              SNP.impute = SNP.impute,
              name.of.trait = name.of.trait,
              Create.indicator = Create.indicator,
              Major.allele.zero = Major.allele.zero
            )
          
          Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "GWAS")
          Memory = GAPIT.Memory(Memory = Memory, Infor = "GWAS")
          print("--------------EMMAxP3D with the optimum setting done------------------")
          
        }#end of if(!byPass)
      }#end of if(numSetting>1 & hasGenotype & !SNP.test)
      
      #print("Screening wiht the optimum setting done")
      
      if (byPass)
      {
        print("---------------Sandwich bottom with grilled burger---------------------")
        print("---------------Sandwich bottom: reload bins ---------------------------")
        
        #SUPER: Final screening
        GK = GK.save
        myBread = GAPIT.Bread(
          Y = Y,
          CV = CV,
          Z = Z,
          GK = GK,
          GD = cbind(as.data.frame(GT), as.data.frame(GD)),
          GM = GI,
          method = sangwich.bottom,
          GTindex = GTindex,
          LD = LD,
          file.output = file.output
        )
        
        print("SUPER saving results...")
        
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "GWAS")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "GWAS")
        
        
      }   #end of if(byPass)
      
      print("--------------------Final results presentations------------------------")
      
      
      
      #Plotting optimum group kinship
      if (!byPass)
      {
        if (length(bk$KW) > 1 &
            length(bk$KW) < length(KI) & length(bk$KW) < 1000 & GAPIT3.output) {
          pdf(
            paste("GAPIT.", name.of.trait, ".Kin.Optimum.pdf", sep = ""),
            width = 12,
            height = 12
          )
          par(mar = c(25, 25, 25, 25))
          heatmap.2(
            as.matrix(bk$KW),
            cexRow = .2,
            cexCol = 0.2,
            col = rev(heat.colors(256)),
            scale = "none",
            symkey = FALSE,
            trace = "none"
          )
          dev.off()
        }
      }
      
      
      #Merge GWAS resultss from files to update ps,maf and nobs in p3d
      if (byFile & !fullGD)
      {
        print("Loading GWAS results from file...")
        for (file in file.from:file.to)
        {
          #Initicalization
          frag = 1
          numSNP = file.fragment
          
          while (numSNP == file.fragment) {
            #this is problematic if the read end at the last line
            
            #Initicalization GI to detect reading empty line
            #theGI=NULL
            #theP=NULL
            #theMAF=NULL
            #thenobs=NULL
            
            
            #reload results from files
            print(paste("Current file ", file, "Fragment: ", frag))
            
            theGI <-
              try(read.table(
                paste(
                  "GAPIT.TMP.GI.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                head = TRUE
              )   , silent = TRUE)
            theP <-
              try(read.table(
                paste(
                  "GAPIT.TMP.ps.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                head = FALSE
              )   , silent = TRUE)
            theMAF <-
              try(read.table(
                paste(
                  "GAPIT.TMP.maf.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                head = FALSE
              ), silent = TRUE)
            thenobs <-
              try(read.table(
                paste(
                  "GAPIT.TMP.nobs.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                head = FALSE
              ), silent = TRUE)
            thersquare_base <-
              try(read.table(
                paste(
                  "GAPIT.TMP.rsquare.base.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                head = FALSE
              ), silent = TRUE)
            thersquare <-
              try(read.table(
                paste(
                  "GAPIT.TMP.rsquare.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                head = FALSE
              ), silent = TRUE)
            thedf  <-
              try(read.table(
                paste(
                  "GAPIT.TMP.df.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                head = FALSE
              ), silent = TRUE)
            thetvalue  <-
              try(read.table(
                paste(
                  "GAPIT.TMP.tvalue.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                head = FALSE
              ), silent = TRUE)
            thestderr  <-
              try(read.table(
                paste(
                  "GAPIT.TMP.stderr.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                head = FALSE
              ), silent = TRUE)
            theeffect.est <-
              try(read.table(
                paste(
                  "GAPIT.TMP.effect.est.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                head = FALSE
              ), silent = TRUE)
            
            if (inherits(theGI, "try-error"))  {
              #if(nrow(theGI)<1){
              numSNP = 0
              #print("This fragment is empty.")
            } else{
              #print("Records loaded for this fragment.")
              numSNP = nrow(theGI)
              colnames(theP) = "P"
              colnames(theMAF) = "MAF"
              colnames(thenobs) = "nobs"
              colnames(thersquare_base) = "Base.Model.R.square"
              colnames(thersquare) = "Model.R.square"
              colnames(thedf) = "Model.DF"
              colnames(thetvalue) = "Model.tvalue"
              colnames(thestderr) = "Model.stderr"
              colnames(theeffect.est) = "Effect.Est"
              colnames(theGI) = colnames(GI)
              
              
              
              
              #Merge results
              if (file == file.from & frag == 1) {
                GI = theGI
                #print(dim(GI))
                allP = theP
                #print(head(theP))
                allMAF = theMAF
                allnobs = thenobs
                allrsquare_base = thersquare_base
                allrsquare = thersquare
                alldf = thedf
                alltvalue = thetvalue
                allstderr = thestderr
                alleffect.est = theeffect.est
                
              } else{
                allP = as.data.frame(rbind(as.matrix(allP), as.matrix(theP)))
                allMAF = as.data.frame(rbind(as.matrix(allMAF), as.matrix(theMAF)))
                allnobs = as.data.frame(rbind(as.matrix(allnobs), as.matrix(thenobs)))
                allrsquare_base = as.data.frame(rbind(
                  as.matrix(allrsquare_base),
                  as.matrix(thersquare_base)
                ))
                allrsquare = as.data.frame(rbind(as.matrix(allrsquare), as.matrix(thersquare)))
                alldf = as.data.frame(rbind(as.matrix(alldf), as.matrix(thedf)))
                alltvalue = as.data.frame(rbind(as.matrix(alltvalue), as.matrix(thetvalue)))
                allstderr = as.data.frame(rbind(as.matrix(allstderr), as.matrix(thestderr)))
                alleffect.est = as.data.frame(rbind(
                  as.matrix(alleffect.est),
                  as.matrix(theeffect.est)
                ))
                #print("!!!!!!!!!!!!!!!")
                #print(dim(GI))
                #print(dim(theGI))
                GI = as.data.frame(rbind(as.matrix(GI), as.matrix(theGI)))
              }
              
            }#end of  if(inherits(theGI, "try-error")) (else section)
            
            #setup for next fragment
            frag = frag + 1   #Progress to next fragment
            
          }#end of loop on fragment: while(numSNP==file.fragment)
        }#end of loop on file
        
        #update p3d with components from files
        
        p3d$ps = allP
        p3d$maf = allMAF
        p3d$nobs = allnobs
        p3d$rsquare_base = allrsquare_base
        p3d$rsquare = allrsquare
        p3d$df = alldf
        p3d$tvalue = alltvalue
        p3d$stderr = allstderr
        p3d$effect.est = alleffect.est
        
        #Delete all the GAPIT.TMP files
        theFile = paste("GAPIT.TMP.", name.of.trait, ".*")
        system('cmd /c del "GAPIT.TMP*.*"')
        system('cmd /c del "GAPIT.TMP*.*"')
        print("GWAS results loaded from all files succesfully!")
      } #end of if(byFile)
      
      #--------------------------------------------------------------------------------------------------------------------#
      #Final report
      print("Generating summary")
      GWAS = NULL
      GPS = NULL
      rm(zc)
      gc()
      
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Final")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "Final")
      
      #genomic prediction
      print("Genomic Breeding Values (GBV) ...")
      #print(p3d$BLUP)
      gs = NULL
      if (!byPass)
      {
        if (length(bk$KW) > ncol(X0)) {
          gs <-
            GAPIT.GS(
              KW = bk$KW,
              KO = bk$KO,
              KWO = bk$KWO,
              GAU = bk$GAU,
              UW = cbind(p3d$BLUP, p3d$PEV)
            )
        }
        
        print("Writing GBV and Acc...")
        
        GPS = NULL
        if (length(bk$KW) > ncol(X0))
          GPS = gs$BLUP
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "GPS")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "GPS")
        
        #Make heatmap for distribution of BLUP and PEV
        print("GBV and accuracy distribution...")
        if (length(bk$KW) > ncol(X0) & file.output) {
          GAPIT.GS.Visualization(
            gsBLUP = gs$BLUP,
            BINS = BINS,
            name.of.trait = name.of.trait
          )
        }
        
        #Make a plot Summarzing the Compression Results, if more than one "compression level" has been assessed
        print("Compression portfolios...")
        #print(Compression)
        if (file.output)
          GAPIT.Compression.Visualization(Compression = Compression, name.of.trait = name.of.trait)
        print("Compression Visualization done")
        
        if (length(Compression) < 1) {
          h2.opt = NULL
        } else{
          print(Compression)
          if (length(Compression) < 6)
            Compression = t(as.matrix(Compression[which(Compression[, 4] != "NULL" |
                                                          Compression[, 4] != "NaN"), ]))
          if (length(Compression) == 6)
            Compression = matrix(Compression, 1, 6)
          if (length(Compression) > 6)
            Compression = Compression[which(Compression[, 4] != "NULL" |
                                              Compression[, 4] != "NaN"), ]
          Compression.best = Compression[1, ]
          variance = as.numeric(Compression.best[5:6])
          varp = variance / sum(variance)
          h2.opt = varp[1]
        }
        
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Compression.Visualization")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "Compression.Visualization")
        # print("$$$$$")
        # print(str(p3d))
        
        ps = p3d$ps
        nobs = p3d$nobs
        maf = p3d$maf
        rsquare_base = p3d$rsquare_base
        rsquare = p3d$rsquare
        df = p3d$df
        tvalue = p3d$tvalue
        stderr = p3d$stderr
        effect.est = p3d$effect.est
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Extract p3d results")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "Extract p3d results")
        print("p3d objects transfered")
        
        #where does it start: 936
      } else{
        #byPass
        #print("The head of myBread$GWAS is")
        #print(head(myBread$GWAS))
        GPS = myBread$BLUP
        ps = myBread$GWAS[, 4]
        nobs = myBread$GWAS[, 6]
        #print(dim(GI))
        #print(head())
        Bread_index = match(as.character(myBread$GWAS[, 1]), as.character(GI[, 1]))
        #print(GD[1:5,1:5])
        Bread_X = GD[, Bread_index]
        #print(dim(Bread_X))
        maf = apply(Bread_X, 2, function(one)
          abs(1 - sum(one) / (2 * nrow(Bread_X))))
        maf[maf > 0.5] = 1 - maf[maf > 0.5]
        rsquare_base = rep(NA, length(ps))
        rsquare = rep(NA, length(ps))
        df = rep(NA, length(nobs))
        tvalue = rep(NA, length(nobs))
        stderr = rep(NA, length(nobs))
        effect.est = rep(NA, length(nobs))
        
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Extract bread results")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "Extract bread results")
        
      }
      print("Merge BLUP and BLUE")
      #print(head(ps))
      #Merge BLUP and BLUE
      Pred = NULL
      if ((!byPass) & (!Model.selection)) {
        print("GAPIT before BLUP and BLUE")
        #print(dim(p3d$BLUE))
        BLUE = data.frame(cbind(data.frame(CV.taxa), data.frame(p3d$BLUE)))
        colnames(BLUE) = c("Taxa", "BLUE")
        
        #Initial BLUP as BLUe and add additional columns
        gs.blup = cbind(BLUE, NA, NA, 0, NA)
        
        if (!is.null(gs))
          gs.blup = gs$BLUP

        BB=suppressWarnings (merge(gs.blup, BLUE, by.x = "Taxa", by.y = "Taxa"))
        #BB= (merge(gs.blup, BLUE, by.x = "Taxa", by.y = "Taxa"))
        if (is.null(my_allCV)) {
          my_allX = matrix(1, length(my_taxa), 1)
        } else{
          # my_allX=as.matrix(my_allCV[,-1])
          my_allX = cbind(1, as.matrix(my_allCV[, -1]))
        }
        
        #print(dim(my_allX))
        #print(head(my_allX))
        #print(dim(BB))
        #print(CV.Inheritance)
        if (is.null(CV.Inheritance))
          
        {
          Prediction = BB[, 5] + BB[, 7]
          Pred_Heritable = Prediction
        }
        if (!is.null(CV.Inheritance))
        {
          #inher_CV=my_allX[,1:(1+CV.Inheritance)]
          #beta.Inheritance=p3d$effect.cv[1:(1+CV.Inheritance)]
          #print(beta.Inheritance)
          #if(length(beta)==1)CV=X
          all_BLUE = try(my_allX %*% p3d$effect.cv, silent = T)
          if (inherits(BLUE, "try-error"))
            all_BLUE = NA
          
          
          Pred_Heritable = BB[, 5] + BB[, 7]
          Prediction = BB[, 5] + all_BLUE
        }
        #print("@@@@@@@@@@")
        #print(dim(CVI))
        #print(BB)
        #CV.Inheritance
        #Pred_Heritable=p3d$effect.cv[CV.Inheritance]%*%CVI[CV.Inheritance]+BB[,7]
        Pred = data.frame(cbind(BB, data.frame(Prediction)), data.frame(Pred_Heritable))
        if (noCV)
        {
          if (NOBLUP)
          {
            Pred = NA
          } else{
            BLUE = Pred$BLUE[1]
            prediction = as.matrix(GPS$BLUP) + (BLUE)
            Pred = cbind(GPS, BLUE, prediction)
            colnames(Pred) = c("Taxa",
                               "Group",
                               "RefInf",
                               "ID",
                               "BLUP",
                               "PEV",
                               "BLUE",
                               "Prediction")
          }#end NOBLUP
        }#end noCV
        print("GAPIT after BLUP and BLUE")
      }
      
      #Export BLUP and PEV
      if (!byPass & GAPIT3.output)
      {
        print("Exporting BLUP and Pred")
        #try(write.table(gs$BLUP, paste("GAPIT.", name.of.trait,".BLUP.csv" ,sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE))
        try(write.table(
          Pred,
          paste("GAPIT.", name.of.trait, ".PRED.csv" , sep = ""),
          quote = FALSE,
          sep = ",",
          row.names = FALSE,
          col.names = TRUE
        ))
      }
      
      if (byPass)
      {
        theK.back = NULL
      } else{
        theK.back = cp$KG
      }
      if (byPass)
        Compression[1, 4] = 0 #create a fake value to aloow output of SUPER
      
      #Export GWAS results
      PWI.Filtered = NULL
      if (hasGenotype &
          SNP.test & !is.na(Compression[1, 4]))
        #require not NA REML
      {
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Extract GWAS start")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "Extract GWAS start")
        
        
        #print("Filtering SNPs with MAF..." )
        #index=maf>=SNP.MAF
        
        PWI.Filtered = cbind(GI, ps, maf, nobs, rsquare_base, rsquare, effect.est)#[index,]
        #print(dim(PWI.Filtered))
        colnames(PWI.Filtered) = c(
          "SNP",
          "Chromosome",
          "Position ",
          "P.value",
          "maf",
          "nobs",
          "Rsquare.of.Model.without.SNP",
          "Rsquare.of.Model.with.SNP",
          "effect"
        )
        
        if (!byPass) {
          if (Create.indicator) {
            #Add a counter column for GI
            GI.counter <- cbind(GI, seq(1:nrow(GI)))
            
            #Turn GI and effect.est into data frames
            GI.counter.data.frame <- data.frame(GI.counter)
            colnames(GI.counter.data.frame) <- c("X1", "X2", "X3", "X4")
            
            effect.est.data.frame <- data.frame(effect.est)
            colnames(effect.est.data.frame) <- c("X1", "X2", "X3")
            print(head(GI.counter.data.frame))
            print(head(effect.est.data.frame))
            #Do a merge statement
            GWAS.2 <-
              merge(
                GI.counter.data.frame,
                effect.est.data.frame,
                by.x = "X4",
                by.y = "X1"
              )
            
            #Remove the counter column
            GWAS.2 <- GWAS.2[, -1]
            
            #Add column names
            colnames(GWAS.2) <-
              c("SNP",
                "Chromosome",
                "Position ",
                "Genotype",
                "Allelic Effect Estimate")
            
            
          }
          if (!Create.indicator) {
            GWAS.2 <- PWI.Filtered[, c(1:3, 9)]
            colnames(GWAS.2) <-
              c("SNP",
                "Chromosome",
                "Position ",
                "Allelic Effect Estimate")
          }
        }
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "MAF filtered")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "MAF filtered")
        
        #print("SNPs filtered with MAF")
        
        
        if (!is.null(PWI.Filtered))
        {
          #Run the BH multiple correction procedure of the results
          #Create PWIP, which is a table of SNP Names, Chromosome, bp Position, Raw P-values, FDR Adjusted P-values
          #print("Calculating FDR..." )
          
          PWIP <-
            GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(PWI = PWI.Filtered,
                                                               FDR.Rate = FDR.Rate,
                                                               FDR.Procedure = "BH")
          
          Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Multiple Correction")
          Memory = GAPIT.Memory(Memory = Memory, Infor = "Multiple Correction")
          
          
          #QQ plots
          #print("QQ plot..." )
          if (file.output)
            GAPIT.QQ(
              P.values = PWIP$PWIP[, 4],
              name.of.trait = name.of.trait,
              DPP = DPP
            )
          
          
          Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "QQ plot")
          Memory = GAPIT.Memory(Memory = Memory, Infor = "QQ plot")
          
          
          #Manhattan Plots
          
          
          #print("Manhattan plot (Genomewise)..." )
          #  if(file.output) GAPIT.Manhattan(GI.MP = PWIP$PWIP[,2:4], name.of.trait = name.of.trait, DPP=DPP, plot.type = "Genomewise",cutOff=cutOff)
          #  if(file.output) GAPIT.Manhattan(GI.MP = PWIP$PWIP[,2:4], name.of.trait = name.of.trait, DPP=DPP, plot.type = "Genomewise",cutOff=cutOff,seqQTN=QTN.position)  #QTN does not work with sorted P
          if (file.output)
            GAPIT.Manhattan(
              GI.MP = PWI.Filtered[, 2:4],
              name.of.trait = name.of.trait,
              DPP = DPP,
              plot.type = "Genomewise",
              cutOff = cutOff,
              seqQTN = QTN.position,
              plot.style = plot.style,
              plot.bin = plot.bin,
              chor_taxa = chor_taxa
            )
          
          #print("Manhattan plot (Chromosomewise)..." )
          
          #if(file.output) GAPIT.Manhattan(GI.MP = PWIP$PWIP[,2:4], name.of.trait = name.of.trait, DPP=DPP, plot.type = "Chromosomewise",cutOff=cutOff)
          if (file.output &
              SNP.fraction == 1)
            GAPIT.Manhattan(
              GI.MP = PWI.Filtered[, 2:4],
              GD = GD,
              CG = CG,
              name.of.trait = name.of.trait,
              DPP = DPP,
              plot.type = "Chromosomewise",
              cutOff = cutOff,
              plot.bin = plot.bin,
              chor_taxa = chor_taxa
            )
          
          Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Manhattan plot")
          Memory = GAPIT.Memory(Memory = Memory, Infor = "Manhattan plot")
          
          
          #Association Table
          #print("Association table..." )
          #print(dim(PWIP$PWIP))
          #GAPIT.Table(final.table = PWIP$PWIP, name.of.trait = name.of.trait,SNP.FDR=SNP.FDR)
          #print(head(PWIP$PWIP))
          GWAS = PWIP$PWIP[PWIP$PWIP[, 9] <= SNP.FDR, ]
          #print("Joining tvalue and stderr" )
          
          DTS = cbind(GI, df, tvalue, stderr, effect.est)
          colnames(DTS) = c("SNP",
                            "Chromosome",
                            "Position",
                            "DF",
                            "t Value",
                            "std Error",
                            "effect")
          
          #print("Creating ROC table and plot" )
          if (file.output)
            myROC = GAPIT.ROC(
              t = tvalue,
              se = stderr,
              Vp = var(ys),
              trait = name.of.trait
            )
          #print("ROC table and plot created" )
          
          #MAF plots
          #print("MAF plot..." )
          if (file.output)
            myMAF1 = GAPIT.MAF(
              MAF = GWAS[, 5],
              P = GWAS[, 4],
              E = NULL,
              trait = name.of.trait
            )
          
          
          #print(dim(GWAS))
          
          if (file.output) {
            write.table(
              GWAS,
              paste("GAPIT.", name.of.trait, ".GWAS.Results.csv", sep = ""),
              quote = FALSE,
              sep = ",",
              row.names = FALSE,
              col.names = TRUE
            )
            write.table(
              DTS,
              paste(
                "GAPIT.",
                name.of.trait,
                ".Df.tValue.StdErr.csv",
                sep = ""
              ),
              quote = FALSE,
              sep = ",",
              row.names = FALSE,
              col.names = TRUE
            )
            if (!byPass)
              write.table(
                GWAS.2,
                paste(
                  "GAPIT.",
                  name.of.trait,
                  ".Allelic_Effect_Estimates.csv",
                  sep = ""
                ),
                quote = FALSE,
                sep = ",",
                row.names = FALSE,
                col.names = TRUE
              )
          }
          
          
          
        } #end of if(!is.null(PWI.Filtered))
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Extract GWAS end")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "Extract GWAS end")
        
        
      } #end of if(hasGenotype )
      
      #Log
      if (GAPIT3.output)
        log = GAPIT.Log(
          Y = Y,
          KI = KI,
          Z = Z,
          CV = CV,
          SNP.P3D = SNP.P3D,
          group.from = group.from ,
          group.to = group.to ,
          group.by = group.by ,
          kinship.cluster = kinship.cluster,
          kinship.group = kinship.group,
          ngrid = ngrid ,
          llin = llin ,
          ulim = ulim ,
          esp = esp ,
          name.of.trait = name.of.trait
        )
      #Memory usage
      #GAPIT.Memory.Object(name.of.trait=name.of.trait)
      
      #Timming
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Report")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "Report")
      if (file.output) {
        file = paste("GAPIT.", name.of.trait, ".Timming.csv" , sep = "")
        write.table(
          Timmer,
          file,
          quote = FALSE,
          sep = ",",
          row.names = FALSE,
          col.names = TRUE
        )
        
        file = paste("GAPIT.", name.of.trait, ".Memory.Stage.csv" , sep = "")
        write.table(
          Memory,
          file,
          quote = FALSE,
          sep = ",",
          row.names = FALSE,
          col.names = TRUE
        )
      }
      print(paste(name.of.trait, "has been analyzed successfully!"))
      print(paste("The results are saved in the directory of ", getwd()))
      
      
      
      #print("==========================================================================================")
      TV <- list()
      TV$ps = ps
      TV$nobs = nobs
      TV$maf = maf
      TV$rsquare_base = rsquare_base
      TV$rsquare = rsquare
      TV$df = df
      TV$tvalue = tvalue
      TV$stderr = stderr
      TV$effect.est = effect.est
      #print("!!!!!!!!!!!!!")
      #print(head(effect.est))
      #print(head(DTS[,7]))
      #print(ys)
      if (byPass | Model.selection)
        Pred <- NA
      print("before ending GAPIT.Main")
      #print(dim(Compression))
      return (
        list(
          Timmer = Timmer,
          Compression = Compression,
          kinship.optimum = theK.back,
          kinship = KI,
          PC = PC,
          GWAS = PWI.Filtered,
          GPS = GPS,
          Pred = Pred,
          REMLs = Compression[count, 4],
          Timmer = Timmer,
          Memory = Memory,
          SUPER_GD = SUPER_GD,
          P = ps,
          effect.snp = DTS[, 7],
          effect.cv = p3d$effect.cv,
          h2 = h2.opt,
          TV = TV
        )
      )
    } #end if non-SUPER.GS situation, this is a long if statement, structure needs improvement
  }#The function GAPIT.Main ends here


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.RemoveDuplicate` <-
  function(Y) {
    #Object: NA
    #Output: NA
    #Authors: Zhiwu Zhang
    # Last update: Augus 30, 2011
    ##############################################################################################
    return (Y[match(unique(Y[, 1]), Y[, 1], nomatch = 0),])
  }

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.QC` <-
  function(Y = NULL,
           KI = NULL,
           GT = NULL,
           CV = NULL,
           Z = NULL,
           GK = NULL) {
    #Object: to do data quality control
    #Output: Y, KI, GD, CV, Z, flag
    #Authors: Zhiwu Zhang and Alex Lipka
    # Last update: April 14, 2011
    ##############################################################################################
    #Remove duplicates
    print("Removing duplicates...")
    #print(dim(CV))
    Y = GAPIT.RemoveDuplicate(Y)
    CV = GAPIT.RemoveDuplicate(CV)
    GK = GAPIT.RemoveDuplicate(GK)
    if (!is.null(Z))
      Z = GAPIT.RemoveDuplicate(Z)
    
    #Remove missing phenotype
    print("Removing NaN...")
    Y = Y[which(Y[, 2] != "NaN"), ]
    
    # Remove duplicates for GT
    # GT row wise, Z column wise, and KI both direction.
    print("Remove duplicates for GT...")
    #print(dim(GT))
    if (!is.null(GT))
    {
      if (is.null(dim(GT)))
        taxa.kept = unique(GT)
      if (!is.null(dim(GT)))
        taxa.kept = unique(GT[, 1])
      
    } else{
      taxa.kept = unique(Y[, 1])
    }
    
    # Remove duplicates for KI
    print("Remove duplicates for KI...")
    # improve speed: remove t() and use cbind
    if (!is.null(KI))
    {
      taxa.all = KI[, 1]
      taxa.uniqe = unique(taxa.all)
      position = match(taxa.uniqe, taxa.all, nomatch = 0)
      position.addition = cbind(1, t(1 + position))
      KI = KI[position, position.addition]
    }
    
    #Sort KI
    if (!is.null(KI))
    {
      taxa.all = KI[, 1]
      position = order(taxa.all)
      position.addition = cbind(1, t(1 + position))
      KI = KI[position, position.addition]
    }
    
    # Remove duplicates for Z rowwise
    print("Remove duplicates for Z (column wise)...")
    if (!is.null(Z))
    {
      taxa.all = as.matrix(Z[1, ])
      taxa.uniqe = intersect(taxa.all, taxa.all)
      position = match(taxa.uniqe, taxa.all, nomatch = 0)
      Z = Z[, position]
    }
    
    
    #Remove the columns of Z if they are not in KI/GT. KI/GT are allowed to have individuals not in Z
    print("Maching Z with Kinship colwise...")
    if (!is.null(KI))
    {
      taxa.all = KI[, 1]
      taxa.kinship = unique(taxa.all)
    }
    
    if (!is.null(Z) & !is.null(KI))
    {
      #get common taxe between KI and Z
      taxa.Z = as.matrix(Z[1, ])
      #taxa.Z=colnames(Z) #This does not work for names starting with numerical or "-"   \
      if (is.null(KI)) {
        taxa.Z_K_common = taxa.Z
      } else{
        taxa.Z_K_common = intersect(taxa.kinship, taxa.Z)
      }
      Z <- cbind(Z[, 1], Z[, match(taxa.Z_K_common, taxa.Z, nomatch = 0)])
      
      #Remove the rows of Z if all the ellements sum to 0
      #@@@ improve speed: too many Zs
      print("Maching Z without origin...")
      Z1 = Z[-1, -1]
      Z2 = data.frame(Z1)
      Z3 = as.matrix(Z2)
      Z4 = as.numeric(Z3) #one dimemtion
      Z5 = matrix(data = Z4,
                  nrow = nrow(Z1),
                  ncol = ncol(Z1))
      RS = rowSums(Z5) > 0
      #The above process could be simplified!
      Z <- Z[c(TRUE, RS), ]
      
      #make individuals the same in Z, Y, GT and CV
      print("Maching GT and CV...")
      if (length(Z) <= 1)
        stop("GAPIT says: there is no place to match IDs!")
    }# end of  if(!is.null(Z) & !is.null(K))
    
    # get intersect of all the data
    taxa = intersect(Y[, 1], Y[, 1])
    if (!is.null(Z))
      taxa = intersect(Z[-1, 1], taxa)
    if (!is.null(GT))
      taxa = intersect(taxa, taxa.kept)
    if (!is.null(CV))
      taxa = intersect(taxa, CV[, 1])
    if (!is.null(GK))
      taxa = intersect(taxa, GK[, 1])
    if (length(taxa) <= 1)
      stop("GAPIT says: There is no individual ID matched to covariate. Please check!")
    
    
    if (!is.null(Z))
    {
      #Remove taxa in Z that are not in others, columnwise
      t = c(TRUE, Z[-1, 1] %in% taxa)
      if (length(t) <= 2)
        stop("GAPIT says: There is no individual ID matched among data. Please check!")
      Z <- Z[t, ]
      
      #Remove the columns of Z if all the ellements sum to 0
      print("QC final process...")
      #@@@ improve speed: too many Zs
      Z1 = Z[-1, -1]
      Z2 = data.frame(Z1)
      Z3 = as.matrix(Z2)
      Z4 = as.numeric(Z3) #one dimemtion
      Z5 = matrix(data = Z4,
                  nrow = nrow(Z1),
                  ncol = ncol(Z1))
      CS = colSums(Z5) > 0
      #The above process could be simplified!
      Z <- Z[, c(TRUE, CS)]
    }
    
    #Filtering with comman taxa
    Y <- Y[Y[, 1] %in% taxa, ]
    if (!is.null(CV))
      CV = CV[CV[, 1] %in% taxa, ]
    if (!is.null(GK))
      GK = GK[GK[, 1] %in% taxa, ]
    if (!is.null(GT))
      taxa.kept = data.frame(taxa.kept[taxa.kept %in% taxa])
    #Y <- Y[Y[,1]%in%taxa.kept,]
    
    #To sort Y, GT, CV and Z
    Y = Y[order(Y[, 1]), ]
    CV = CV[order(CV[, 1]), ]
    if (!is.null(GK))
      GK = GK[order(GK[, 1]), ]
    if (!is.null(Z))
      Z = Z[c(1, 1 + order(Z[-1, 1])), ]
    
    #get position of taxa.kept in GT
    #position=match(taxa.kept[,1], GT[,1],nomatch = 0)
    if (is.null(dim(GT)))
      position = match(taxa.kept, GT, nomatch = 0)
    if (!is.null(dim(GT)))
      position = match(taxa.kept[, 1], GT[, 1], nomatch = 0)
    
    
    if (is.null(dim(taxa.kept)))
      order.taxa.kept = order(taxa.kept)
    if (!is.null(dim(taxa.kept)))
      order.taxa.kept = order(taxa.kept[, 1])
    
    GTindex = position[order.taxa.kept]
    flag = nrow(Y) == nrow(Z) - 1 & nrow(Y) == nrow(GT) &
      nrow(Y) == nrow(CV)
    
    print("GAPIT.QC accomplished successfully!")
    
    #print(dim(Y))
    #print(dim(CV))
    #print(dim(KI))
    return(list(
      Y = Y,
      KI = KI,
      GT = GT,
      CV = CV,
      Z = Z,
      GK = GK,
      GTindex = GTindex,
      flag = flag
    ))
  }#The function GAPIT.QC ends here

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.Compress` <-
  function(KI,
           kinship.cluster = "average",
           kinship.group = "Mean",
           GN = nrow(KI),
           Timmer,
           Memory) {
    #Object: To cluster individuals into groups based on kinship
    #Output: GA, KG
    #Authors: Alex Lipka and Zhiwu Zhang
    # Last update: April 14, 2011
    ##############################################################################################
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "CP start")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "cp start")
    
    # Extract the line names
    line.names <- KI[, 1]
    
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Does this change memory0")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "Does this change memory0")
    
    # Remove the first column of the kinship matrix, which is the line names
    KI <- KI[, -1]
    
    # Convert kinship to distance
    #distance.matrix <- 2 - KI
    
    
    #distance.matrix.as.dist <- as.dist(distance.matrix)
    #distance.matrix.as.dist <- as.dist(2 - KI)
    
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "CP distance")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "cp distance")
    
    #print(paste("The value of kinship.cluster is ", kinship.cluster, sep = ""))
    
    
    
    # hclust() will perform the hiearchical cluster analysis
    #cluster.distance.matrix <- hclust(distance.matrix.as.dist, method = kinship.cluster)
    #cluster.distance.matrix <- hclust(as.dist(2 - KI), method = kinship.cluster)
    distance.matrix = dist(KI, upper = TRUE) #Jiabo Wang modified ,the dist is right function for cluster
    cluster.distance.matrix = hclust(distance.matrix, method = kinship.cluster)
    #cutree(out_hclust,k=3)
    
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "CP cluster")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "cp cluster")
    
    # Cutree will assign lines into k clusters
    group.membership <- cutree(cluster.distance.matrix, k = GN)
    compress_z = table(group.membership, paste(line.names))  #build compress z with group.membership
    
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "CP cutree")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "cp cutree")
    
    #calculate group kinship
    if (kinship.group == "Mean") {
      #This matrix ooperation is much faster than tapply function for  "Mean"
      x = as.factor(group.membership)
      #b = model.matrix(~x-1)
      n = max(as.numeric(as.vector(x)))
      b = diag(n)[x, ]
      
      KG = t(b) %*% as.matrix(KI) %*% b
      CT = t(b) %*% (0 * as.matrix(KI) + 1) %*% b
      KG = as.matrix(KG / CT)
      rownames(KG) = c(1:nrow(KG))
      colnames(KG) = c(1:ncol(KG))
      
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "CP calculation original")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "cp calculation original")
      
      
      
    } else{
      gm = as.factor(group.membership)
      kv = as.numeric(as.matrix(KI))
      kvr = rep(gm, ncol(KI))
      kvc = as.numeric(t(matrix(kvr, nrow(KI), ncol(KI))))
      
      kInCol = t(rbind(kv, kvr, kvc))
      
      rm(gm)
      rm(kv)
      rm(kvr)
      rm(kvc)
      rm(KI)
      gc()
      
      
      
      #This part does not work yet
      #if(kinship.group == "Mean")
      #    KG<- tapply(kInCol[,1], list(kInCol[,2], kInCol[,3]), mean)
      if (kinship.group == "Max")
        KG <- tapply(kInCol[, 1], list(kInCol[, 2], kInCol[, 3]), max)
      if (kinship.group == "Min")
        KG <- tapply(kInCol[, 1], list(kInCol[, 2], kInCol[, 3]), min)
      if (kinship.group == "Median")
        KG <- tapply(kInCol[, 1], list(kInCol[, 2], kInCol[, 3]), median)
    } #this is end of brancing "Mean" and the rest
    
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "CP calculation")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "cp calculation")
    
    # add line names
    #GA <- data.frame(group.membership)
    GA <-
      data.frame(cbind(as.character(line.names), as.numeric(group.membership)))
    
    #Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="CP Final")
    #Memory=GAPIT.Memory(Memory=Memory,Infor="CP Final")
    
    #write.table(KG, paste("KG_from_", kinship.group, "_Method.txt"), quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)
    
    #print("GAPIT.Compress accomplished successfully!")
    return(list(
      GA = GA,
      KG = KG,
      Timmer = Timmer,
      Memory = Memory
    ))
  }#The function GAPIT.Compress ends here

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.Block` <-
  function(Z, GA, KG) {
    #Object: To split a group kinship into two blocks containing individuals with and without phenotype
    #Output: GAU,KW,KO,KWO
    #Authors: Zhiwu Zhang and Alex Lipka
    # Last update: April 14, 2011
    ##############################################################################################
    # To separate group kiship into two blocks: with and without phenotype.
    # A group goes to with phenotype as loog as it has one phenotyped individual.
    
    #find position in group assignment (GA) for the individual associate with phenotype (specified by Z)
    #taxa=unique(intersect(as.matrix(Z[1,-1]),GA[,1]))
    
    taxa.Z = as.matrix(Z[1, -1])
    taxa.GA = as.matrix(GA[, 1])
    position = taxa.GA %in% taxa.Z
    
    #Initial block as 2
    GAU = cbind(GA, 2)
    
    #Assign block as 1 if the individual has phenotype
    GAU[position, 3] = 1
    
    #Modify the non-phenotyped individuals if they in a group with phenotyped individuals
    #To find the groups with phenotyped individuals
    #update block assignment for all these groups
    #get list of group that should be block 1
    
    #grp.12=as.matrix(unique(GAU[,2]))
    #grp.1=as.matrix(unique(GAU[which(GAU[,3]==1),2]))
    #grp.2= as.matrix(setdiff(grp.12,grp.1))
    
    grp.12 = as.matrix(as.vector(unique(GAU[, 2]))) #unique group
    grp.1 = as.matrix(as.vector(unique(GAU[which(GAU[, 3] == 1), 2]))) #unique phenotyped group
    grp.2 = as.matrix(as.vector(setdiff(grp.12, grp.1))) #unique unphenotyped group
    
    numWithout = length(grp.2)
    
    order.1 = 1:length(grp.1)
    order.2 = 1:length(grp.2)
    if (numWithout > 0)
      grpblock = as.matrix(rbind(cbind(grp.1, 1, order.1), cbind(grp.2,   2,    order.2)))
    if (numWithout == 0)
      grpblock = as.matrix(cbind(grp.1, 1, order.1),)
    
    order.block = order(as.matrix(GAU[, 3]))
    colnames(grpblock) = c("grp", "block", "ID")
    
    #Indicators: 1-Phenotype, 1.5- unphenotyped but in a group with other phenotyped, 2-rest  (Zhiwu, Dec 7,2012)
    #GAU0 <- merge(GAU[order.block,-3], grpblock, by.x = "X2", by.y = "grp")
    #GAU=GAU0[,c(2,1,3,4)]
    #print(head(GAU))
    GAU1 <-
      merge(GAU[order.block, ], grpblock, by.x = "X2", by.y = "grp")
    #print(GAU1)
    GAU1[, 4] = (as.numeric(GAU1[, 3]) + as.numeric(GAU1[, 4])) / 2
    #print(GAU1)
    
    GAU = GAU1[, c(2, 1, 4, 5)]
    KW = KG[grp.1, grp.1]
    KO = KG[grp.2, grp.2]
    KWO = KG[grp.1, grp.2]
    
    #write.table(GAU, "GAU.txt", quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)
    
    #print("GAPIT.Block accomplished successfully!")
    
    return(list(
      GAU = GAU,
      KW = KW,
      KO = KO,
      KWO = KWO
    ))
  }#The function GAPIT.Block ends here

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.ZmatrixCompress` <-
  function(Z, GAU) {
    #Object: To assign the fraction of a individual belonging to a group
    #Output: Z
    #Authors: Zhiwu Zhang
    # Last update: April 14, 2011
    ##############################################################################################
    #Extraction of GAU coresponding to Z, sort GAU rowwise to mach columns of Z, and make design matrix
    #print("GAPIT.ZmatrixCompress")
    #print(dim(Z))
    #print(dim(GAU))
    
    effect.Z = as.matrix(Z[1, -1])
    effect.GAU = as.matrix(GAU[, 1])
    taxa = as.data.frame(Z[-1, 1])
    
    GAU0 = GAU[effect.GAU %in% effect.Z, ]
    order.GAU = order(GAU0[, 1])
    GAU1 <- GAU0[order.GAU, ]
    #id.1=GAU1[which(GAU1[,3]==1),4]
    id.1 = GAU1[which(GAU1[, 3] < 2), 4]
    n = max(as.numeric(as.vector(id.1)))
    x = as.numeric(as.matrix(GAU1[, 4]))
    DS = diag(n)[x, ]
    
    #sort Z column wise
    order.Z = order(effect.Z)
    Z = Z[-1, -1]
    Z <- Z[, order.Z]
    
    #Z matrix from individual to group
    #Z1.numeric <- as.numeric(as.matrix(Z))
    Z <-
      matrix(as.numeric(as.matrix(Z)),
             nrow = nrow(Z),
             ncol = ncol(Z))
    Z = Z %*% DS
    
    #Z3=data.frame(cbind(as.character(Z[-1,1]),Z2))
    Z = data.frame(cbind(taxa, Z))
    
    #Z=Z3[order(Z3[,1]),]
    
    # Z=Z[order(as.matrix(taxa)),]
    
    
    #print("GAPIT.ZmatrixCompress accomplished successfully!")
    return(list(Z = Z))
  }#The function GAPIT.ZmatrixCompress ends here

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.EMMAxP3D` <-
  function(ys,
           xs,
           K = NULL,
           Z = NULL,
           X0 = NULL,
           CVI = NULL,
           CV.Inheritance = NULL,
           GI = NULL,
           GP = NULL,
           file.path = NULL,
           file.from = NULL,
           file.to = NULL,
           file.total = 1,
           genoFormat = "Hapmap",
           file.fragment = NULL,
           byFile = FALSE,
           fullGD = TRUE,
           SNP.fraction = 1,
           file.G = NULL,
           file.Ext.G = NULL,
           GTindex = NULL,
           file.GD = NULL,
           file.GM = NULL,
           file.Ext.GD = NULL,
           file.Ext.GM = NULL,
           SNP.P3D = TRUE,
           Timmer,
           Memory,
           optOnly = TRUE,
           SNP.effect = "Add",
           SNP.impute = "Middle",
           SNP.permutation = FALSE,
           ngrids = 100,
           llim = -10,
           ulim = 10,
           esp = 1e-10,
           name.of.trait = NULL,
           Create.indicator = FALSE,
           Major.allele.zero = FALSE) {
    #Object: To esimate variance component by using EMMA algorithm and perform GWAS with P3D/EMMAx
    #Output: ps, REMLs, stats, dfs, vgs, ves, BLUP,  BLUP_Plus_Mean, PEV
    #Authors: Feng Tian, Alex Lipka and Zhiwu Zhang
    # Last update: April 6, 2016
    # Library used: EMMA (Kang et al, Genetics, Vol. 178, 1709-1723, March 2008)
    # Note: This function was modified from the function of emma.REML.t from the library
    ##############################################################################################
    #print("EMMAxP3D started...")
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "P3D Start")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "P3D Start")
    
    
    #When numeric genotypes are selected, impute the missing SNPs with the allele indicated by the "SNP.impute" value
    if (!optOnly) {
      if (SNP.impute == "Major")
        xs[which(is.na(xs))] = 2
      if (SNP.impute == "Minor")
        xs[which(is.na(xs))] = 0
      if (SNP.impute == "Middle")
        xs[which(is.na(xs))] = 1
    }
    
    
    #--------------------------------------------------------------------------------------------------------------------<
    #Change data to matrix format if they are not
    if (is.null(dim(ys)) ||
        ncol(ys) == 1)
      ys <- matrix(ys, 1, length(ys))
    if (is.null(X0))
      X0 <- matrix(1, ncol(ys), 1)
    
    #handler of special Z and K
    if (!is.null(Z)) {
      if (ncol(Z) == nrow(Z))
        Z = NULL
    }
    if (!is.null(K)) {
      if (length(K) <= 1)
        K = NULL
    }
    
    #Extract dimension information
    g <- nrow(ys) #number of traits
    n <- ncol(ys) #number of observation
    
    q0 <- ncol(X0)#number of fixed effects
    q1 <- q0 + 1  #Nuber of fixed effect including SNP
    
    nr = n
    if (!is.null(K))
      tv = ncol(K)
    
    #decomposation without fixed effect
    #print("Caling emma.eigen.L...")
    if (!is.null(K))
      eig.L <-
      emma.eigen.L(Z, K) #this function handle both NULL Z and non-NULL Z matrix
    
    #eig.L$values[eig.L$values<0]=0
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "eig.L")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "eig.L")
    
    #decomposation with fixed effect (SNP not included)
    #print("Calling emma.eigen.R.w.Z...")
    X <-  X0 #covariate variables such as population structure
    if (!is.null(Z) &
        !is.null(K))
      eig.R <-
      try(emma.eigen.R.w.Z(Z, K, X), silent = TRUE)
    #This will be used to get REstricted ML (REML)
    if (is.null(Z)  &
        !is.null(K))
      eig.R <-
      try(emma.eigen.R.wo.Z(K, X), silent = TRUE)
    #This will be used to get REstricted ML (REML)
    
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "eig.R")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "eig.R")
    
    #eig.R$values[eig.R$values<0]=0
    #print(labels(eig.R))
    #print(length(eig.R$values))
    #print(dim(eig.R$vectors))
    #print("emma.eigen.R.w.Z called!!!")
    #Handler of error in emma
    #print("!!!!!!")
    if (!is.null(K)) {
      if (inherits(eig.R, "try-error"))
        return(
          list(
            ps = NULL,
            REMLs = NA,
            stats = NULL,
            effect.est = NULL,
            dfs = NULL,
            maf = NULL,
            nobs = NULL,
            Timmer = Timmer,
            Memory = Memory,
            vgs = NA,
            ves = NA,
            BLUP = NULL,
            BLUP_Plus_Mean = NULL,
            PEV = NULL,
            BLUE = NULL
          )
        )
      
      #print("@@@@@")
    }
    #-------------------------------------------------------------------------------------------------------------------->
    #print("Looping through traits...")
    #Loop on Traits
    for (j in 1:g)
    {
      if (optOnly) {
        #REMLE <- GAPIT.emma.REMLE(ys[j,], X, K, Z, ngrids, llim, ulim, esp, eig.R)
        #vgs <- REMLE$vg
        #ves <- REMLE$ve
        #REMLs <- REMLE$REML
        #REMLE_delta=REMLE$delta
        
        if (!is.null(K)) {
          REMLE <-
            GAPIT.emma.REMLE(ys[j, ], X, K, Z, ngrids, llim, ulim, esp, eig.R)
          
          Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "REML")
          Memory = GAPIT.Memory(Memory = Memory, Infor = "REML")
          
          rm(eig.R)
          gc()
          Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "eig.R removed")
          Memory = GAPIT.Memory(Memory = Memory, Infor = "eig.R removed")
          
          vgs <- REMLE$vg
          ves <- REMLE$ve
          REMLs <- REMLE$REML
          REMLE_delta = REMLE$delta
          
          rm(REMLE)
          gc()
        }
        
        
        vids <- !is.na(ys[j, ])
        yv <- ys[j, vids]
        
        if (!is.null(Z) &
            !is.null(K))
          U <-
          eig.L$vectors * matrix(c(sqrt(1 / (
            eig.L$values + REMLE_delta
          )), rep(sqrt(
            1 / REMLE_delta
          ), nr - tv)), nr, ((nr - tv) + length(eig.L$values)), byrow = TRUE)
        if (is.null(Z) &
            !is.null(K))
          U <-
          eig.L$vectors * matrix(sqrt(1 / (eig.L$values + REMLE_delta)), nr, length(eig.L$values), byrow =
                                   TRUE)
        
        if (!is.null(Z) &
            !is.null(K))
          eig.full.plus.delta <-
          as.matrix(c((eig.L$values + REMLE_delta), rep(REMLE_delta, (nr - tv))))
        if (is.null(Z) &
            !is.null(K))
          eig.full.plus.delta <- as.matrix((eig.L$values + REMLE_delta))
        
        if (!is.null(K)) {
          if (length(which(eig.L$values < 0)) > 0) {
            #print("---------------------------------------------------The group kinship matrix at this compression level is not positive semidefinite. Please select another compression level.---------------------------------------------------")
            #return(list(ps = NULL, REMLs = 999999, stats = NULL, effect.est = NULL, dfs = NULL,maf=NULL,nobs = NULL,Timmer=Timmer,Memory=Memory,
            #vgs = 1.000, ves = 1.000, BLUP = NULL, BLUP_Plus_Mean = NULL,
            #PEV = NULL, BLUE=NULL))
          }
        }
        
        
        #Calculate the log likelihood function for the intercept only model
        
        X.int <- matrix(1, nrow(as.matrix(yv)), ncol(as.matrix(yv)))
        iX.intX.int <- solve(crossprod(X.int, X.int))
        iX.intY <- crossprod(X.int, as.matrix(as.matrix(yv)))
        beta.int <-
          crossprod(iX.intX.int, iX.intY)  #Note: we can use crossprod here becase iXX is symmetric
        X.int.beta.int <- X.int %*% beta.int
        
        
        logL0 <- 0.5 * ((-length(yv)) * log(((2 * pi) / length(yv))
                                            * crossprod((yv - X.int.beta.int), (yv - X.int.beta.int)
                                            ))
                        - length(yv))
        
        #print(paste("The value of logL0 inside of the optonly template is is",logL0, sep = ""))
        
        
        
        
        #print(paste("The value of nrow(as.matrix(ys[!is.na(ys)])) is ",nrow(as.matrix(ys[!is.na(ys)])), sep = ""))
        
        
        
        if (!is.null(K)) {
          yt <- yt <- crossprod(U, yv)
          X0t <- crossprod(U, X0)
          
          X0X0 <- crossprod(X0t, X0t)
          X0Y <- crossprod(X0t, yt)
          XY <- X0Y
          
          iX0X0 <- try(solve(X0X0), silent = TRUE)
          if (inherits(iX0X0, "try-error")) {
            iX0X0 <- ginv(X0X0)
            print(
              "At least two of your covariates are linearly dependent. Please reconsider the covariates you are using for GWAS and GPS"
            )
          }
          iXX <- iX0X0
        }
        
        if (is.null(K)) {
          iXX <- try(solve(crossprod(X, X)), silent = TRUE)
          if (inherits(iXX, "try-error"))
            iXX <- ginv(crossprod(X, X))
          XY = crossprod(X, yv)
        }
        beta <-
          crossprod(iXX, XY) #Note: we can use crossprod here because iXX is symmetric
        X.beta <- X %*% beta
        
        beta.cv = beta
        BLUE = X.beta
        
        if (!is.null(K)) {
          U.times.yv.minus.X.beta <- crossprod(U, (yv - X.beta))
          logLM <-
            0.5 * (-length(yv) * log(((2 * pi) / length(yv)) * crossprod(U.times.yv.minus.X.beta, U.times.yv.minus.X.beta)
            )
            - sum(log(eig.full.plus.delta)) - length(yv))
        }
        
        
        if (is.null(K)) {
          U.times.yv.minus.X.beta <- yv - X.beta
          logLM <-
            0.5 * (-length(yv) * log(((2 * pi) / length(yv)) * crossprod(U.times.yv.minus.X.beta, U.times.yv.minus.X.beta)
            ) - length(yv))
        }
        
      }#End if(optOnly)
      
      
      
      #--------------------------------------------------------------------------------------------------------------------<
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Trait")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "Trait")
      
      if (!is.null(K)) {
        REMLE <-
          GAPIT.emma.REMLE(ys[j, ], X, K, Z, ngrids, llim, ulim, esp, eig.R)
        
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "REML")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "REML")
        
        rm(eig.R)
        gc()
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "eig.R removed")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "eig.R removed")
        
        vgs <- REMLE$vg
        ves <- REMLE$ve
        REMLs <- REMLE$REML
        REMLE_delta = REMLE$delta
        
        rm(REMLE)
        gc()
      }
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "REMLE removed")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "REMLE removed")
      
      if (!is.null(Z) &
          !is.null(K))
        U <-
        eig.L$vectors * matrix(c(sqrt(1 / (
          eig.L$values + REMLE_delta
        )), rep(sqrt(1 / REMLE_delta), nr - tv)), nr, ((nr - tv) + length(eig.L$values)), byrow =
          TRUE)
      if (is.null(Z) &
          !is.null(K))
        U <-
        eig.L$vectors * matrix(sqrt(1 / (eig.L$values + REMLE_delta)), nr, length(eig.L$values), byrow =
                                 TRUE)
      
      if (!is.null(Z) &
          !is.null(K))
        eig.full.plus.delta <-
        as.matrix(c((eig.L$values + REMLE_delta), rep(REMLE_delta, (nr - tv))))
      if (is.null(Z) &
          !is.null(K))
        eig.full.plus.delta <- as.matrix((eig.L$values + REMLE_delta))
      
      
      
      if (!is.null(K)) {
        if (length(which(eig.L$values < 0)) > 0) {
          #print("---------------------------------------------------The group kinship matrix at this compression level is not positive semidefinite. Please select another compression level.---------------------------------------------------")
          #return(list(ps = NULL, REMLs = 999999, stats = NULL, effect.est = NULL, dfs = NULL,maf=NULL,nobs = NULL,Timmer=Timmer,Memory=Memory,
          #vgs = 1.000, ves = 1.000, BLUP = NULL, BLUP_Plus_Mean = NULL,
          #PEV = NULL, BLUE=NULL))
        }
      }
      
      
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "U Matrix")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "U Matrix")
      
      if (SNP.P3D == TRUE)
        rm(eig.L)
      gc()
      
      Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "eig.L removed")
      Memory = GAPIT.Memory(Memory = Memory, Infor = "eig.L removed")
      
      #-------------------------------------------------------------------------------------------------------------------->
      
      #The cases that go though multiple file once
      file.stop = file.to
      if (optOnly)
        file.stop = file.from
      if (fullGD)
        file.stop = file.from
      if (!fullGD & !optOnly) {
        print("Screening SNPs from file...")
      }
      
      #Add loop for genotype data files
      for (file in file.from:file.stop)
      {
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "New Genotype file")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "New Genotype file")
        
        
        
        frag = 1
        numSNP = file.fragment
        myFRG = NULL
        while (numSNP == file.fragment) {
          #this is problematic if the read end at the last line
          
          
          
          #initial previous SNP storage
          x.prev <- vector(length = 0)
          
          #force to skip the while loop if optOnly
          if (optOnly)
            numSNP = 0
          
          #Determine the case of first file and first fragment: skip read file
          if (file == file.from & frag == 1 & SNP.fraction < 1) {
            firstFileFirstFrag = TRUE
          } else{
            firstFileFirstFrag = FALSE
          }
          
          #In case of xs is not full GD, replace xs from file
          if (!fullGD & !optOnly & !firstFileFirstFrag)
          {
            Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Clean myFRG")
            Memory = GAPIT.Memory(Memory = Memory, Infor = "Clean myFRG")
            
            #update xs for each file
            rm(xs)
            rm(myFRG)
            gc()
            print(paste("Current file: ", file, " , Fragment: ", frag, sep = ""))
            
            Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Read file fragment")
            Memory = GAPIT.Memory(Memory = Memory, Infor = "Read file fragment")
            
            myFRG = GAPIT.Fragment(
              file.path = file.path,
              file.total = file.total,
              file.G = file.G,
              file.Ext.G = file.Ext.G,
              seed = seed,
              SNP.fraction = SNP.fraction,
              SNP.effect = SNP.effect,
              SNP.impute = SNP.impute,
              genoFormat = genoFormat,
              file.GD = file.GD,
              file.Ext.GD = file.Ext.GD,
              file.GM = file.GM,
              file.Ext.GM = file.Ext.GM,
              file.fragment = file.fragment,
              file = file,
              frag = frag,
              Create.indicator = Create.indicator,
              Major.allele.zero = Major.allele.zero
            )
            
            Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Genotype file converted")
            Memory = GAPIT.Memory(Memory = Memory, Infor = "Genotype file converted")
            
            #print("-----------------------------------------------------------------")
            
            if (is.null(myFRG$GD)) {
              xs = NULL
            } else{
              xs = myFRG$GD[GTindex, ]
            }
            
            
            if (!is.null(myFRG$GI))    {
              colnames(myFRG$GI) = c("SNP", "Chromosome", "Position")
              GI = as.matrix(myFRG$GI)
            }
            
            
            if (!is.null(myFRG$GI))    {
              numSNP = ncol(myFRG$GD)
            }  else{
              numSNP = 0
            }
            if (is.null(myFRG))
              numSNP = 0  #force to end the while loop
          } # end of if(!fullGD)
          
          if (fullGD)
            numSNP = 0  #force to end the while loop
          
          #Skip REML if xs is from a empty fragment file
          if (!is.null(xs))  {
            if (is.null(dim(xs)) ||
                nrow(xs) == 1)
              xs <- matrix(xs, length(xs), 1)
            
            xs <- as.matrix(xs)
            
            if (length(which(is.na(xs))) > 0) {
              #for the case where fragments are read in
              if (SNP.impute == "Major")
                xs[which(is.na(xs))] = 2
              if (SNP.impute == "Minor")
                xs[which(is.na(xs))] = 0
              if (SNP.impute == "Middle")
                xs[which(is.na(xs))] = 1
            }
            
            
            m <- ncol(xs) #number of SNPs
            t <- nrow(xs) #number of individuals
            
            Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Before cleaning")
            Memory = GAPIT.Memory(Memory = Memory, Infor = "Before cleaning")
            #allocate spaces for SNPs
            rm(dfs)
            rm(stats)
            rm(effect.est)
            rm(ps)
            rm(nobs)
            rm(maf)
            rm(rsquare_base)
            rm(rsquare)
            rm(df)
            rm(tvalue)
            rm(stderr)
            gc()
            
            Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "After cleaning")
            Memory = GAPIT.Memory(Memory = Memory, Infor = "After cleaning")
            
            dfs <- matrix(nrow = m, ncol = g)
            stats <- matrix(nrow = m, ncol = g)
            if (!Create.indicator)
              effect.est <- matrix(nrow = m, ncol = g)
            if (Create.indicator)
              effect.est <- NULL
            ps <- matrix(nrow = m, ncol = g)
            nobs <- matrix(nrow = m, ncol = g)
            maf <- matrix(nrow = m, ncol = g)
            rsquare_base <- matrix(nrow = m, ncol = g)
            rsquare <- matrix(nrow = m, ncol = g)
            df <- matrix(nrow = m, ncol = g)
            tvalue <- matrix(nrow = m, ncol = g)
            stderr <- matrix(nrow = m, ncol = g)
            #print(paste("Memory allocated.",sep=""))
            Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Memory allocation")
            Memory = GAPIT.Memory(Memory = Memory, Infor = "Memory allocation")
            
            if (optOnly)
              mloop = 0
            if (!optOnly)
              mloop = m
            
            #Loop on SNPs
            #print(paste("Number of SNPs is ",mloop," in genotype file ",file, sep=""))
            
            #set starting point of loop
            if (file == file.from & frag == 1) {
              loopStart = 0
            } else{
              loopStart = 1
            }
            
            for (i in loopStart:mloop) {
              #print(i)
              #--------------------------------------------------------------------------------------------------------------------<
              normalCase = TRUE
              
              if ((i > 0) &
                  (floor(i / 1000) == i / 1000)) {
                print(paste("Genotype file: ", file, ", SNP: ", i, " ", sep = ""))
              }
              # To extract current snp. It save computation for next one in case they are identical
              if (i == 0 & file == file.from & frag == 1) {
                #For the model without fitting SNP
                vids <- !is.na(ys[j, ]) #### Feng changed
                xv <- ys[j, vids] * 0 + 1 #### Feng changed
              }
              
              if (i > 0 | file > file.from | frag > 1) {
                if (Create.indicator) {
                  #I need create indicators and then calculate the minor allele frequency
                  condition.temp <- unique(xs[vids, i])
                  #Define what a bit is
                  
                  bit = nchar(as.character(xs[vids[1], i]))
                  
                  #Expand on the "which" statement below to include all instances of missing data
                  
                  if (bit == 1)
                    condition <-  condition.temp[-which(condition.temp == "N")]
                  if (bit == 2)
                    condition <-  condition.temp[-which(condition.temp == "NN")]
                  
                  #print("condition.temp is ")
                  #print(condition.temp)
                  
                  #print("condition is")
                  #print(condition)
                  
                  #print(paste("The value of i is ", i, sep = ""))
                  
                  
                  if (length(condition) <= 1) {
                    dfs[i,] <- rep(NA, g)
                    stats[i,] <- rep(NA, g)
                    effect.est <- rbind(effect.est, c(i, rep(NA, g), rep(NA, g)))
                    ps[i,] = rep(1, g)
                    rsquare[i,] <- rep(NA, g)
                    rsquare_base[i,] <- rep(NA, g)
                    maf[i,] <- rep(0, g)
                    df[i,] <- rep(NA, g)
                    tvalue[i,] <- rep(NA, g)
                    stderr[i,] <- rep(NA, g)
                    normalCase = FALSE
                    x.prev = vector(length = 0)
                  }
                  
                }
                if (normalCase) {
                  #print("The head of xs[vids,i] is")
                  #print(head(xs[vids,i]))
                  
                  if (Create.indicator) {
                    #I need create indicators and then calculate the minor allele frequency
                    
                    indicator <-
                      GAPIT.Create.Indicator(xs[vids, i], SNP.impute = SNP.impute)
                    xv <- indicator$x.ind
                    vids <- !is.na(xv[, 1]) #### Feng changed
                    
                    vids.TRUE = which(vids == TRUE)
                    vids.FALSE = which(vids == FALSE)
                    ns = nrow(xv)
                    ss = sum(xv[, ncol(xv)])
                    
                    maf[i] = min(ss / ns, 1 - ss / ns)
                    nobs[i] = ns
                    
                    q1 <-
                      q0 + ncol(xv)    # This is done so that parameter estimates for all indicator variables are included
                    
                    
                    #These two matrices need to be reinitiated for each SNP.
                    Xt <- matrix(NA, nr, q1)
                    iXX = matrix(NA, q1, q1)
                  }
                }
                
                if (!Create.indicator) {
                  #### Feng changed
                  #print(xs[1:10,1:10])
                  
                  xv <- xs[vids, i]
                  vids <- !is.na(xs[, i]) #### Feng changed
                  
                  vids.TRUE = which(vids == TRUE)
                  vids.FALSE = which(vids == FALSE)
                  ns = length(xv)
                  #print(xv))
                  ss = sum(xv)
                  
                  maf[i] = min(.5 * ss / ns, 1 - .5 * ss / ns)
                  nobs[i] = ns
                }
                
                nr <- sum(vids)
                if (i == 1 & file == file.from & frag == 1 &
                    !Create.indicator) {
                  Xt <- matrix(NA, nr, q1)
                  iXX = matrix(NA, q1, q1)
                }
                
              }
              
              #Situation of no variation for SNP except the fisrt one(synthetic for EMMAx/P3D)
              if ((min(xv) == max(xv)) & (i > 0 | file > file.from | frag > 1))
              {
                dfs[i,] <- rep(NA, g)
                stats[i,] <- rep(NA, g)
                if (!Create.indicator)
                  effect.est[i, ] <- rep(NA, g)
                if (Create.indicator)
                  effect.est <- rbind(effect.est, c(i, rep(NA, g), rep(NA, g)))
                ps[i,] = rep(1, g)
                rsquare[i,] <- rep(NA, g)
                rsquare_base[i,] <- rep(NA, g)
                df[i,] <- rep(NA, g)
                tvalue[i,] <- rep(NA, g)
                stderr[i,] <- rep(NA, g)
                normalCase = FALSE
              } else if (identical(x.prev, xv))
                #Situation of the SNP is identical to previous
              {
                if (i > 1 | file > file.from | frag > 1) {
                  dfs[i,] <- dfs[i - 1,]
                  stats[i,] <- stats[i - 1,]
                  if (!Create.indicator)
                    effect.est[i,] <- effect.est[i - 1,]
                  if (Create.indicator)
                    effect.est <-
                    rbind(effect.est, c(i, rep(NA, g), rep(NA, g))) #If the previous SNP is idnetical, indicate this by "NA"
                  ps[i,] <- ps[i - 1,]
                  rsquare[i,] <- rsquare[i - 1,]
                  rsquare_base[i,] <- rsquare_base[i - 1,]
                  df[i,] <- df[i - 1,]
                  tvalue[i,] <- tvalue[i - 1,]
                  stderr[i,] <- stderr[i - 1,]
                  normalCase = FALSE
                }
              }
              #-------------------------------------------------------------------------------------------------------------------->
              if (i == 0 & file == file.from & frag == 1) {
                #Calculate the log likelihood function for the intercept only model
                
                #vids <- !is.na(ys[j,])
                yv <- ys[j, vids]
                
                X.int <- matrix(1, nrow(as.matrix(yv)), ncol(as.matrix(yv)))
                iX.intX.int <- solve(crossprod(X.int, X.int))
                iX.intY <- crossprod(X.int, as.matrix(as.matrix(yv)))
                beta.int <-
                  crossprod(iX.intX.int, iX.intY)  #Note: we can use crossprod here becase iXX is symmetric
                X.int.beta.int <- X.int %*% beta.int
                
                
                
                
                #X.int <- matrix(1,nrow(as.matrix(ys[!is.na(ys)])),ncol(as.matrix(ys[!is.na(ys)])))
                #iX.intX.int <- solve(crossprod(X.int, X.int))
                #iX.intY <- crossprod(X.int, as.matrix(ys[!is.na(ys)]))
                #beta.int <- crossprod(iX.intX.int, iX.intY)  #Note: we can use crossprod here becase iXX is symmetric
                #X.int.beta.int <- X.int%*%beta.int
                
                logL0 <- 0.5 * ((-length(yv)) * log(((2 * pi) / length(yv))
                                                    * crossprod((yv - X.int.beta.int), (yv - X.int.beta.int)
                                                    ))
                                - length(yv))
                
                
                #logL0 <- 0.5*((-nrow(as.matrix(ys[!is.na(ys)])))*log(((2*pi)/nrow(ys))
                # *crossprod(((as.matrix(ys[!is.na(ys)]))-X.int.beta.int),((as.matrix(ys[!is.na(ys)]))-X.int.beta.int)))
                # -nrow(as.matrix(ys[!is.na(ys)])))
                
                #print(paste("The value of logL0 inside of the calculating SNPs loop is", logL0, sep = ""))
              }
              
              #Normal case
              if (normalCase)
              {
                #--------------------------------------------------------------------------------------------------------------------<
                #nv <- sum(vids)
                yv <- ys[j, vids] #### Feng changed
                nr <- sum(vids) #### Feng changed
                if (!is.null(Z) & !is.null(K))
                {
                  r <-
                    ncol(Z) ####Feng, add a variable to indicate the number of random effect
                  vran <-
                    vids[1:r] ###Feng, add a variable to indicate random effects with nonmissing genotype
                  tv <- sum(vran)  #### Feng changed
                }
                
                
                
                #-------------------------------------------------------------------------------------------------------------------->
                
                #--------------------------------------------------------------------------------------------------------------------<
                
                
                
                if (i > 0 | file > file.from | frag > 1)
                  dfs[i, j] <- nr - q1
                if (i > 0 | file > file.from | frag > 1) {
                  if (!Create.indicator)
                    X <- cbind(X0[vids, , drop = FALSE], xs[vids, i])
                  if (Create.indicator) {
                    X <- cbind(X0[vids, , drop = FALSE], xv)
                    #if(i == 1) {print("the head of X for running GWAS is")}
                    #if(i == 1) {print(head(X))}
                  }
                  
                }
                #Recalculate eig and REML if not using P3D  NOTE THIS USED TO BE BEFORE the two solid lines
                if (SNP.P3D == FALSE & !is.null(K))
                {
                  if (!is.null(Z))
                    eig.R <-
                      emma.eigen.R.w.Z(Z, K, X) #This will be used to get REstricted ML (REML)
                  if (is.null(Z))
                    eig.R <-
                      emma.eigen.R.wo.Z(K, X) #This will be used to get REstricted ML (REML)
                  if (!is.null(Z))
                    REMLE <-
                      GAPIT.emma.REMLE(ys[j, ], X, K, Z, ngrids, llim, ulim, esp, eig.R)
                  if (is.null(Z))
                    REMLE <-
                      GAPIT.emma.REMLE(ys[j, ], X, K, Z = NULL, ngrids, llim, ulim, esp, eig.R)
                  if (!is.null(Z) &
                      !is.null(K))
                    U <-
                      eig.L$vectors * matrix(c(sqrt(
                        1 / (eig.L$values + REMLE$delta)
                      ), rep(
                        sqrt(1 / REMLE$delta), nr - tv
                      )), nr, ((nr - tv) + length(eig.L$values)), byrow = TRUE)
                  if (is.null(Z) &
                      !is.null(K))
                    U <-
                      eig.L$vectors * matrix(sqrt(1 / (
                        eig.L$values + REMLE$delta
                      )), nr, length(eig.L$values), byrow = TRUE)
                  
                  vgs <- REMLE$vg
                  ves <- REMLE$ve
                  REMLs <- REMLE$REML
                  REMLE_delta = REMLE$delta
                  
                }
                
                if (n == nr)
                {
                  if (!is.null(K))
                  {
                    yt <- crossprod(U, yv)
                    if (i == 0 & file == file.from & frag == 1) {
                      X0t <- crossprod(U, X0)
                      Xt <- X0t
                    }
                    if (i > 0 | file > file.from | frag > 1) {
                      #if(i ==1 & file==file.from&frag==1) Xt <- matrix(NA,nr, q1)
                      
                      if (Create.indicator) {
                        xst <- crossprod(U, X[, (q0 + 1):q1])
                        Xt[1:nr, 1:q0] <- X0t
                        Xt[1:nr, (q0 + 1):q1] <- xst
                        
                      }
                      
                      #print(paste("i:",i,"q0:",q0,"q1:",q1,"nt:",nr,"XT row",nrow(Xt),"XT col",ncol(Xt),sep=" "))
                      if (!Create.indicator) {
                        xst <- crossprod(U, X[, ncol(X)])
                        Xt[1:nr, 1:q0] <- X0t
                        Xt[1:nr, q1] <- xst
                      }
                    }
                  } else{
                    yt = yv
                    if (i == 0 & file == file.from & frag == 1)
                      X0t <- X0
                    if (i > 0 | file > file.from | frag > 1)
                      xst <- X[, ncol(X)]
                  }
                  
                  if (i == 0 & file == file.from & frag == 1) {
                    X0X0 <- crossprod(X0t, X0t)
                    #XX <- X0X0
                  }
                  if (i > 0 | file > file.from | frag > 1) {
                    #if(i == 1)XX=matrix(NA,q1,q1)
                    
                    
                    X0Xst <- crossprod(X0t, xst)
                    XstX0 <- t(X0Xst)
                    xstxst <- crossprod(xst, xst)
                    # if(i == 1){
                    # Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Calculate_X0Xst_XstX0_xstxst")
                    # Memory=GAPIT.Memory(Memory=Memory,Infor="Calculate_X0Xst_XstX0_xstxst")
                    # }
                    #XX <- rbind(cbind(X0X0, X0Xst), cbind(XstX0, xstxst))
                    
                    #XX[1:q0,1:q0] <- X0X0
                    #XX[q1,1:q0] <- X0Xst
                    #XX[1:q0,q1] <- X0Xst
                    #XX[q1,q1] <- xstxst
                    
                    
                  }
                  
                  
                  if (X0X0[1, 1] == "NaN")
                  {
                    Xt[which(Xt == "NaN")] = 0
                    yt[which(yt == "NaN")] = 0
                    XX = crossprod(Xt, Xt)
                  }
                  if (i == 0 & file == file.from & frag == 1) {
                    X0Y <- crossprod(X0t, yt)
                    XY <- X0Y
                  }
                  if (i > 0 | file > file.from | frag > 1) {
                    xsY <- crossprod(xst, yt)
                    XY <- c(X0Y, xsY)
                    # if(i == 1){
                    # Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Calculate_xsY_X0Y")
                    # Memory=GAPIT.Memory(Memory=Memory,Infor="Calculate_xsY_X0Y")
                    # }
                  }
                  #XY = crossprod(Xt,yt)
                }
                
                #Missing SNP
                if (n > nr)
                {
                  UU = crossprod(U, U)
                  A11 = UU[vids.TRUE, vids.TRUE]
                  A12 = UU[vids.TRUE, vids.FALSE]
                  A21 = UU[vids.FALSE, vids.TRUE]
                  A22 = UU[vids.FALSE, vids.FALSE]
                  A22i = try(solve(A22), silent = TRUE)
                  if (inherits(A22i, "try-error"))
                    A22i <- ginv(A22)
                  
                  F11 = A11 - A12 %*% A22i %*% A21
                  XX = crossprod(X, F11) %*% X
                  XY = crossprod(X, F11) %*% yv
                }
                if (i == 0 & file == file.from & frag == 1) {
                  iX0X0 <- try(solve(X0X0), silent = TRUE)
                  if (inherits(iX0X0, "try-error")) {
                    iX0X0 <- ginv(X0X0)
                    print(
                      "At least two of your covariates are linearly dependent. Please reconsider the covariates you are using for GWAS and GPS"
                    )
                  }
                  iXX <- iX0X0
                }
                if (i > 0 | file > file.from | frag > 1) {
                  #if(i ==1 &file==file.from &frag==1) iXX=matrix(NA,q1,q1)
                  if (Create.indicator) {
                    B22 <- xstxst - XstX0 %*% iX0X0 %*% X0Xst
                    invB22 <- solve(B22)
                    B21 <- tcrossprod(XstX0, iX0X0)
                    NeginvB22B21 <- crossprod(-invB22, B21)
                    B11 <- iX0X0 + as.numeric(invB22) * crossprod(B21, B21)
                    
                    
                    
                    iXX[1:q0, 1:q0] = B11
                    iXX[(q0 + 1):q1, (q0 + 1):q1] = solve(B22)
                    iXX[(q0 + 1):q1, 1:q0] = NeginvB22B21
                    iXX[1:q0, (q0 + 1):q1] = t(NeginvB22B21)
                    
                  }
                  
                  
                  if (!Create.indicator) {
                    B22 <- xstxst - XstX0 %*% iX0X0 %*% X0Xst
                    invB22 <- 1 / B22
                    #B12 <- crossprod(iX0X0,X0Xst)
                    B21 <- tcrossprod(XstX0, iX0X0)
                    NeginvB22B21 <- crossprod(-invB22, B21)
                    #B11 <- iX0X0 + B12%*%invB22%*%B21
                    B11 <- iX0X0 + as.numeric(invB22) * crossprod(B21, B21)
                    #iXX <- rbind(cbind(B11,t(NeginvB22B21)), cbind(NeginvB22B21,invB22))
                    
                    iXX[1:q0, 1:q0] = B11
                    iXX[q1, q1] = 1 / B22
                    iXX[q1, 1:q0] = NeginvB22B21
                    iXX[1:q0, q1] = NeginvB22B21
                    
                    
                  }
                  #if(i == 1){
                  # Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Calculate_iXX")
                  # Memory=GAPIT.Memory(Memory=Memory,Infor="Calculate_iXX")
                  #}
                  
                }
                
                if (is.null(K)) {
                  iXX <- try(solve(crossprod(X, X)), silent = TRUE)
                  if (inherits(iXX, "try-error"))
                    iXX <- ginv(crossprod(X, X))
                  XY = crossprod(X, yv)
                }
                
                #iXX <- try(solve(XX))
                #if(inherits(iXX, "try-error")) iXX <- ginv(crossprod(Xt, Xt))
                #print("The dimension if iXX is")
                #print(dim(iXX))
                #print("The length of XY is")
                #print(length(XY))
                
                beta <-
                  crossprod(iXX, XY) #Note: we can use crossprod here becase iXX is symmetric
                #print("beta was estimated")
                
                #-------------------------------------------------------------------------------------------------------------------->
                
                
                #--------------------------------------------------------------------------------------------------------------------<
                if (i == 0 & file == file.from & frag == 1 & !is.null(K))
                {
                  Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "ReducedModel")
                  Memory = GAPIT.Memory(Memory = Memory, Infor = "ReducdModel")
                  
                  #beta.cv=beta
                  
                  
                  
                  XtimesBetaHat <- X %*% beta
                  
                  YminusXtimesBetaHat <- ys[j, ] - XtimesBetaHat
                  vgK <- vgs * K
                  Dt <- crossprod(U, YminusXtimesBetaHat)
                  
                  if (!is.null(Z))
                    Zt <- crossprod(U, Z)
                  if (is.null(Z))
                    Zt <- t(U)
                  
                  if (X0X0[1, 1] == "NaN")
                  {
                    Dt[which(Dt == "NaN")] = 0
                    Zt[which(Zt == "NaN")] = 0
                  }
                  
                  BLUP <-
                    K %*% crossprod(Zt, Dt) #Using K instead of vgK because using H=V/Vg
                  
                  #print("!!!!")
                  #Clean up the BLUP starf to save memory
                  Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "before Dt clean")
                  Memory = GAPIT.Memory(Memory = Memory, Infor = "before Dt clean")
                  rm(Dt)
                  gc()
                  Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Dt clean")
                  Memory = GAPIT.Memory(Memory = Memory, Infor = "Dt clean")
                  
                  
                  
                  
                  
                  grand.mean.vector <- rep(beta[1], length(BLUP))
                  BLUP_Plus_Mean <- grand.mean.vector + BLUP
                  Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "BLUP")
                  Memory = GAPIT.Memory(Memory = Memory, Infor = "BLUP")
                  
                  #PEV
                  C11 = try(vgs * solve(crossprod(Xt, Xt)), silent = TRUE)
                  if (inherits(C11, "try-error"))
                    C11 = vgs * ginv(crossprod(Xt, Xt))
                  
                  C21 = -K %*% crossprod(Zt, Xt) %*% C11
                  Kinv = try(solve(K)  , silent = TRUE)
                  if (inherits(Kinv, "try-error"))
                    Kinv = ginv(K)
                  
                  if (!is.null(Z))
                    term.0 = crossprod(Z, Z) / ves
                  if (is.null(Z))
                    term.0 = diag(1 / ves, nrow(K))
                  
                  term.1 = try(solve(term.0 + Kinv / vgs) , silent = TRUE)
                  if (inherits(term.1, "try-error"))
                    term.1 = ginv(term.0 + Kinv / vgs)
                  
                  term.2 = C21 %*% crossprod(Xt, Zt) %*% K
                  C22 = (term.1 - term.2)
                  PEV = as.matrix(diag(C22))
                  #print(paste("The value of is.na(CVI) is", is.na(CVI),  sep = ""))
                  if (!is.na(CVI)) {
                    XCV = as.matrix(cbind(1, data.frame(CVI[, -1])))
                    
                    #CV.Inheritance specified
                    beta.Inheritance = beta
                    if (!is.null(CV.Inheritance)) {
                      XCV = XCV[, 1:(1 + CV.Inheritance)]
                      beta.Inheritance = beta[1:(1 + CV.Inheritance)]
                    }
                    #Interception only
                    if (length(beta) == 1)
                      XCV = X
                    
                    BLUE = try(XCV %*% beta.Inheritance, silent = TRUE)
                    if (inherits(BLUE, "try-error"))
                      BLUE = NA
                    #print("GAPIT just after BLUE")
                    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "PEV")
                    Memory = GAPIT.Memory(Memory = Memory, Infor = "PEV")
                    
                  }#end of if(i ==0&file==file.from   & !is.null(K))
                  if (is.na(CVI))
                    BLUE = NA
                }#end if(!is.na(CVI))
                #-------------------------------------------------------------------------------------------------------------------->
                
                #--------------------------------------------------------------------------------------------------------------------<
                if (i == 0 & file == file.from & frag == 1 & is.null(K))
                {
                  YY = crossprod(yt, yt)
                  ves = (YY - crossprod(beta, XY)) / (n - q0)
                  r = yt - X %*% iXX %*% XY
                  REMLs = -.5 * (n - q0) * log(det(ves)) - .5 * n - .5 * (n - q0) *
                    log(2 * pi)
                  # REMLs=-.5*n*log(det(ves)) -.5*log(det(iXX)/ves) -.5*crossprod(r,r)/ves -.5*(n-q0)*log(2*pi)
                  vgs = 0
                  BLUP = 0
                  BLUP_Plus_Mean = NaN
                  PEV = ves
                  #print(paste("X row:",nrow(X)," col:",ncol(X)," beta:",length(beta),sep=""))
                  XCV = as.matrix(cbind(1, data.frame(CVI[, -1])))
                  
                  #CV.Inheritance specified
                  beta.Inheritance = beta
                  if (!is.null(CV.Inheritance)) {
                    XCV = XCV[, 1:(1 + CV.Inheritance)]
                    beta.Inheritance = beta[1:(1 + CV.Inheritance)]
                  }
                  #Interception only
                  if (length(beta) == 1)
                    XCV = X
                  
                  
                  #BLUE=XCV%*%beta.Inheritance   modified by jiabo wang 2016.11.21
                  BLUE = try(XCV %*% beta.Inheritance, silent = TRUE)
                  if (inherits(BLUE, "try-error"))
                    BLUE = NA
                  
                }
                
                
                #Clean up the BLUP stuff to save memory
                if (i == 0 & file == file.from & frag == 1 & !is.null(K))
                {
                  Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "K normal")
                  Memory = GAPIT.Memory(Memory = Memory, Infor = "K normal")
                  if (SNP.P3D == TRUE)
                    K = 1  #NOTE: When SNP.P3D == FALSE, this line will mess up the spectral decomposition of the kinship matrix at each SNP.
                  rm(Dt)
                  rm(Zt)
                  rm(Kinv)
                  rm(C11)
                  rm(C21)
                  rm(C22)
                  
                  gc()
                  Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "K set to 1")
                  Memory = GAPIT.Memory(Memory = Memory, Infor = "K set to 1")
                }
                
                if (i == 0 & file == file.from & frag == 1) {
                  beta.cv = beta
                  X.beta <- X %*% beta
                  
                  if (!is.null(K)) {
                    U.times.yv.minus.X.beta <- crossprod(U, (yv - X.beta))
                    logLM_Base <-
                      0.5 * (-length(yv) * log(((2 * pi) / length(yv)) * crossprod(
                        U.times.yv.minus.X.beta,
                        U.times.yv.minus.X.beta
                      )
                      )
                      - sum(log(eig.full.plus.delta)) - length(yv))
                    
                  }
                  if (is.null(K)) {
                    U.times.yv.minus.X.beta <- yv - X.beta
                    logLM_Base <-
                      0.5 * (-length(yv) * log(((2 * pi) / length(yv)) * crossprod(
                        U.times.yv.minus.X.beta,
                        U.times.yv.minus.X.beta
                      )
                      ) - length(yv))
                  }
                  rsquare_base_intitialized <-
                    1 - exp(-(2 / length(yv)) * (logLM_Base - logL0))
                  
                }
                
                
                #print(Create.indicator)
                #calculate t statistics and P-values
                if (i > 0 | file > file.from | frag > 1)
                {
                  if (!Create.indicator) {
                    #if(i<5)print(beta[q1])
                    #if(i<5)print(iXX[q1, q1])
                    if (!is.null(K))
                      stats[i, j] <- beta[q1] / sqrt(iXX[q1, q1] * vgs)
                    if (is.null(K))
                      stats[i, j] <- beta[q1] / sqrt(iXX[q1, q1] * ves)
                    effect.est[i,] <- beta[q1]
                    ps[i,] <-
                      2 * pt(abs(stats[i,]), dfs[i,], lower.tail = FALSE)
                    if (is.na(ps[i, ]))
                      ps[i, ] = 1
                    #print(c(i,ps[i,],stats[i,],beta[q1],iXX[q1, q1]))
                  }
                  if (Create.indicator) {
                    F.num.first.two <-
                      crossprod(beta[(q0 + 1):q1], solve(iXX[(q0 + 1):q1, (q0 + 1):q1]))
                    if (!is.null(K))
                      stats[i, j] <-
                        (F.num.first.two %*% beta[(q0 + 1):q1]) / (length((q0 + 1):q1) * vgs)
                    if (is.null(K))
                      stats[i, j] <-
                        (F.num.first.two %*% beta[(q0 + 1):q1]) / (length((q0 + 1):q1) * ves)
                    effect.est <-
                      rbind(effect.est,
                            cbind(
                              rep(i, length((
                                q0 + 1
                              ):q1)),
                              indicator$unique.SNPs,
                              beta[(q0 + 1):q1]
                            )) #Replace with rbind
                    ps[i,] <-
                      pf(
                        stats[i, j],
                        df1 = length((q0 + 1):q1),
                        df2 = (nr - ncol(X)),
                        lower.tail = FALSE
                      ) #Alex, are these denominator degrees of freedom correct?
                    dfs[i, ] <- nr - nrow(X)
                    
                  }
                  #Calculate the maximum full likelihood function value and the r square
                  
                  X.beta <- X %*% beta
                  if (!is.null(K)) {
                    U.times.yv.minus.X.beta <- crossprod(U, (yv - X.beta))
                    logLM <-
                      0.5 * (-length(yv) * log(((2 * pi) / length(yv)) * crossprod(
                        U.times.yv.minus.X.beta,
                        U.times.yv.minus.X.beta
                      )
                      )
                      - sum(log(eig.full.plus.delta)) - length(yv))
                  }
                  if (is.null(K)) {
                    U.times.yv.minus.X.beta <- yv - X.beta
                    logLM <-
                      0.5 * (-length(yv) * log(((2 * pi) / length(yv)) * crossprod(
                        U.times.yv.minus.X.beta,
                        U.times.yv.minus.X.beta
                      )
                      ) - length(yv))
                  }
                  
                  rsquare_base[i,] <- rsquare_base_intitialized
                  rsquare[i,] <- 1 - exp(-(2 / length(yv)) * (logLM - logL0))
                  
                  #Calculate df, t value and standard error _xiaolei changed
                  df[i, ] <- dfs[i, ]
                  
                  tvalue[i, ] <- stats[i, j]
                  stderr[i, ] <- beta[ncol(CVI) + 1] / stats[i, j]
                  #stderr[i,] <- sqrt(vgs)
                  # modified by Jiabo at 20191115
                }
                #print("!!!!!!!!!!!!!!!")
                #print(Create.indicator)
                #-------------------------------------------------------------------------------------------------------------------->
                
              } # End of if(normalCase)
              x.prev = xv #update SNP
              
            } # End of loop on SNPs
            
            Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Screening SNPs")
            Memory = GAPIT.Memory(Memory = Memory, Infor = "Screening SNPs")
            # print(head(tvalue))
            # print(head(stderr))
            # print(head(effect.est))
            #output p value for the genotype file
            if (!fullGD)
            {
              #print("!!!!!!!!!!")
              #print(dim(GI))
              write.table(
                GI,
                paste(
                  "GAPIT.TMP.GI.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = TRUE
              )
              write.table(
                ps,
                paste(
                  "GAPIT.TMP.ps.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE
              )
              write.table(
                maf,
                paste(
                  "GAPIT.TMP.maf.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE
              )
              write.table(
                nobs,
                paste(
                  "GAPIT.TMP.nobs.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE
              )
              write.table(
                rsquare_base,
                paste(
                  "GAPIT.TMP.rsquare.base.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE
              )
              write.table(
                rsquare,
                paste(
                  "GAPIT.TMP.rsquare.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE
              )
              write.table(
                df,
                paste(
                  "GAPIT.TMP.df.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE
              )
              write.table(
                tvalue,
                paste(
                  "GAPIT.TMP.tvalue.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE
              )
              write.table(
                stderr,
                paste(
                  "GAPIT.TMP.stderr.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE
              )
              write.table(
                effect.est,
                paste(
                  "GAPIT.TMP.effect.est.",
                  name.of.trait,
                  file,
                  ".",
                  frag,
                  ".txt",
                  sep = ""
                ),
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE
              )
              
              #rm(dfs,stats,ps,nobs,maf,GI)   #This cause problem on return
              #gc()
            }
            
            frag = frag + 1   #Progress to next fragment
            
          } #end of if(!is.null(X))
          
        } #end of repeat on fragment
        
        
        
      } # Ebd of loop on file
    } # End of loop on traits
    
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "GWAS done for this Trait")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "GWAS done for this Trait")
    #print("GAPIT.EMMAxP3D accomplished successfully!")
    
    return(
      list(
        ps = ps,
        REMLs = -2 * REMLs,
        stats = stats,
        effect.est = effect.est,
        rsquare_base = rsquare_base,
        rsquare = rsquare,
        dfs = dfs,
        df = df,
        tvalue = tvalue,
        stderr = stderr,
        maf = maf,
        nobs = nobs,
        Timmer = Timmer,
        Memory = Memory,
        vgs = vgs,
        ves = ves,
        BLUP = BLUP,
        BLUP_Plus_Mean = BLUP_Plus_Mean,
        PEV = PEV,
        BLUE = BLUE,
        logLM = logLM,
        effect.snp = effect.est,
        effect.cv = beta.cv
      )
    )
    
  }#end of GAPIT.EMMAxP3D function

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure` <-
  function(PWI = PWI,
           FDR.Rate = 0.05,
           FDR.Procedure = "BH") {
    #Object: Conduct the Benjamini-Hochberg FDR-Controlling Procedure
    #Output: PWIP, number.of.significant.SNPs
    #Authors: Alex Lipka and Zhiwu Zhang
    # Last update: May 5, 2011
    ##############################################################################################
    #Make sure that your compouter has the latest version of Bioconductor (the "Biobase" package) and multtest
    
    if (is.null(PWI))
    {
      PWIP = NULL
      number.of.significant.SNPs = 0
    }
    
    if (!is.null(PWI))
    {
      #library(multtest)
      
      if (dim(PWI)[1] == 1) {
        PWIP <- cbind(PWI, PWI[4])
        colnames(PWIP)[9] <- "FDR_Adjusted_P-values"
      }
      
      if (dim(PWI)[1] > 1) {
        #mt.rawp2adjp Performs the Simes procedure.  The output should be two columns, Left column: originial p-value
        #Right column: Simes corrected p-value
        res <- mt.rawp2adjp(PWI[, 4], FDR.Procedure)
        
        #This command should order the p-values in the order of the SNPs in the data set
        adjp <- res$adjp[order(res$index),]
        
        #round(adjp[1:7,],4)
        #Logical statment: 0, if Ho is not rejected; 1, if  Ho is rejected, by the Simes corrected p-value
        #  temp <- mt.reject(adjp[,2], FDR.Rate)
        
        #Lists all number of SNPs that were rejected by the BY procedure
        #temp$r
        
        #Attach the FDR adjusted p-values to AS_Results
        
        PWIP <- cbind(PWI, adjp[, 2])
        
        #Sort these data by lowest to highest FDR adjusted p-value
        PWIP <- PWIP[order(PWIP[, 4]), ]
        
        colnames(PWIP)[9] <- "FDR_Adjusted_P-values"
        #  number.of.significant.SNPs = temp$r
      }
      #print("GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure accomplished successfully!")
    }
    #return(list(PWIP=PWIP, number.of.significant.SNPs = number.of.significant.SNPs))
    return(list(PWIP = PWIP))
  }#GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure ends here

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.Power.Compare` <-

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.Power` <-
  function(WS = c(1e0, 1e3, 1e4, 1e5, 1e6, 1e7),
           GM = NULL,
           seqQTN = NULL,
           GWAS = NULL,
           maxOut = 100,
           alpha = c(.01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1),
           MaxBP = 1e10) {
    #Object: To evaluate power and FDR for the top (maxOut) positive interval defined by WS
    #Input: WS- window size
    #Input: GM - m by 3  matrix for SNP name, chromosome and BP
    #Input: seqQTN - s by 1 vecter for index of QTN on GM (+1 for GDP column wise)
    #Input: GWAS- SNP,CHR,BP,P,MAF
    #maxOut: maximum number of rows to report
    #Requirement: None
    #Output: Table and Plots
    #Authors: Zhiwu Zhang
    # Date  start: April 2, 2013
    # Last update: April 2, 2013
    ##############################################################################################
    #print("GAPIT.Power Started")
    if (is.null(seqQTN) |
        is.null(GM) |
        is.null(GWAS))
      return(list(
        FDR = NULL,
        Power = NULL,
        Power.Alpha = NULL,
        alpha = NULL
      ))
    
    #-----------------FDR and Power analysis-------------------------
    #Information needed: myGAPIT$GWAS,myGM and QTN(r)
    nWin = matrix(NA, length(WS), 1)
    
    format_GWAS = cbind(GWAS[, 1:4], NA, NA, NA)
    
    names(format_GWAS) <-
      c("SNP",
        "Chromosome",
        "Position",
        "P.value",
        "maf",
        "nobs",
        "FDR_Adjusted_P-values")
    myGM = GM
    
    #loop window size here
    
    theWS = 1
    for (theWS in 1:length(WS)) {
      ws = WS[theWS]
      #Label QTN intervals
      #Restore original order
      #QTNList=r-1
      QTNList = seqQTN
      myGM2 = cbind(myGM, rep(0, nrow(myGM)), 1:nrow(myGM), NA) #Initial QTN status as 0
      
      
      #Extract QTN positions
      myGM2[, 6] = floor((
        as.numeric(as.character(myGM2[, 2])) * MaxBP + as.numeric(as.character(myGM2[, 3]))
      ) / ws) #Label QTN as 1
      
      QTNInterval = myGM2[QTNList, 6]
      thePosition = myGM2[, 6] %in% QTNInterval
      
      myGM2[thePosition, 4] = 1 #Label QTN as 1
      names(myGM2) <- c("SNP", "Chromosome", "Position", "QTN", "Seq")
      
      #Merge to P vlaues
      #GWAS<- merge(myGAPIT$GWAS[,1:7],myGM2[,c(1,4,5)],by="SNP")
      GWAS <-
        merge(format_GWAS[, 1:7], myGM2[, c(1, 4, 5)], by = "SNP")#xiaoalei changed
      
      #checking
      #zw=GWAS[order(GWAS[,4],decreasing = FALSE),]
      #zw=GWAS[order(GWAS[,8],decreasing = TRUE),]
      #head(zw)
      
      #Creat windows
      myQTN = GAPIT.Specify(
        GI = GWAS[, 1:3],
        GP = GWAS,
        bin.size = ws,
        MaxBP = MaxBP
      )
      QTN = GWAS[myQTN$index, ]
      
      #Calculate alpha
      qtnLoc = which(QTN[, 8] == 1) #get the position of QTN
      P.QTN = QTN[qtnLoc, 4] #p value of QTN
      P.marker = QTN[-qtnLoc, 4] #p value of non qtn (marker)
      cutOff = matrix(quantile(P.marker, alpha, na.rm = TRUE), ncol = 1)#xiaoalei changed
      myPower.Alpha = apply(cutOff, 1, function(x) {
        Power = length(which(P.QTN < x)) / length(P.QTN)
      })
      
      
      #Sort on P
      #QTN=QTN[order(as.numeric(as.character(QTN[,3])),decreasing = FALSE),]
      #QTN=QTN[order(as.numeric(as.character(QTN[,2])),decreasing = FALSE),]
      QTN = QTN[order(as.numeric(as.character(QTN[, 4])), decreasing = FALSE), ]
      names(QTN) <-
        c("SNP",
          "Chromosome",
          "Position",
          "P",
          "FDR",
          "Power",
          "Order",
          "QTN",
          "Seq")
      
      #calculate power
      QTN[, 7] = 1:nrow(QTN)
      QTN[, 5] = cumsum(1 - QTN[, 8]) / QTN[, 7]   #FDR
      QTN[, 6] = cumsum(QTN[, 8]) / sum(QTN[, 8]) #Power
      
      #Save results
      if (theWS == 1) {
        nWin = matrix(NA, length(WS), 1)
        FDR = array(NA, dim = c(nrow(QTN), length(WS)))
        Power = array(NA, dim = c(nrow(QTN), length(WS)))
        Power.Alpha = array(NA, dim = c(length(alpha), length(WS)))
      }
      
      nWin[theWS] = nrow(QTN)
      FDR[1:nWin[theWS], theWS] = QTN[, 5]
      Power[1:nWin[theWS], theWS] = QTN[, 6]
      Power.Alpha[, theWS] = myPower.Alpha
      
    }#end of window size loop
    nOut = min(maxOut, max(nWin))
    index = 1:nOut
    return(list(
      FDR = FDR[index, ],
      Power = Power[index, ],
      Power.Alpha = Power.Alpha,
      alpha = alpha
    ))
  }#end of GAPIT.Power

#-----------------------------------------------------------------------
# Renamed rm base fuction to avoid warnmings
#-----------------------------------------------------------------------
rm <- function (...) {
	suppressWarnings (base::rm(...))
}

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
emma.eigen.L <- function(Z, K, complete = TRUE) {
  if (is.null(Z)) {
    return(emma.eigen.L.wo.Z(K))
  }
  else {
    return(emma.eigen.L.w.Z(Z, K, complete))
  }
}

emma.eigen.L.wo.Z <- function(K) {
  eig <- eigen(K, symmetric = TRUE)
  return(list(values = eig$values, vectors = eig$vectors))
}

emma.eigen.L.w.Z <- function(Z, K, complete = TRUE) {
  if (complete == FALSE) {
    vids <- colSums(Z) > 0
    Z <- Z[, vids]
    K <- K[vids, vids]
  }
  eig <- eigen(K %*% crossprod(Z, Z),
               symmetric = FALSE,
               EISPACK = TRUE)
  return(list(
    values = eig$values,
    vectors = qr.Q(qr(Z %*% eig$vectors), complete = TRUE)
  ))
}

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
`GAPIT.GS` <-
  function(KW, KO, KWO, GAU, UW) {
    #Object: to derive BLUP for the individuals without phenotype
    #UW:BLUP and PEV of ID with phenotyp
    #Output: BLUP
    #Authors: Zhiwu Zhang
    # Last update: Oct 22, 2015  by Jiabo Wang
    ##############################################################################################
    #print(dim(UW))
    UO = try(t(KWO) %*% solve(KW) %*% UW, silent = TRUE)
    #print(dim(KWO)) #kinship without inference
    #print(dim(KW))  #kinship within inference
    # print(dim(UW))  #BLUP AND PEV of reference
    
    if (inherits(UO, "try-error"))
    {
      GTT = try(t(KWO) %*% ginv(as.matrix(KW)) %*% UW)
      if (inherits(GTT, "try-error"))
      {
        write.csv(KW, "KW.csv", quote = F, row.names = F)
        KW = read.csv("KW.csv", head = T)
        UO = t(KWO) %*% ginv(as.matrix(KW)) %*% UW
        system("rm KW.csv")
      } else{
        UO = GTT
      }
    }
    n = ncol(UW) #get number of columns, add additional for individual name
    
    BLUP = data.frame(as.matrix(GAU[, 1:4]))
    BLUP.W = BLUP[which(GAU[, 3] < 2), ]
    W_BLUP = BLUP.W[order(as.numeric(as.matrix(BLUP.W[, 4]))), ]
    UW = UW[which(rownames(UW) == colnames(KW)), ] # get phenotype groups order
    
    ID.W = as.numeric(as.matrix(W_BLUP[, 4]))
    n.W = max(ID.W)
    DS.W = diag(n.W)[ID.W, ]
    # print(dim(DS.W))
    # print(dim(UW))
    ind.W = DS.W %*% UW
    
    all.W = cbind(W_BLUP, ind.W)
    all = all.W
    
    BLUP.O = BLUP[which(GAU[, 3] == 2), ]
    O_BLUP = BLUP.O[order(as.numeric(as.matrix(BLUP.O[, 4]))), ]
    #print(dim(O_BLUP))
    if (nrow(O_BLUP) > 0) {
      ID.O = as.numeric(as.matrix(O_BLUP[, 4]))
      n.O = max(ID.O)
      DS.O = diag(n.O)[ID.O, ]
      ind.O = DS.O %*% UO
      all.O = cbind(O_BLUP, ind.O)
      all = rbind(all.W, all.O)
    }
    
    colnames(all) = c("Taxa", "Group", "RefInf", "ID", "BLUP", "PEV")
    
    print("GAPIT.GS accomplished successfully!")
    return(list(BLUP = all))
  }#The function GAPIT.GS ends here

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#main()
