#!/usr/bin/Rscript
#
# Create multiGWAS profile according to current directory
# 

MULTIGWAS_HOME = strsplit (getwd (), "/install")[[1]][1]
message ("MULTIGWAS_HOME=",MULTIGWAS_HOME)

message ("Creating multiGWAS profile...")
sink ("multiGWAS_profile.sh", append=F)
writeLines ("\n#------------------- multiGWAS.R profile ---------------------")
writeLines (paste0 ("export MULTIGWAS_HOME=", MULTIGWAS_HOME))
writeLines ("MULTIGWAS_TOOLS=$MULTIGWAS_HOME/opt/tools")
writeLines ("MULTIGWAS_MAIN=$MULTIGWAS_HOME/main")
writeLines ("export PATH=$PATH:$MULTIGWAS_TOOLS:$MULTIGWAS_MAIN")
sink()

# Write into .bashrc
profileFile = paste0 (path.expand ("~"), "/.bashrc")
sink (profileFile, append=T)
writeLines ("\n#------------------- multiGWAS.R tool profile ---------------------")
writeLines (paste0 (". ", MULTIGWAS_HOME, "/multiGWAS_profile.sh"))
sink ()

message ("\nMultiGWAS is ready to use, right after installed!\n")

#---------------------- # Install R libraries----------------------------------
libPath    = paste0  (MULTIGWAS_HOME, "/opt/Rlibs")
libSources = paste0  (MULTIGWAS_HOME, "/install/repo/src/contrib")
libRepo    = sprintf ("file://%s/install/repo", MULTIGWAS_HOME)
if (!dir.exists (libPath))
	dir.create (libPath)

.libPaths (libPath)
message ("\n\nInstalling R libraries for MultiGWAS into: ", libPath)


multigwas_packages = c("rrBLUP", "parallel","config","dplyr", "stringi","qqman", "VennDiagram", "RColorBrewer","circlize",
					   "gplots", "rmarkdown", "kableExtra" ,"doParallel", "ldsep", "yaml", "BiocManager")
packagesInstalled    = rownames(installed.packages())
packagesNotInstalled = setdiff(multigwas_packages, rownames(installed.packages()))  

message ("Packages to install: ", paste (packagesNotInstalled));Sys.sleep (3)
install.packages (packagesNotInstalled, repos=libRepo)

if (!"multtest" %in% packagesInstalled)
	BiocManager::install("multtest", site_repository="file:///opt/miniCRAN")
if (!"GWASpoly" %in% packagesInstalled)
	install.packages(paste0(libSources,'/GWASpoly_1.3.tar.gz'), lib=libPath, repos=NULL, type="source") 
#------------------------------------------------------------------------------





.libPaths (libPath)

# Copy multiGWAS profile to multiGWAS home
x=file.copy ("multiGWAS_profile.sh", "..", overwrite=T)

message ("\n------------------------------------------\n")
message ("Close the terminal to finish the installation process")
message ("Then, open a new terminal and write: ")
message ("multiGWAS")
message ("\n------------------------------------------\n")




