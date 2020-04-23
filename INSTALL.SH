#!/bin/bash
#
# Create multiGWAS profile according to current directory
# 
#
echo "Creating multiGWAS profile..."

echo ""
echo "#------------------- multiGWAS HOME ---------------------" >  multiGWAS-profile.sh
cmm1="export MULTIGWAS_HOME=$PWD" 
echo $cmm1 >> multiGWAS-profile.sh
eval $cmm1
cmm2="MULTIGWAS_TOOLS=$MULTIGWAS_HOME/tools" 
echo $cmm2 >> multiGWAS-profile.sh
eval $cmm2

cmm3="MULTIGWAS_SOURCES=$MULTIGWAS_HOME/sources"
echo $cmm3 >> multiGWAS-profile.sh
eval $cmm3

cmm4="export PATH=$PATH:$MULTIGWAS_TOOLS:$MULTIGWAS_SOURCES"
echo $cmm4 >> multiGWAS-profile.sh
eval $cmm4

echo ""
echo "MultiGWAS is ready to use, right after installed!"
echo ""
echo ""
echo "#------------------- multiGWAS.R tool profile ---------------------" >> $HOME/.profile
echo "source $PWD/multiGWAS-profile.sh" >> $HOME/.profile

