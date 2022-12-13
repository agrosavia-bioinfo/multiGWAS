#!/usr/bin/env bash
# Remove previous installations
"sed '/multiGWAS/d' -i ~/.bashrc"

# Start installation
MULTIGWAS_HOME=`dirname $PWD`
echo "MULTIGWAS_HOME=$MULTIGWAS_HOME"

echo "Creating multiGWAS profile..."
echo "\n#------------------- multiGWAS.R profile ---------------------" > $MULTIGWAS_HOME/multiGWAS_profile.sh
echo "export MULTIGWAS_HOME=$MULTIGWAS_HOME" >> $MULTIGWAS_HOME/multiGWAS_profile.sh
echo "MULTIGWAS_TOOLS=\$MULTIGWAS_HOME/opt/tools" >> $MULTIGWAS_HOME/multiGWAS_profile.sh
echo "MULTIGWAS_MAIN=\$MULTIGWAS_HOME/main" >> $MULTIGWAS_HOME/multiGWAS_profile.sh
echo "export PATH=\$PATH:\$MULTIGWAS_TOOLS:\$MULTIGWAS_MAIN" >> $MULTIGWAS_HOME/multiGWAS_profile.sh

# Write into .bashrc
echo "\n#------------------- multiGWAS.R tool profile ---------------------" >> ~/.bashrc
echo ". $MULTIGWAS_HOME/multiGWAS_profile.sh" >> ~/.bashrc

echo "\nMultiGWAS is ready to use, right after installed!\n"

echo "Decompressing zip file Rlibs-Ubuntu22-R41.zip into ../opt folder..."
unzip ../opt/Rlibs-Ubuntu22-R41.zip -d ../opt/
