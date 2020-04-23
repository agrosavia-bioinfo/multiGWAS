
# Create multiGWAS profile according to current directory
echo "Creating multiGWAS profile..."

echo "export MULTIGWAS_HOME=$PWD" > multiGWAS-profile.sh
echo 'MULTIGWAS_TOOLS=$MULTIGWAS_HOME/tools' >> multiGWAS-profile.sh
echo 'MULTIGWAS_SOURCES=$MULTIGWAS_HOME/sources' >> multiGWAS-profile.sh
echo 'export PATH=$PATH:$MULTIGWAS_TOOLS:$MULTIGWAS_SOURCES'>> multiGWAS-profile.sh

# Update user profile to execute multiGWAS-profile at startup
echo "Updating user profile..."

echo "\n#------------------- multiGWAS.R tool profi profile ---------------------" >> $HOME/.profile
echo "source $PWD/multiGWAS-profile.sh" >> $HOME/.profile

