# Modify this Snakemake call to your needs.
# If you run the analysis on a cluster, we recommend using a cluster configuration via --cluster-config

echo "####################################################################"
echo "#              Thank you for the interest in diffTF!               #"
echo "#    If you have questions or comments, feel free to contact us.   #"
echo "# We will be happy to answer any questions related to this project #"
echo "#    as well as questions related to the software implementation.  #"
echo "####################################################################"


echo "\nThis wrapper script executes Snakemake to start diffTF. Modify the Snakemake call to your needs.\n"

echo "#########################################################"
echo "#   NOTE THAT THIS ANALYSIS MAY TAKE A WHILE TO FINISH  #"
echo "# INCREASING THE NUMBER OF CORES SPEEDS UP THE ANALYSIS #"
echo "#########################################################"

# Real run, using 16 cores
snakemake --snakefile /home/prvaldes/scratch/diffTF/src/Snakefile --cores 16 \
--configfile /home/prvaldes/scratch/diffTF/input/config_V717I_cisBP.json \
--use-singularity --singularity-args "--bind /home/prvaldes/scratch/diffTF/ --bind /tmp"
