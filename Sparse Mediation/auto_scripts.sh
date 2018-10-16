### How I automatically generate x number of scripts to qsub for the sparse_qtl_scan.R
#       Don't need the #PBS lines.

#PBS -l nodes=1:ppn=4
#PBS -q batch
#PBS -l walltime=72:00:00

for i in {1..50}
do
   echo '#PBS -l nodes=1:ppn=4' >> "run_4_$i.sh"
   echo '#PBS -q batch' >> "run_4_$i.sh"
   echo '#PBS -l walltime=72:00:00' >> "run_4_$i.sh"
   echo 'module load R/3.5.1' >> "run_4_$i.sh"
   echo "Rscript sparse_qtl_scan.R 4 972 $i 98" >> "run_4_$i.sh"
   qsub "run_4_$i.sh"
done
