#PBS -l nodes=1:ppn=1
#PBS -q batch
#PBS -l walltime=72:00:00

for i in {101..218}
do
   echo "#PBS -l nodes=1:ppn=1
#PBS -q batch
#PBS -l walltime=24:00:00
module load R/3.5.1

Rscript sparse_matrix_mediation_run.R 1500_2000 $i 50000 dataset.islet.proteins dataset.islet.mrna protein_id gene_id" >> "sparse_qtl_protein_mrna_mediation_1500_2000_$i.sh"
qsub "sparse_qtl_protein_mrna_mediation_1500_2000_$i.sh"
done
