#!/usr/bin/env bash

# USAGE: bash magma.sh <working dir> <binary file prefix> <gene location file> <path to sample dir> <phenotype file> <gene set file> <output prefix> <directory of the IMML package> <number of iterations>

wdir=$1
bfile=$2
geneloc=$3
sample_dir=$4
phenofile=$5
genesetfile=$6
output=$7
package_dir=$8
n_iter=$9

cd ${wdir}

dirlist=(subsets annotates magma_gene magma_geneset)
for dir in  ${dirlist[@]}
do
    if [[ -d ${dir} ]]
    then
        echo "Folder ${dir} exists, making new folder..."
        rm -r ${dir}
        mkdir ${dir}
    else
        echo Making ${dir}
        mkdir ${dir}
    fi
done

echo Annotating SNPs to genes...  | sed G

magma --annotate window=2,0.5 --snp-loc ${bfile}.bim --gene-loc ${geneloc} --out ./annotates/annotates

echo Starting geneset analysis | sed G

for ((i=1; i<=${n_iter}; i++))
do
	###Removing training samples
	plink --silent --bfile ${bfile} --keep ${sample_dir}/samples_${i}.txt --make-bed --out ./subsets/${output}_${i}

	echo Running gene-level analysis... | sed G
  # Gene-level analysis for 1 bootstrapped subset
	magma --bfile ./subsets/${output}_${i} --gene-annot ./annotates/annotates.genes.annot --pheno file=${phenofile} --out ./magma_gene/${output}_${i}

	echo Running gene set-level analysis... | sed G
  # Gene-set level analysis for 1 bootstrapped subset
	magma --gene-results ./magma_gene/${output}_${i}.genes.raw --model alpha=0.3 --settings gene-info --set-annot ${genesetfile} --out ./magma_geneset/${output}_${i}

done

# Extract entire feature selection set
plink --silent --bfile ${bfile} --keep ${sample_dir}/samples_fs.txt --make-bed --out fs_${output}

# Gene analysis for entire feature selection set
magma --bfile ./fs_${output} --gene-annot ./annotates/annotates.genes.annot --pheno file=${phenofile} --out ./magma_gene/${output}

# Create the gene_annot.txt file
awk -v OFS=' ' '/^[^#]/ {for (i=2; i<=NF; i++) if ($i ~ /^[0-9]+:[0-9]+$/) print $1, $i}' ./annotates/annotates.genes.annot > ./annotates/gene_annot.txt

# Final GSEA to extract the leading edge SNPs
Rscript ${package_dir}/final_gsea.R ${output} ${n_iter}

# Extract redcued bed file for prunning
plink --silent -bfile ${bfile} --extract leadingEdge_snps_${output}.txt --keep ${sample_dir}/samples_fs.txt --make-bed --out tmp

# Look for LD SNPs
plink --bfile tmp --indep 50 5 2

# Make final bed file
plink --silent -bfile ${bfile} --extract plink.prune.in --make-bed --out genomics_${output}

