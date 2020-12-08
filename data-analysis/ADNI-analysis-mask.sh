#!/bin/bash

dirmask="/Users/bingkaiwang/Dropbox/research/graphical-model/hierarchical-lasso/compositional-hierarchical-tree-regression/data-analysis/mask/"

dirbrainmap="/Users/bingkaiwang/Dropbox/research/graphical-model/hierarchical-lasso/compositional-hierarchical-tree-regression/data-analysis/brainmap/"

cd $dirbrainmap
pwd

#out_alpha_allgroup="alpha_allgroup.nii.gz"
#fslmaths $dirmask"mask1.nii.gz" -mul 0 $out_alpha_allgroup
#while read -r rowname alpha index
#do
#	echo $index
#	fslmaths $dirmask"mask"$index".nii.gz" -mul $alpha maskmultmp.nii.gz
#	fslmaths $out_alpha_allgroup -add maskmultmp.nii.gz $out_alpha_allgroup
#done < ../alpha-allgroup.txt
#rm maskmultmp.nii.gz

#out_alpha_disease="alpha_disease.nii.gz"
#fslmaths $dirmask"mask1.nii.gz" -mul 0 $out_alpha_disease
#while read -r rowname alpha index
#do
#	echo $index
#	fslmaths $dirmask"mask"$index".nii.gz" -mul $alpha maskmultmp.nii.gz
#	fslmaths $out_alpha_disease -add maskmultmp.nii.gz $out_alpha_disease
#done < ../alpha-disease.txt
#rm maskmultmp.nii.gz

beta_group=(1 2 3 4)
for i in "${beta_group[@]}"
do
	input="beta"$i".txt"
	output="beta"$i".nii.gz"
	fslmaths $dirmask"mask1.nii.gz" -mul 0 $output
	while read -r rowname beta index
	do
		echo $index
		fslmaths $dirmask"mask"$index".nii.gz" -mul $beta maskmultmp.nii.gz
		fslmaths $output -add maskmultmp.nii.gz $output
	done < ../$input
	rm maskmultmp.nii.gz
done

