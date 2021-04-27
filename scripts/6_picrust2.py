#Copyright (c) 2021 Insitute of Evolutionary Ecology and Conservation Genomics, Ulm University, Germany. See the LICENSE file at the top-level directory of this distribution (https://github.com/gfackelmann/human-encroachment-into-wildlife-gut-microbiomes/blob/main/LICENSE). 
#Author: Gloria Fackelmann (gloria.fackelmann@uni-ulm.de; https://www.researchgate.net/profile/Gloria-Fackelmann)

#Python 3.6.11
conda create -n picrust2 -c bioconda -c conda-forge picrust2=2.1.4_b
conda activate picrust2

###############################################################################################
#first, generate unstrat table on whole dataset at once:
###############################################################################################

biom head -i otu_biom_whole_dataset.biom

less rep-seqs.fna

mkdir picrust2_out_pipeline_whole_dataset

cd picrust2_out_pipeline_whole_dataset/

place_seqs.py -s ../rep-seqs.fna -o out.tre -p 1 --intermediate intermediate/place_seqs

hsp.py -i 16S -t out.tre -o marker_predicted_and_nsti.tsv.gz -p 1 -n

hsp.py -i EC -t out.tre -o EC_predicted.tsv.gz -p 1

zless -S marker_predicted_and_nsti.tsv.gz

zless -S EC_predicted.tsv.gz

metagenome_pipeline.py -i ../otu_biom_whole_dataset.biom -m marker_predicted_and_nsti.tsv.gz -f EC_predicted.tsv.gz -o EC_metagenome_out

zless -S EC_metagenome_out/pred_metagenome_unstrat.tsv.gz

pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -o pathways_out -p 1

add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o pathways_out/path_abun_unstrat_descrip.tsv.gz

#unzip 'path_abun_unstrat_descrip.tsv.gz' file and replace spaces with underscore.

###############################################################################################
#second, generate strat table by dividing dataset into 4 by samples (not by ASVs):
###############################################################################################

#part1:
biom summarize-table -i otu_biom_whole_dataset_part1.biom
#Num samples: 96
#Num observations: 2.662
#Total count: 2.634.703

less rep-seqs.fna

mkdir picrust2_out_pipeline_whole_part1

cd picrust2_out_pipeline_whole_part1/

place_seqs.py -s ../rep-seqs.fna -o out.tre -p 1 --intermediate intermediate/place_seqs

hsp.py -i 16S -t out.tre -o marker_predicted_and_nsti.tsv.gz -p 1 -n

hsp.py -i EC -t out.tre -o EC_predicted.tsv.gz -p 1

zless -S marker_predicted_and_nsti.tsv.gz

zless -S EC_predicted.tsv.gz

metagenome_pipeline.py -i ../otu_biom_whole_dataset_part1.biom -m marker_predicted_and_nsti.tsv.gz -f EC_predicted.tsv.gz -o EC_metagenome_out --strat_out

zless -S EC_metagenome_out/pred_metagenome_unstrat.tsv.gz

zless -S EC_metagenome_out/pred_metagenome_strat.tsv.gz

pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_strat.tsv.gz -o pathways_out -p 1

add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o pathways_out/path_abun_unstrat_descrip.tsv.gz

#part2:
biom summarize-table -i otu_biom_whole_dataset_part2.biom
#Num samples: 96
#Num observations: 2.777
#Total count: 2.285.363

mkdir picrust2_out_pipeline_whole_part2

cd picrust2_out_pipeline_whole_part2/

place_seqs.py -s ../rep-seqs.fna -o out.tre -p 1 --intermediate intermediate/place_seqs

hsp.py -i 16S -t out.tre -o marker_predicted_and_nsti.tsv.gz -p 1 -n

hsp.py -i EC -t out.tre -o EC_predicted.tsv.gz -p 1

zless -S marker_predicted_and_nsti.tsv.gz

zless -S EC_predicted.tsv.gz

metagenome_pipeline.py -i ../otu_biom_whole_dataset_part2.biom -m marker_predicted_and_nsti.tsv.gz -f EC_predicted.tsv.gz -o EC_metagenome_out --strat_out

zless -S EC_metagenome_out/pred_metagenome_unstrat.tsv.gz

zless -S EC_metagenome_out/pred_metagenome_strat.tsv.gz

pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_strat.tsv.gz -o pathways_out -p 1

add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o pathways_out/path_abun_unstrat_descrip.tsv.gz

#part3:
biom summarize-table -i otu_biom_whole_dataset_part3.biom
#Num samples: 96
#Num observations: 2.866
#Total count: 2.298.655

mkdir picrust2_out_pipeline_whole_part3

cd picrust2_out_pipeline_whole_part3/

place_seqs.py -s ../rep-seqs.fna -o out.tre -p 1 --intermediate intermediate/place_seqs

hsp.py -i 16S -t out.tre -o marker_predicted_and_nsti.tsv.gz -p 1 -n

hsp.py -i EC -t out.tre -o EC_predicted.tsv.gz -p 1

zless -S marker_predicted_and_nsti.tsv.gz

zless -S EC_predicted.tsv.gz

metagenome_pipeline.py -i ../otu_biom_whole_dataset_part3.biom -m marker_predicted_and_nsti.tsv.gz -f EC_predicted.tsv.gz -o EC_metagenome_out --strat_out

zless -S EC_metagenome_out/pred_metagenome_unstrat.tsv.gz

zless -S EC_metagenome_out/pred_metagenome_strat.tsv.gz

pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_strat.tsv.gz -o pathways_out -p 1

add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o pathways_out/path_abun_unstrat_descrip.tsv.gz

#part4:
biom summarize-table -i otu_biom_whole_dataset_part4.biom
#Num samples: 96
#Num observations: 2.990
#Total count: 2.677.770

mkdir picrust2_out_pipeline_whole_part4

cd picrust2_out_pipeline_whole_part4/

place_seqs.py -s ../rep-seqs.fna -o out.tre -p 1 --intermediate intermediate/place_seqs

hsp.py -i 16S -t out.tre -o marker_predicted_and_nsti.tsv.gz -p 1 -n

hsp.py -i EC -t out.tre -o EC_predicted.tsv.gz -p 1

zless -S marker_predicted_and_nsti.tsv.gz

zless -S EC_predicted.tsv.gz

metagenome_pipeline.py -i ../otu_biom_whole_dataset_part4.biom -m marker_predicted_and_nsti.tsv.gz -f EC_predicted.tsv.gz -o EC_metagenome_out --strat_out

zless -S EC_metagenome_out/pred_metagenome_unstrat.tsv.gz

zless -S EC_metagenome_out/pred_metagenome_strat.tsv.gz

pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_strat.tsv.gz -o pathways_out -p 1

add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o pathways_out/path_abun_unstrat_descrip.tsv.gz

#Unzip all the 'pathways_out/path_abun_strat.tsv.gz' files.