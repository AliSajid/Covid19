RS=Rscript
DATA_DIR=dataset/
RES_DIR=results/
FIG_DIR=figures/

.PHONY: complete
	
complete: clean all

clean:
	find $(DATA_DIR) -name "*.tsv" -delete
	find $(RES_DIR) -name "*.tsv" -delete
	find $(RES_DIR) -name "*.csv" -delete
	find $(FIG_DIR) -name "*.png" -delete
	find $(FIG_DIR) -name "*.jpg" -delete
	find $(FIG_DIR) -name "*.pdf" -delete

process_drugs: process_drug*.R
	Rscript process_drug_data.R
	Rscript process_drug_groups.R

process_disease: process_disease_data.R
	Rscript process_disease_data.R

analyse: analyse_data.R
	Rscript analyse_data.R

visualise: figure*.R
	Rscript figure3.R
	Rscript figure3v2.R
	Rscript figure3v3.R
	Rscript generate_histograms.R	

all: process_drugs process_disease analyse visualise
