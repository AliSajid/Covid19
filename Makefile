RS=Rscript
DATA_DIR=data/
RES_DIR=results/
FIG_DIR=figures/

.PHONY: complete

complete: clean all

clean:
	rm -rfv $(DATA_DIR)
	rm -rfv $(RES_DIR)
	rm -rfv $(FIG_DIR)
	$(RS) 0-setup.R

process_drugs: 1-download_drug_data.R
	$(RS) 1-download_drug_data.R

process_groups: 2-generate_group_data.R
	$(RS) 2-generate_group_data.R

process_diseases: 3-generate_disease_data.R
	$(RS) 3-generate_disease_data.R

analyse: 4-generate_common_perturbagens.R
	$(RS) 4-generate_common_perturbagens_sars.R
	$(RS) 5-generate_common_perturbagens_sars2.R

visualise:
	echo "Visualizations To Be Generated"

all: process_drugs process_groups process_disease analyse visualise
