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

process_groups: 1-generate_group_data.R
	$(RS) 1-generate_group_data.R

process_diseases: 2-*.R 3-*.R
	$(RS) 2-download_influenza_data.R
	$(RS) 2-process_mers_data.R
	$(RS) 2-process_sars_data.R
	$(RS) 3-generate_disease_data.R

process_new: 6-*.R
	$(RS) 6-download_new_drug_data.R
	$(RS) 6-generate_new_group_data.R

analyse: 4-*.R
	$(RS) 4-generate_common_perturbagens_sars.R
	$(RS) 4-generate_common_perturbagens_sars2.R
	$(RS) 4-generate_sars_sars2_combined_list.R

visualise: f-*.R
	$(RS) f-generate_concordance_scatterplot-sars.R
	$(RS) f-generate_concordance_scatterplot-sars2.R
	$(RS) f-generate_threshold_histograms.R

all: process_drugs process_groups process_diseases analyse visualise
