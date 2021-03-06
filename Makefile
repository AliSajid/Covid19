RS=Rscript
DATA_DIR=data/
RES_DIR=results/
FIG_DIR=figures/
MAP_DIR=maps/

.PHONY: complete

complete: clean all

clean:
	rm -rfv $(DATA_DIR)
	rm -rfv $(RES_DIR)
	rm -rfv $(FIG_DIR)
	rm -rfv $(MAP_DIR)
	$(RS) 0-setup.R

make_lists: 1-explore_drug_sig_lists.R
	$(RS) 1-explore_drug_sig_lists.R

process_drugs: 2-download_drug_data.R
	$(RS) 2-download_drug_data.R

process_groups: 2-generate_group_data.R
	$(RS) 2-generate_group_data.R

process_diseases: 3-*.R 4-*.R
	$(RS) 3-download_influenza_data.R
	$(RS) 3-process_mers_data.R
	$(RS) 3-process_sars_data.R
	$(RS) 4-generate_disease_data.R

analyse: 5-*.R 6-*.R
	$(RS) 5-generate_common_perturbagens_sars.R
	$(RS) 5-generate_common_perturbagens_sars2.R
	$(RS) 5-generate_common_perturbagens_covidc.R
	$(RS) 5-generate_common_perturbagens_covidm.R
	$(RS) 5-generate_common_perturbagens_ace2.R
	$(RS) 6-generate_sars_sars2_combined_list.R
	$(RS) 6-output_summarized_dataset_sars2.R
	$(RS) 6-output_signature_normality.R
	$(RS) 6-output_summarized_dataset_covidc.R
	$(RS) 6-output_summarized_dataset_covidm.R
	$(RS) 6-output_summarized_dataset_ace2.R
	$(RS) 6-output_threshold_tables.R
	$(RS) 6-output_disease_threshold_tables.R

visualise: f-*.R
	$(RS) f-generate_concordance_scatterplot-sars.R
	$(RS) f-generate_concordance_scatterplot-sars2.R
	$(RS) f-generate_threshold_histograms.R
	$(RS) f-generate_correlation_plot_covic_covidm.R
	$(RS) f-generate_mean_sd_plot.R
	$(RS) f-generate_upset_figure.R

all: make_lists process_drugs process_groups process_diseases analyse visualise
