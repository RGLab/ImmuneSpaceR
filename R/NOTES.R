# globalVariables to remove RCC NOTES due to data.table scoping
globalVariables(c("biosample_accession",
                "study_time_collected",
                "study_time_collected_unit",
                "arm_name",
                "response",
                "analyte",
                "analyte_name",
                "value_reported",
                "spot_number_reported",
                "cell_number_reported",
                "threshold_cycles",
                "entrez_gene_id",
                "mfi",
                "concentration_value",
                "population_cell_number",
                "population_name_reported",
                "x", "y", "err", #ggplot2 aes
                "ID", #qpHeatmap
                "stc", "stcu", "time_str", #standardize_time
                "virus_strain", "cohort",
                "participant_id"
                ))
