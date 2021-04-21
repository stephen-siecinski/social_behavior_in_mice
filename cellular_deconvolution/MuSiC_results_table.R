library(gt)
library(webshot)

options(digits = 0, scipen = 999)

# Making a table of the cell counts for the Allen Reference dataset ####
allen_cell_table <- allen_metadata %>%
  filter(!is.na(subclass_label) & subclass_label != "") %>%
  group_by(external_donor_name_label, donor_sex_label) %>%
  count(subclass_label)

annotation_df <- filter(distinct(allen_metadata, 
                                 subclass_label, .keep_all = TRUE)) %>%
  filter(!is.na(subclass_label)) %>%
  dplyr::select(class_label, subclass_label)


allen_cell_table_males <- allen_cell_table %>%
  filter(donor_sex_label == "M") %>%
  group_by(subclass_label) %>%
  summarise(total_m = sum(n), 
            mean_m = mean(n),
            sd_m = sd(n), 
            min_m = min(n), 
            max_m = max(n),
            `n mice m` = length(n != 0))

allen_cell_table_females <- allen_cell_table %>%
  filter(donor_sex_label == "F") %>%
  group_by(subclass_label) %>%
  summarise(total = sum(n), 
            mean = mean(n),
            sd = sd(n), 
            min = min(n), 
            max = max(n),
            `n mice` = length(n != 0))

allen_cell_table_2 <- full_join(allen_cell_table_females, allen_cell_table_males)
allen_cell_table_2 <- allen_cell_table_2 %>% replace(is.na(.), 0) %>%
  arrange(desc(total)) %>% 
  left_join(annotation_df, by = "subclass_label") %>%
  dplyr::select(class_label, subclass_label, everything())



library(gt)

allen_cell_table_pretty <- allen_cell_table_2 %>% 
  gt(groupname_col = "class_label", auto_align = TRUE) %>%
  tab_style(style = list(cell_text(weight = "bold")),
            locations = cells_column_labels(everything())) %>%
  opt_table_lines(extent = "default") %>%
  tab_header(
    title = md("Sub-class Cell Counts by Sex")) %>%
  tab_spanner(label = "Females",
              columns = 3:8) %>%
  tab_spanner(label = "Males",
              columns = 9:14) %>%
  cols_label(subclass_label = "cell classification",
             total_m = "total",
             mean_m = "mean",
             sd_m = "sd",
             min_m = "min",
             max_m = "max",
             'n mice m' = "n mice") %>%
  tab_options(
    column_labels.border.top.color = "black",
    column_labels.border.top.width = px(0),
    column_labels.border.bottom.color = "black",
    table_body.hlines.color = "white",
    table.border.bottom.color = "white",
    table.border.bottom.width = px(0),
    data_row.padding = px(0),
    column_labels.font.size = 12,
    table.font.size = 10,
    column_labels.font.weight = "bold",
    row_group.background.color = "black",
    row.striping.include_table_body = TRUE,
    row.striping.background_color = "#eaeaea"
    #table.font.size = 2
    #container.width = px(1200),
    #row_group.padding = px(-3)
    
    #container.height = pct(100)
  )


gtsave(allen_cell_table_pretty, 
       filename = "20210209_allen_cell_counts_table.png", 
       path = "/Users/sks36/Desktop/R/mouse_project/MuSiC/methods/outputs")


# Making a table of the cell cell estimate results ####

results_table <- results_combined %>%
  dplyr::select(sample_id, strain, brain_region, sex, sociability, 2:14)

results_table <- results_table[1:20,] %>% 
  gt(auto_align = TRUE) %>%
  #cols_width(vars("Global AF") ~ pct(10),
  #           vars("Association") ~ pct(60)) %>%
  tab_style(style = list(cell_text(weight = "bold")),
            locations = cells_column_labels(everything())) %>%
  #tab_style(style = list(cell_fill(color = "#b1fffb")),
  #          locations = cells_body(rows = c(5,16,34,38))) %>%
  opt_table_lines(extent = "default") %>%
  tab_header(
    title = md("MuSiC Results Table")) %>%
  #tab_spanner(label = "Females",
  #            columns = 3:8) %>%
  # tab_spanner(label = "Males",
  #             columns = 9:14) %>%
  # cols_label(subclass_label = "cell classification",
  #            total_m = "total",
  #            mean_m = "mean",
  #            sd_m = "sd",
  #            min_m = "min",
  #            max_m = "max",
  #            'n mice m' = "n mice") %>%
  tab_options(
    column_labels.border.top.color = "black",
    column_labels.border.top.width = px(0),
    column_labels.border.bottom.color = "black",
    table_body.hlines.color = "white",
    table.border.bottom.color = "white",
    table.border.bottom.width = px(0),
    data_row.padding = px(0),
    column_labels.font.size = 12,
    table.font.size = 10,
    column_labels.font.weight = "bold",
    row_group.background.color = "black",
    row.striping.include_table_body = TRUE,
    row.striping.background_color = "#eaeaea"
    #table.font.size = 2
    #container.width = px(1200),
    #row_group.padding = px(-3)
    
    #container.height = pct(100)
  )


gtsave(results_table, 
       filename = "20210219_MuSiC_short_results_table.png", 
       path = outputs)

