library(readxl)
library(lubridate)
library(openxlsx)
library(tidyverse)
library(shiny)


# Clean population names (some contain PLATE or PANEL)
clean_pop <- function(s) {
  s <- str_to_upper(s)
  
  case_when(
    str_detect(s, "BARBADOS") ~ "AFRICAN ANCESTRY FROM BARBADOS IN THE CARIBBEAN",
    str_detect(s, "SOUTHWEST") ~ "AFRICAN ANCESTRY IN SOUTHWEST USA",
    str_detect(s, "BENGALI") ~ "BENGALI IN BANGLADESH",
    str_detect(s, "BRITISH") ~ "BRITISH FROM ENGLAND AND SCOTLAND, UK",
    str_detect(s, "XISHUANGBANNA") ~ "CHINESE DAI IN XISHUANGBANNA, CHINA",
    str_detect(s, "COLOMBIAN") ~ "COLOMBIAN IN MEDELLIN, COLOMBIA",
    str_detect(s, "ESAN") ~ "ESAN FROM NIGERIA",
    str_detect(s, "FINNISH") ~ "FINNISH IN FINLAND",
    str_detect(s, "GAMBIAN") ~ "GAMBIAN IN WESTERN DIVISION, THE GAMBIA",
    str_detect(s, "GUJARATI") ~ "GUJARATI INDIANS IN HOUSTON, TEXAS, USA",
    str_detect(s, "BEIJING") ~ "HAN CHINESE IN BEIJING, CHINA",
    str_detect(s, "CHINESE SOUTH") ~ "HAN CHINESE SOUTH, CHINA", 
    str_detect(s, "IBERIAN") ~ "IBERIAN POPULATIONS IN SPAIN",
    str_detect(s, "TELUGU") ~ "INDIAN TELUGU IN THE UK",
    str_detect(s, "JAPANESE") ~ "JAPANESE IN TOKYO, JAPAN",
    str_detect(s, "KINH") ~ "KINH IN HO CHI MINH CITY, VIETNAM",
    str_detect(s, "LUHYA") ~ "LUHYA IN WEBUYE, KENYA",
    str_detect(s, "MAASAI") ~ "MAASAI IN KINYAWA, KENYA",
    str_detect(s, "MENDE") ~ "MENDE IN SIERRA LEONE",
    str_detect(s, "MEXICAN") ~ "MEXICAN ANCESTRY IN LOS ANGELES, CALIFORNIA, USA",
    str_detect(s, "PERUVIAN") ~ "PERUVIAN IN LIMA, PERU",
    str_detect(s, "PUERTO RICAN") ~ "PUERTO RICAN IN PUERTO RICO",
    str_detect(s, "PUNJABI") ~ "PUNJABI IN LAHORE, PAKISTAN",
    str_detect(s, "SRI LANKAN") ~ "SRI LANKAN TAMIL IN THE UK",
    str_detect(s, "TOSCANI") ~ "TOSCANI IN ITALIA",
    str_detect(s, "YORUBA") ~ "YORUBA IN IBADAN, NIGERIA",
    str_detect(s, "DENVER") ~ "CHINESE IN METROPOLITAN DENVER CO USA",
    TRUE ~ "UNKNOWN_DESCRIPTION"
  )
}


# Create the header style
hs <- createStyle(fgFill = "#DCE6F2",
                  halign = "left", 
                  valign = "center", 
                  textDecoration = "Bold",
                  border = "TopBottomLeftRight")

# Create the title style for merged cells
title_style <- createStyle(fontName = "Calibri", 
                           fontSize = 14, 
                           halign = "center", 
                           valign = "center")

# Create cell style for lay_df cells
lay_style <- createStyle(halign = "left",
                         valign = "center",
                         wrapText = TRUE)

# create cell style for # Ordered column
order_style <- createStyle(halign = "center")

# custom function for setting row heights (used on Lay Summary column)
get_row_height <- function(s, col_width = 60, scale_factor = 20) {
  h <- str_length(s) / col_width * scale_factor
  h[is.na(h)] <- 15
  h
}

# custom function setting column widths
get_col_widths <- function(df, char_buff = 2) {
  apply(df, 2, function(x) max(nchar(as.character(x)) + char_buff, na.rm = TRUE))
}


# function for creating excel sheets
create_excel <- function(report_df, lay_df, title, from_date, to_date, filename) {
  LAY_COL_WIDTH <- 60
  DATE1 <- ymd(from_date)
  DATE2 <- ymd(to_date)
  
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "Research Intent")
  addWorksheet(wb, sheetName = "Lay Summaries")
  
  # Create the title cells
  writeData(wb, 
            sheet = 1, 
            x = paste0("NHGRI Sample Repository for Human Genetic Research Community Report\nCovering ", 
                       month(DATE1, label = TRUE, abbr = FALSE), " ", day(DATE1), ", ", year(DATE1),
                       " through ", 
                       month(DATE2, label = TRUE, abbr = FALSE), " ", day(DATE2), ", ", year(DATE2), 
                       "\n", title), 
            startCol = 1, 
            startRow = 1)
  setRowHeights(wb, sheet = 1, rows = 1, heights = 65)
  addStyle(wb, sheet = 1, style = title_style, cols = 1, rows = 1)
  mergeCells(wb, sheet = 1, cols = 1:7, rows = 1)
  setColWidths(wb, sheet = 1, cols = 1:4, ignoreMergedCells = TRUE, widths = get_col_widths(report_df, char_buff = 1)[1:4])
  setColWidths(wb, sheet = 1, cols = 7, ignoreMergedCells = TRUE, widths = get_col_widths(report_df, char_buff = 1)[7])
  setColWidths(wb, sheet = 1, cols = 5, widths = 8)
  setColWidths(wb, sheet = 1, cols = 6, widths = 10)
  
  # write the data to the spreadsheet
  writeData(wb, sheet = 1, x = report_df, headerStyle = hs, startRow = 3, borders = "all")
  writeData(wb, sheet = 2, x = lay_df, headerStyle = hs, startRow = 1, borders = "all")
  
  # Merge last and first columns into single Investigator column
  writeData(wb, sheet = 1, startRow = 3, startCol = 1, x = "Investigator")
  writeData(wb, sheet = 2, startRow = 1, startCol = 1, x = "Investigator")
  mergeCells(wb, sheet = 1, cols = 1:2, rows = 3)
  mergeCells(wb, sheet = 2, cols = 1:2, rows = 1)
  
  # apply styling to the lay summary sheet
  setColWidths(wb, sheet = 2, cols = 1:3,  widths = get_col_widths(lay_df)[1:3], ignoreMergedCells = TRUE)
  setColWidths(wb, sheet = 2, cols = 4, widths = 10)
  setColWidths(wb, sheet = 2, cols = 5, width = LAY_COL_WIDTH, ignoreMergedCells = TRUE)
  addStyle(wb, sheet = 2, style = lay_style, rows = 1:nrow(lay_df)+1, cols = 1:5, gridExpand = TRUE, stack = TRUE)
  
  setRowHeights(wb, sheet = 2, rows = 1:nrow(lay_df)+1, 
                heights = get_row_height(lay_df[["Lay Summary"]], 
                                         col_width = LAY_COL_WIDTH))
  
  # apply centering to # ordered column of report df
  addStyle(wb, sheet = 1, style = order_style, rows = 1:nrow(report_df) + 3, cols = 6, gridExpand = TRUE, stack = TRUE)

  saveWorkbook(wb, file = filename, overwrite = TRUE)
}


# clean and split the SSRS dump file dataframe
clean_and_split_ssrs <- function(ssrs_df) {
  # Extract JAPANESE IN TOKYO, JAPAN AND HAN CHINESE IN BEIJING, CHINA rows from df
  # these need to be duplicated, one for JAPAN and ONE for CHINA
  # find the indexes of these rows
  dup_idxs <- which(grepl("JAPANESE IN TOKYO, JAPAN AND HAN CHINESE IN BEIJING, CHINA", ssrs_df$diag_desc))
  if (!rlang::is_empty(dup_idxs)) {
    dup_rows <- ssrf_df[rep(dup_idxs, each = 2), ]
    dup_rows$Population <- rep(c("JAPANESE IN TOKYO, JAPAN", "HAN CHINESE IN BEIJING, CHINA"), length(dup_idxs))
    ssrs_df <- ssrs_df[-dup_idxs, ]
    ssrs_df <- bind_rows(ssrf_df, dup_rows)
  }

  # perform cleaning, summarization, and splitting on each population
  ssrs_df %>%
    filter(!is.na(ref) & !is.na(diag_desc)) %>%  # if the sample name and population are blank, remove
    mutate(population = clean_pop(diag_desc)) %>%
    group_by(population,
             customer_id,
             name,
             institution,
             country,
             product,
             r_intent_type) %>%
    summarize(number_ordered = sum(quantity_ordered),
              .groups = "drop") %>%
    select(population, name, institution, country, product, number_ordered, r_intent_type) %>%
    rename(Investigator = name,
           Institution = institution,
           Country = country,
           Product = product,
           `# Ordered` = number_ordered,
           `Research Intent` = r_intent_type,
           Population = population) %>%
    split(f = as.factor(.$Population))
}


# function for cleaning and unnesting the lay_df
clean_and_unnest_lay <- function(lay_df) {
  df <- lay_df %>%
    mutate(name = str_remove(customer, ",")) %>%
    select(name, population, lay_summary) %>%
    mutate(pops = str_split(population, pattern = "\n")) %>%
    select(name, lay_summary, pops) %>%
    unnest(pops) %>%
    mutate(population = clean_pop(pops)) %>%
    filter(population != "UNKNOWN_DESCRIPTION") %>% # assume these are from bad line breaks
    rename(Investigator = name,
           `Lay Summary` = lay_summary,
           Population = population)
  
  # same procedure for detecting and duplicating JAPANESE IN TOKYO, JAPAN AND HAN CHINESE IN BEIJING, CHINA rows
  dup_idxs <- which(grepl("JAPANESE IN TOKYO, JAPAN AND HAN CHINESE IN BEIJING, CHINA", df$pops))
  if (!rlang::is_empty(dup_idxs)) {
    dup_rows <- df[rep(dup_idxs, each = 2), ]
    dup_rows$Population <- rep(c("JAPANESE IN TOKYO, JAPAN", "HAN CHINESE IN BEIJING, CHINA"), length(dup_idxs))
    df <- df[-dup_idxs, ]
    df <- bind_rows(df, dup_rows)
  }
  df
}

# Join the lay df onto the report dfs
join_lay_onto_reports <- function(report_dfs, lay_unnested) {
  report_dfs %>%
    map(left_join, lay_unnested) %>%
    map(select, Investigator, Institution, Product, `Lay Summary`) %>%
    map(arrange, Investigator) %>%
    map(distinct)
}

# Split Investigator into first and last name and clean whitespace
split_investigator <- function(dfs) {
  dfs %>% 
    map(separate, Investigator, into = c("last", "first"), sep = "\\s", extra = "merge") %>% 
    map(mutate, first = str_trim(first, "both")) %>% 
    map(mutate, last = str_trim(last, "both"))
}

# check the report and lay dfs for missing values
check_processing <- function(reports, lays) {
  check_missing <- function(df) {
    missing <- sapply(df, function(x) {
      any(is.na(x))
    })
    names(missing[missing])
  }
  
  reports_missings <- lapply(reports, check_missing)
  reports_missings <- reports_missings[!sapply(reports_missings, rlang::is_empty)] %>% lapply(paste, collapse = ",")
  lay_missings <- lapply(lays, check_missing)
  lay_missings <- lay_missings[!sapply(lay_missings, rlang::is_empty)] %>% lapply(paste, collapse = ",")
  
  messages <- c("RESEARCH INTENT MESSAGES ----------")
  if (!rlang::is_empty(reports_missings)) {
    for (p in names(reports_missings)) {
      msg <- paste(reports_missings[[p]], "column(s) contain missing values in", p, "Research Intent table. Check data.")
      messages <- append(messages, msg)
    }
  }
  
  messages <- append(messages, "LAY SUMMARY MESSAGES ----------")
  if (!rlang::is_empty(lay_missings)) {
    for (p in names(lay_missings)) {
      msg <- paste(lay_missings[[p]], "column(s) contain missing values in", p, "Lay Summary table. Check data.")
      messages <- append(messages, msg)
    }
  }
  messages
}