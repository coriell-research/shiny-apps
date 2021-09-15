library(data.table)
library(openxlsx)


# Processing Functions ---------------------------------------------------------
# Function for reading input files
readInput <- function(fpath, skip_rows) {
  dt <- fread(fpath, encoding = "UTF-8", skip = skip_rows)
  return(dt)
}

# Function for renaming the populations
cleanPopulation <- function(s) {
  s <- toupper(s)
  clean_pop <- fcase(
    grepl("BARBADOS", s), "AFRICAN ANCESTRY FROM BARBADOS IN THE CARIBBEAN",
    grepl("SOUTHWEST", s), "AFRICAN ANCESTRY IN SOUTHWEST USA",
    grepl("BENGALI", s), "BENGALI IN BANGLADESH",
    grepl("BRITISH", s), "BRITISH FROM ENGLAND AND SCOTLAND, UK",
    grepl("XISHUANGBANNA", s), "CHINESE DAI IN XISHUANGBANNA, CHINA",
    grepl("COLOMBIAN", s), "COLOMBIAN IN MEDELLIN, COLOMBIA",
    grepl("ESAN", s), "ESAN FROM NIGERIA",
    grepl("FINNISH", s), "FINNISH IN FINLAND",
    grepl("GAMBIAN", s), "GAMBIAN IN WESTERN DIVISION, THE GAMBIA",
    grepl("GUJARATI", s), "GUJARATI INDIANS IN HOUSTON, TEXAS, USA",
    grepl("BEIJING", s), "HAN CHINESE IN BEIJING, CHINA",
    grepl("CHINESE SOUTH", s), "HAN CHINESE SOUTH, CHINA", 
    grepl("IBERIAN", s), "IBERIAN POPULATIONS IN SPAIN",
    grepl("TELUGU", s), "INDIAN TELUGU IN THE UK",
    grepl("JAPANESE", s), "JAPANESE IN TOKYO, JAPAN",
    grepl("KINH", s), "KINH IN HO CHI MINH CITY, VIETNAM",
    grepl("LUHYA", s), "LUHYA IN WEBUYE, KENYA",
    grepl("MAASAI", s), "MAASAI IN KINYAWA, KENYA",
    grepl("MENDE", s), "MENDE IN SIERRA LEONE",
    grepl("MEXICAN", s), "MEXICAN ANCESTRY IN LOS ANGELES, CALIFORNIA, USA",
    grepl("PERUVIAN", s), "PERUVIAN IN LIMA, PERU",
    grepl("PUERTO RICAN", s), "PUERTO RICAN IN PUERTO RICO",
    grepl("PUNJABI", s), "PUNJABI IN LAHORE, PAKISTAN",
    grepl("SRI LANKAN", s), "SRI LANKAN TAMIL IN THE UK",
    grepl("TOSCANI", s), "TOSCANI IN ITALIA",
    grepl("YORUBA", s), "YORUBA IN IBADAN, NIGERIA",
    grepl("DENVER", s), "CHINESE IN METROPOLITAN DENVER CO USA",
    default = "UNKNOWN_POPULATION"
  )
  return(clean_pop)
}

# Function for detecting "JAPANESE IN TOKYO, JAPAN AND HAN CHINESE IN BEIJING, CHINA" and then duplicating the rows
dupRows <- function(dt, pop_col) {
  dt <- copy(dt)
  dup_idxs <- which(grepl("JAPANESE IN TOKYO, JAPAN AND HAN CHINESE IN BEIJING, CHINA", dt[[pop_col]]))
  
  if (!rlang::is_empty(dup_idxs)) {
    dup_rows <- dt[rep(dup_idxs, each = 2), ]
    dup_rows$Population <- rep(c("JAPANESE IN TOKYO, JAPAN", "HAN CHINESE IN BEIJING, CHINA"), length(dup_idxs))
    dt <- dt[-dup_idxs, ]
    dt <- rbind(dt, dup_rows)
  }
  return(dt)
}

# Function for cleaning, grouping, and splitting SSRS dt
cleanGroupSplitSSRS <- function(dt) {
  # duplicate rows from combined panel
  dt <- dupRows(dt, pop_col = "DiagDesc")
  
  # check for and remove any all NA rows
  dt <- dt[!apply(dt, 1, function(x) all(is.na(x))), ]
  
  # clean populations and group and summarize orders
  dt <- dt[Order_Type != "Replacement"][, 
           Population := cleanPopulation(DiagDesc)][,
           .("# Ordered" = sum(Quantity)), 
           by = c("Population", "Name", "Institution", "Country", "Product", "RIntentType")]
  
  # rename columns to match desired output
  setnames(
    dt, 
    old = c("# Ordered", "Population", "Name", "Institution", "Country", "Product", "RIntentType"),
    new = c("# Ordered", "Population", "Investigator", "Institution", "Country", 
            "Product", "Research Intent"))
  
  dt <- dt[, .(Investigator, Institution, Country, Product, `# Ordered`, `Research Intent`, Population)]
  
  # split the grouped data.table into a list of dts by population
  ssrs_dts <- split(dt, by = "Population")
  
  return(ssrs_dts)
}

# Function for cleaning and unnesting the lay summary file
cleanAndUnestLay <- function(dt) {
  dt <- copy(dt)
  
  # remove commas from names and unnnest based on population
  dt[, Customer := trimws(Customer)]
  dt <- dt[, Customer := gsub(",", "", Customer), ][, 
            .(Customer, Population, Lay_Summary1)][, 
            .(Population = strsplit(Population, "\r\n|\n", perl = TRUE), Customer, Lay_Summary1)][, 
                                                                                                                                                                                 .(Population = as.character(unlist(Population))), by = c("Customer", "Lay_Summary1")]
  
  # check for duplicates in combined panel
  dt <- dupRows(dt, pop_col = "Population")
  
  dt[Population != "UNKNOWN_POPULATION",
     .(Population = cleanPopulation(Population),
      Investigator = Customer,
      `Lay Summary` = Lay_Summary1)]
}

# Function for joining lay summaries onto SSRS report data.tables
joinLayOntoSSRS <- function(ssrs_list, lay_dt) {
  joined <- lapply(ssrs_list, function(ssrs_dt) {
    merge(x = ssrs_dt, 
          y = lay_dt, 
          by.x = c("Investigator", "Population"), 
          by.y = c("Investigator", "Population"),
          all.x = TRUE,
          all.y = FALSE)
  })
  
  # Select required columns and sort
  selected <- lapply(
    joined, function(dt) {
      dt[, .(Investigator, Institution, Product, `Lay Summary`)][order(Investigator)]
    }
  )
  
  # Remove any possible duplicates
  uniqs <- lapply(selected, function(dt) {unique(dt)})
  
  return(uniqs)
}

# Function to split the `Investigator` column into first and last names
splitInvestigator <- function(dt) {
  dt <- copy(dt)
  dt <- tidyr::separate(dt, col = Investigator, into = c("last", "first"), sep = "\\s", extra = "merge")
  dt[, `:=`(first = trimws(first), last = trimws(last))]
}

# Function for checking for missing values in final reports
checkProcessing <- function(reports, lays) {
  check_missing <- function(df) {
    missing <- sapply(df, function(x) {
      any(is.na(x))
    })
    names(missing[missing])
  }
  
  reports_missings <- lapply(reports, check_missing)
  reports_missings <- lapply(reports_missings[!sapply(reports_missings, rlang::is_empty)], paste, collapse = ",")
  lay_missings <- lapply(lays, check_missing)
  lay_missings <- lapply(lay_missings[!sapply(lay_missings, rlang::is_empty)], paste, collapse = ",")
  
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

# Excel Styling Functions ------------------------------------------------------
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
  h <- nchar(s) / col_width * scale_factor
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
  DATE1 <- lubridate::ymd(from_date)
  DATE2 <- lubridate::ymd(to_date)
  
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "Research Intent")
  addWorksheet(wb, sheetName = "Lay Summaries")
  
  # Create the title cells
  writeData(wb, 
            sheet = 1, 
            x = paste0("NHGRI Sample Repository for Human Genetic Research Community Report\nCovering ", 
                       lubridate::month(DATE1, label = TRUE, abbr = FALSE), " ", lubridate::day(DATE1), ", ", lubridate::year(DATE1),
                       " through ", 
                       lubridate::month(DATE2, label = TRUE, abbr = FALSE), " ", lubridate::day(DATE2), ", ", lubridate::year(DATE2), 
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
  setColWidths(wb, sheet = 2, cols = 5, widths = LAY_COL_WIDTH, ignoreMergedCells = TRUE)
  addStyle(wb, sheet = 2, style = lay_style, rows = 1:nrow(lay_df)+1, cols = 1:5, gridExpand = TRUE, stack = TRUE)
  
  setRowHeights(wb, sheet = 2, rows = 1:nrow(lay_df)+1, 
                heights = get_row_height(lay_df[["Lay Summary"]], 
                                         col_width = LAY_COL_WIDTH))
  
  # apply centering to # ordered column of report df
  addStyle(wb, sheet = 1, style = order_style, rows = 1:nrow(report_df) + 3, cols = 6, gridExpand = TRUE, stack = TRUE)
  
  saveWorkbook(wb, file = filename, overwrite = TRUE)
}
