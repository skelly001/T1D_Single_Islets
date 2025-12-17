# ==============================================================================
# Donor Information Summary Table for Publication
# ==============================================================================
# Script: 10b_1_table_1_donor_information.R
# Description: Creates publication Table 1 summarizing donor demographic,
#              clinical, and HLA typing information. Converts HTML table to PNG
#              image using webshot2 for manuscript inclusion.
#
# Input: - data/donor_table.xlsx
# Output: - output/RD10-misc/RD10b_1-table_1_donor_information.png
# ==============================================================================

library(readxl)
library(tidyverse)
library(officer)
library(flextable)
library(webshot2)


# Import -----------------------------------------------------------------

# Load donor metadata
df <- read_excel("data/donor_table.xlsx", sheet = "HLA")


# Process ----------------------------------------------------------------

# Fill missing values from repeated donor rows
df_complete <- df %>%
	mutate(`nPOD ID` = as.character(`nPOD ID`)) %>%
	fill(`nPOD ID`, .direction = "down") %>%
	group_by(`nPOD ID`) %>%
	fill(
		`Case Type`,
		Race,
		`AAb detected`,
		Age,
		Sex,
		HbA1C,
		`Insulitis Status`,
		BMI,
		.direction = "down"
	)

# Reorder columns for publication format
table_data <- df_complete %>%
	relocate(`Insulitis Status`, `AAb detected`, .after = `Case Type`)


## Create Flextable ----
ft <- flextable(table_data) %>%

	# Merge cells
	merge_v(j = c("nPOD ID", "Case Type", "AAb detected", "Age", "Insulitis Status", "BMI")) %>%
	valign(
		j = c("nPOD ID", "Case Type", "AAb detected", "Age", "Insulitis Status", "BMI"),
		valign = "center"
	) %>%
	{
		reduce(seq(1, 12, 2), ~ merge_at(.x, .y:(.y + 1), "Race"), .init = .)
	} %>%
	{
		reduce(seq(1, 12, 2), ~ merge_at(.x, .y:(.y + 1), "Sex"), .init = .)
	} %>%
	{
		reduce(seq(1, 12, 2), ~ merge_at(.x, .y:(.y + 1), "HbA1C"), .init = .)
	} %>%

	# Center align all headers
	align(align = "center", part = "header") %>%
	align(align = "center", part = "body") %>%

	# Add column spanner for HLA types
	add_header_row(
		values = c("", "", "", "", "", "", "", "", "", "HLA"),
		colwidths = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 5)
	) %>%

	autofit(part = "header") %>%
	bold(part = "header") %>%

	# Highlight the HLA group header with background color
	bg(i = 1:2, j = 1:14, bg = "#E8E8E8", part = "header") %>%

	# Remove border between HLA title and column names
	hline(i = 1, j = 1:9, border = fp_border(width = 0), part = "header") %>%

	# Add outer borders
	border_outer(border = fp_border(color = "grey70", width = 2), part = "all") %>%

	# body grid
	border_inner_h(border = fp_border(color = "gray70", width = 1), part = "body") %>%
	border_inner_v(border = fp_border(color = "gray70", width = 1), part = "body")


# Save as HTML (to preserve background)
save_as_html(ft, path = "output/RD10-misc/donor_table.html")

# Use webshot2 to convert HTML to PNG
webshot(
	"output/RD10-misc/donor_table.html",
	"output/RD10-misc/RD10b_1-table_1_donor_information.png",
	vwidth = 1100,
	vheight = 800,
	zoom = 2.5
)

# Cleanup temporary file
file.remove("output/RD10-misc/donor_table.html")
