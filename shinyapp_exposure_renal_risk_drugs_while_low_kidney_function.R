for (p in c("shiny", "dplyr", "RPostgreSQL", "ggplot2", "stringr", "tidyr", "plotly", "lubridate"))
	library(p, character.only = TRUE)

con <- DBI::dbConnect(RPostgreSQL::PostgreSQL(), 
					  db = Sys.getenv("db"),
					  host = Sys.getenv("host"),
					  port = Sys.getenv("port"),
					  user = Sys.getenv("user"),
					  password = Sys.getenv("pass"))
on.exit(dbDisconnect(con))

drugs_of_interest <- c(Digoxin = 1326303,
					   Gabapentin = 797399,
                       Lithium = 19124477,
                       Metformin = 1503297,
                       Methotrexate = 1305058,
                       Morphine = 1110410,
                       Sotalol = 1370109)

df <- dplyr::tbl(con, dbplyr::in_schema("results_schema", 
										"ud_egfr_with_select_drug_exposure")) %>%
  as_tibble() %>%
  filter(drug_concept_id %in% drugs_of_interest, 
         quantity > 0,
         between(year(drug_exposure_date), 2006, 2016)) %>%
  group_by(person_id, drug_exposure_date, drug_name) %>%
  mutate(daily_dose = sum(quantity),
         egfr = mean(egfr)) %>% # if >1 eGFR measurements that day
  distinct(person_id, drug_exposure_date, drug_name, .keep_all = TRUE) %>%
  group_by(drug_name) %>%
  filter(quantity > 0,
         daily_dose <= quantile(daily_dose, 0.75) * 1.5) %>%
  mutate(daily_dose_jittered = jitter(daily_dose, amount = diff(range(daily_dose)) / 100),
         daily_dose_jittered = ifelse(daily_dose_jittered < 0, 0, daily_dose_jittered),
         egfr_jittered = jitter(egfr, amount = 0.5)) %>%
  ungroup() %>%
  mutate(age_group = cut_width(age_at_exposure, width = 10, center = 5),
         year_of_exposure = as.numeric(substring(drug_exposure_date, 1, 4)),
         n_drugs = cut(n_drugs, 
                       breaks = c(0, 1, 5, 8, 12, 99), 
                       labels = c("1", "2-5", "6-8", "9-12", "13+"), 
                       ordered_result = TRUE))

# A stat convex hull for patient-specific areas
# -- see: https://ggplot2.tidyverse.org/articles/extending-ggplot2.html#creating-a-new-stat
StatChull <- ggproto("StatChull", Stat,
  compute_group = function(data, scales) {
    data[chull(data$x, data$y), , drop = FALSE]
  },
  required_aes = c("x", "y")
)
stat_chull <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  layer(
    stat = StatChull, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

# Define UI for application that draws scatter plots of daily doses and eGFRs
ui <- fluidPage(
  # Application title
  # titlePanel("Drug exposure vs. kidney function"),
  br(),
  sidebarLayout(
    sidebarPanel(
    	tags$head(
      	tags$style(type= "text/css", 
      			   "#inline_label label { display: table-cell; table-layout: fixed; width: 50%; white-space: no-wrap; text-align: left; font-weight: normal; vertical-align: middle; } 
									 #inline_label .form-group { display: table-row; }
									 #normal_label label { font-weight: normal; }")),
    	
    	h4("Data"),
    	selectInput("drug", "Select drug(s)", choices = as.list(drugs_of_interest), 
    				selected = 1110410, multiple = TRUE),
    	sliderInput("egfr_bound", "eGFR span", min = 0, max = 90, value = c(0, 60), step = 1, post = " ml/min"),
    	sliderInput("age_bound", "Age span", min = 0, max = 100, value = c(41, 100), step = 1, post = " years"),
    	# hr(),
    	fluidRow(
    		column(6, 
    			   h4("Distribution"),
    			   checkboxInput("show_median", "Medians", value = FALSE),
    			   checkboxInput("show_iqr", "Inter-percentile range", value = TRUE),
    			   tags$div(id = "normal_label",
    			   		 sliderInput("iqr_range", NULL, min = 2.5, max = 97.2, 
    			   		 			value = c(25, 75), step = 2.5, post = "%", ticks = FALSE)),
    			   tags$div(id = "inline_label",
    			   		 numericInput("egfr_bin_width", "Bin width", value = 10, min = 2.5, step = 2.5))),
    		column(6, 
    			   h4("Loess curve"),
    			   checkboxInput("show_smooth", "Show curve", value = TRUE),
    			   checkboxInput("show_smooth_se", "Confidence interval", value = FALSE),
    			   tags$div(id = "normal_label",
    			   		 sliderInput("smooth_se_level", NULL, min = 5, max = 99, 
    			   		 			value = 80, step = 5, post = "%", ticks = FALSE)),
    			   tags$div(id = "inline_label",
    			   		 numericInput("smooth_span", "Loess span", value = 0.75, min = 0.1, step = 0.05)))
    	),
    	# hr(),
    	br(),
    	fluidRow(
    		column(6, 
    			   h4("Patient convex hulls"),
    			   tags$div(id = "inline_label",
    			   		 numericInput("n_obs_per_patient", "Obs. / patient ≥", value = 30, min = 3)),
    			   tags$div(id = "normal_label",
    			   		 sliderInput("chull_alpha", "Opacity", min = 0, max = 100, 
    			   		 			value = 0, step = 5, post = "%", ticks = FALSE))),
    		column(6,
    			   h4("Plot appearance"),
    			   checkboxInput("show_obs", "Show data points", value = TRUE),
    			   checkboxInput("log_dose", "Daily dose on log-scale", value = FALSE),
    			   tags$div(id = "inline_label",
    			   		 numericInput("ncol", "Columns", value = 1, min = 1)))
    	),
    	# hr(),
    	h4("Tabulate by"),
    	fluidRow(
    		column(6, 
    			   selectInput("wrap_row", NULL, selected = "gender",
    			   			choices = list("None" = "none",
    			   						   "No. of drugs" = "n_drugs",
    			   						   "Patient sex" = "gender",
    			   						   "Age group" = "age_group"))),
    		column(6, 
    			   selectInput("wrap_col", NULL,
    			   			choices = list("None" = "none",
    			   						   "No. of drugs" = "n_drugs",
    			   						   "Patient sex" = "gender",
    			   						   "Age group" = "age_group")))
    	)
    ),
    mainPanel(plotOutput("distPlot", height = "800px", width = "100%"))
  )
)

# Define server logic required to draw plots
server <- function(input, output) {
	
	# Might be worthwhile to do plot caching
	
	# SECTION: CREATE REACTIVE EXPRESSIONS FOR COMPUTATIONAL PERFORMANCE
	ldf <- reactive({
		filter(df,
			   between(egfr, input$egfr_bound[1], input$egfr_bound[2]),
			   between(age_at_exposure, input$age_bound[1], input$age_bound[2]),
			   drug_concept_id %in% input$drug)
	})
	
	annotate_df <- reactive({
		df <- group_by(ldf(), drug_name) %>%
			mutate(label_y = max(daily_dose), 
				   q_offset = diff(range(daily_dose)) / 100)
		
		if (input$wrap_col != "none") 
			df <- group_by(df, !!sym(input$wrap_col), add = TRUE)
		if (input$wrap_row != "none") 
			df <- group_by(df, !!sym(input$wrap_row), add = TRUE)
		
		df
	})
	
	chull_df <- reactive({
		df <- group_by(ldf(), drug_name) %>%
			mutate(label_y = max(daily_dose), 
				   q_offset = diff(range(daily_dose)) / 100)
		
		if (input$wrap_col != "none") 
			df <- group_by(df, !!sym(input$wrap_col), add = TRUE)
		if (input$wrap_row != "none") 
			df <- group_by(df, !!sym(input$wrap_row), add = TRUE)
		
		group_by(df, person_id, add = TRUE) %>%
			add_tally() %>%
			filter(n >= input$n_obs_per_patient)
	})
	
	chull_label_df <- reactive({
		grouping_vars <- group_vars(chull_df())[-length(group_vars(chull_df()))]
		group_by(chull_df(), .dots = grouping_vars) %>%
			summarise(n_obs = n(),
					  n_patient = length(unique(person_id)),
					  label_y = label_y[1]) %>%
			mutate(label = sprintf("N = %s (%s chulls)", n_obs, n_patient))
	})
	
	label_df <- reactive({
		summarise(annotate_df(),
				  n_patient = length(unique(person_id)),
				  n_obs = n(),
				  label_y = label_y[1]) %>%
			mutate(label = sprintf("N = %s (%s pt's)", n_obs, n_patient))
	})
	
	ribbon_df <- reactive({
		mutate(annotate_df(), 
			   egfr_group = cut_width(egfr, width = input$egfr_bin_width,
			   					   center = input$egfr_bin_width / 2)) %>%
			group_by(egfr_group, add = TRUE) %>%
			summarise(lo = quantile(daily_dose, input$iqr_range[1] / 100) - q_offset[1], # add margin to cover points
					  median = quantile(daily_dose, 0.5),
					  hi = quantile(daily_dose, input$iqr_range[2] / 100) + q_offset[1]) %>%
			mutate(x = str_replace_all(substring(egfr_group, 2), c("\\)" = "", "\\]" = "")),
				   xmin = as.numeric(str_split_fixed(x, ",", 2)[, 1]),
				   xmax = as.numeric(str_split_fixed(x, ",", 2)[, 2])) %>%
			filter(!is.na(egfr_group)) %>%
			gather(x_pos, x, xmin, xmax) %>%
			arrange(egfr_group, x)
	})
	
	output$distPlot <- renderPlot({
		p <- ggplot(ldf(), aes(x = egfr_jittered, y = daily_dose_jittered)) +
			theme_minimal() +
			theme(text = element_text(size = 12), strip.text = element_text(size = 13),
				  plot.title = element_text(size = 14)) +
			labs(x = "Estimated glomerular filtration rate (eGFR)", y = "Estimated daily dose")
		
		if (input$chull_alpha > 0)
			p <- p + geom_label(aes(x = input$egfr_bound[1], y = label_y, label = label), chull_label_df(),
								hjust = "left", vjust = "top", colour = grey(0.6), size = 12 / ggplot2::.pt)
		else
			p <- p + geom_label(aes(x = input$egfr_bound[1], y = label_y, label = label), label_df(),
								hjust = "left", vjust = "top", colour = grey(0.6), size = 12 / ggplot2::.pt)
		
		if (isTRUE(input$show_obs)) {
			if (input$chull_alpha == 0) # if chull shown, show non-coloured data points
				p <- p + geom_point(aes(colour = drug_name), size = 0.5, alpha = 0.3, show.legend = FALSE) +
					scale_color_brewer(palette = "Set1")
			else
				p <- p + geom_point(size = 0.5, alpha = 0.3)
		}
		
		# How to wrap/grid
		if (input$wrap_col == "none" & input$wrap_row != "none")
			p <- p + facet_grid(drug_name ~ eval(as.name(input$wrap_row)), scales = "free_y")
		else if (input$wrap_col != "none" & input$wrap_row == "none")
			p <- p + facet_grid(drug_name ~ eval(as.name(input$wrap_col)), scales = "free_y")
		else if (input$wrap_col != "none" & input$wrap_row != "none" & length(input$drug) == 1)
			p <- p + facet_grid(eval(as.name(input$wrap_row)) ~ eval(as.name(input$wrap_col)),
								scales = "free_y") +
				labs(title = names(drugs_of_interest[drugs_of_interest == input$drug]))
		else
			p <- p + facet_wrap(~ drug_name, scales = "free_y", ncol = input$ncol)
		
		# Log-scale y axis? If not, keep (0, 0) in plot
		if (isTRUE(input$log_dose))
			p <- p + scale_y_log10()
		else
			p <- p + ylim(0, NA)
		
		# Show IQRs or medians?
		if (isTRUE(input$show_iqr))
			p <- p +
				geom_ribbon(aes(x = x, y = NULL, ymin = lo, ymax = hi, group = egfr_group), ribbon_df(),
							fill = "black", alpha = 0.1)
		if (isTRUE(input$show_median))
			p <- p + geom_line(aes(x = x, y = median, group = egfr_group), ribbon_df(), colour = "black", size = 0.5)
		
		# Show smoother?
		if (isTRUE(input$show_smooth))
			p <- p + stat_smooth(aes(x = egfr, y = daily_dose), method = "loess", se = input$show_smooth_se,
								 fill = "black", alpha = 0.1, span = input$smooth_span, size = 0.5, 
								 colour = "black", level = input$smooth_se_level / 100)
		
		# Show patients' convex hulls?
		if (input$chull_alpha > 0) 
			p <- p + stat_chull(aes(fill = as.factor(person_id)), chull_df(), 
								alpha = input$chull_alpha / 100, show.legend = FALSE) +
				geom_point(aes(colour = as.factor(person_id)), chull_df(), size = 0.5, alpha = 0.3, 
						   show.legend = FALSE) 
		
		# RETURN THE FINAL PLOT OBJECT
		# ggplotly(p) # for later deveopment into even more interactive plots
		p
	})
}

# Run the application 
shinyApp(ui = ui, server = server)

