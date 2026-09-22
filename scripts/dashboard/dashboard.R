# Ancient environmental genomics teaching dashboard
# Run in RStudio with: shiny::runApp("dashboard.R")

required_packages <- c("shiny", "bslib", "data.table", "ggplot2", "plotly", "DT", "scales", "rlang", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
  stop(
    "Install the required packages first: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(data.table)
  library(ggplot2)
  library(plotly)
  library(DT)
})

`%||%` <- function(x, fallback) {
  if (is.null(x) || length(x) == 0L || (is.character(x) && !nzchar(x[[1L]]))) fallback else x
}

app_file <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) NA_character_)
app_dir <- if (is.na(app_file)) normalizePath(getwd(), mustWork = TRUE) else dirname(app_file)
data_file <- file.path(app_dir, "merged_unicorn_metadmg_all_samples.tsv")

if (!file.exists(data_file)) {
  stop(
    "Could not find merged_unicorn_metadmg_all_samples.tsv next to dashboard.R.\n",
    "Put both files in the same directory and try again.",
    call. = FALSE
  )
}

course_data <- data.table::fread(
  data_file,
  sep = "\t",
  na.strings = c("", "NA", "NaN"),
  check.names = FALSE,
  showProgress = FALSE
)

# The course metadata sometimes arrives as character text. Keep this explicit:
course_data[, depth_cm := suppressWarnings(as.numeric(depth_cm))]
course_data[, latitude := suppressWarnings(as.numeric(depth_cm))]
course_data[, longitude := suppressWarnings(as.numeric(depth_cm))]
course_data[, years_bp := suppressWarnings(as.numeric(depth_cm))]

lineage <- as.character(course_data$taxa_path)
course_data[, kingdom := data.table::fcase(
  grepl("Bacteria", lineage, fixed = TRUE), "Bacteria",
  grepl("Archaea", lineage, fixed = TRUE), "Archaea",
  grepl("Fungi", lineage, fixed = TRUE), "Fungi",
  grepl("Viridiplantae", lineage, fixed = TRUE) |
    grepl("Streptophyta", lineage, fixed = TRUE), "Plant",
  grepl("Metazoa", lineage, fixed = TRUE), "Animal",
  default = "Other"
)]
course_data[, kingdom := factor(
  kingdom,
  levels = c("Bacteria", "Archaea", "Fungi", "Plant", "Animal", "Other")
)]
course_data[, .row_key := .I]

damage_columns <- grep("^(fw|bw)(K|N|f|dx|dxConf)[0-9]+$", names(course_data), value = TRUE)
numeric_columns <- names(course_data)[vapply(course_data, is.numeric, logical(1))]
numeric_columns <- setdiff(numeric_columns, c(".row_key", damage_columns))
filter_columns <- setdiff(names(course_data), c(".row_key", damage_columns, "taxa_path", "lca"))
categorical_columns <- filter_columns[!vapply(course_data[, ..filter_columns], is.numeric, logical(1))]

pretty_label <- function(x) {
  labels <- c(
    mean_alnani = "Mean ANI (%)",
    breath_cov = "Breadth of coverage",
    breath_ratio = "Breadth ratio",
    A = "Damage amplitude (A)",
    Zfit_new = "Damage significance (Zfit)",
    depth_cm = "Depth (cm)",
    sample_id = "Sample ID",
    name_unicorn = "Unicorn taxon name",
    rank_unicorn = "Unicorn rank"
  )
  if (x %in% names(labels)) labels[[x]] else tools::toTitleCase(gsub("_", " ", x, fixed = TRUE))
}

named_choices <- function(columns) {
  stats::setNames(columns, columns)
}

format_value <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(is.numeric(x), format(signif(x, 6), trim = TRUE, scientific = FALSE), as.character(x))
  )
}

log10_axis_labels <- scales::trans_format("log10", scales::math_format(10^.x))

compare_numeric <- function(x, operator, value) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(value)) return(rep(TRUE, length(x)))
  answer <- switch(
    operator,
    ">" = x > value,
    ">=" = x >= value,
    "==" = x == value,
    "<=" = x <= value,
    "<" = x < value,
    rep(TRUE, length(x))
  )
  answer[is.na(answer)] <- FALSE
  answer
}

build_damage_profile <- function(record) {
  make_end <- function(prefix, direction) {
    positions <- 0:14
    data.table(
      direction = direction,
      position_0based = positions,
      f = vapply(positions, function(i) suppressWarnings(as.numeric(record[[paste0(prefix, "f", i)]][[1L]])), numeric(1)),
      dx = vapply(positions, function(i) suppressWarnings(as.numeric(record[[paste0(prefix, "dx", i)]][[1L]])), numeric(1)),
      dx_conf = vapply(positions, function(i) suppressWarnings(as.numeric(record[[paste0(prefix, "dxConf", i)]][[1L]])), numeric(1))
    )
  }
  profile <- rbindlist(list(make_end("fw", "forward"), make_end("bw", "backward")))
  profile[is.finite(f) | is.finite(dx)]
}

make_damage_plot <- function(profile, y_max = 1) {
  if (!nrow(profile)) return(NULL)
  y_max <- suppressWarnings(as.numeric(y_max))
  if (!is.finite(y_max) || !(y_max %in% c(1, 0.8, 0.6, 0.4))) y_max <- 1
  profile[, `:=`(
    position = position_0based + 1L,
    lower = pmax(0, dx - dx_conf),
    upper = dx + dx_conf
  )]

  make_end <- function(data, backward = FALSE) {
    plot <- ggplot(data, aes(x = position)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.5) +
      geom_line(aes(y = dx), linetype = "dashed", alpha = 0.5) +
      geom_line(aes(y = f)) +
      geom_point(aes(y = f), shape = 21, fill = "white") +
      coord_cartesian(ylim = c(0, y_max)) +
      labs(x = "Position", y = if (backward) "G to A Frequency" else "C to T Frequency") +
      theme_classic()
    if (backward) {
      plot + scale_x_reverse(breaks = 1:15) + scale_y_continuous(position = "right")
    } else {
      plot + scale_x_continuous(breaks = 1:15)
    }
  }

  forward <- profile[direction == "forward"]
  backward <- profile[direction == "backward"]
  if (!nrow(forward)) return(ggplotly(make_end(backward, TRUE)))
  if (!nrow(backward)) return(ggplotly(make_end(forward, FALSE)))
  subplot(
    ggplotly(make_end(forward, FALSE)),
    ggplotly(make_end(backward, TRUE)),
    nrows = 1,
    shareY = TRUE,
    titleX = TRUE,
    margin = 0.06
  ) |>
    layout(showlegend = FALSE, margin = list(t = 10, b = 65, l = 65, r = 65)) |>
    config(displaylogo = FALSE, responsive = TRUE)
}

panel_handles <- function() {
  tagList(
    tags$span(class = "resize-handle resize-nw"),
    tags$span(class = "resize-handle resize-ne"),
    tags$span(class = "resize-handle resize-sw"),
    tags$span(class = "resize-handle resize-se")
  )
}

numeric_operator_choices <- c(
  "greater than" = ">", "greater than or equal to" = ">=", "equal to" = "==",
  "less than or equal to" = "<=", "less than" = "<"
)

filter_rule_ui <- function(id) {
  ns <- NS(id)
  div(
    class = "filter-sentence",
    selectInput(ns("action"), NULL, c("Keep" = "keep", "Remove" = "remove"), width = "95px"),
    tags$span("rows with"),
    selectInput(
      ns("variable"), NULL,
      choices = c("Choose variable" = "", named_choices(filter_columns)),
      selected = "", width = "240px"
    ),
    uiOutput(ns("rule_controls"), class = "filter-rule-controls"),
    actionButton(ns("remove"), NULL, icon = icon("xmark"), class = "btn-sm btn-outline-secondary remove-filter")
  )
}

filter_rule_server <- function(id, remove_callback) {
  moduleServer(id, function(input, output, session) {
    output$rule_controls <- renderUI({
      variable <- input$variable %||% ""
      if (!nzchar(variable) || !variable %in% names(course_data)) {
        return(tags$span(class = "filter-prompt", "choose a variable to complete this filter"))
      }
      if (variable %in% numeric_columns) {
        tagList(
          selectInput(session$ns("operator"), NULL, numeric_operator_choices, selected = ">=", width = "220px"),
          numericInput(session$ns("value"), NULL, value = NA, width = "130px")
        )
      } else {
        values <- sort(unique(as.character(course_data[[variable]])))
        values <- values[!is.na(values) & nzchar(values)]
        tagList(
          selectInput(
            session$ns("operator"), NULL,
            c("in" = "in", "not in" = "not_in", "equal to" = "equal"),
            selected = "in", width = "115px"
          ),
          selectizeInput(
            session$ns("value"), NULL, choices = values, selected = NULL,
            multiple = TRUE, width = "100%",
            options = list(placeholder = "choose one or more values", maxOptions = 1000, plugins = list("remove_button"))
          )
        )
      }
    })

    observeEvent(input$remove, remove_callback(), ignoreInit = TRUE)

    rule <- reactive({
      variable <- input$variable %||% ""
      list(
        variable = variable,
        action = input$action %||% "keep",
        operator = input$operator %||% if (variable %in% numeric_columns) ">=" else "in",
        value = input$value
      )
    })

    reset <- function() {
      updateSelectInput(session, "action", selected = "keep")
      updateSelectInput(session, "variable", selected = "")
    }

    list(rule = rule, reset = reset)
  })
}

plot_panel_ui <- function(id, number) {
  ns <- NS(id)
  numeric_axis_choices <- c("Choose a numerical column" = "", named_choices(numeric_columns))
  aesthetic_choices <- c("Choose a variable" = "", named_choices(filter_columns))

  card(
    class = "builder-card",
    card_header(
      div(class = "plot-card-title", paste("Plot", number)),
      actionButton(ns("remove_plot"), "Remove plot", class = "btn-sm btn-outline-danger")
    ),
    div(
      class = "plot-panel-body",
      div(
        class = "plot-options-grid",
        div(
          class = "option-section main-options",
          h5("Main plot"),
          selectInput(
            ns("plot_type"), "Plot type",
            c("Choose a plot type" = "", "Scatter plot" = "scatter", "Histogram" = "histogram"),
            selected = ""
          ),
          conditionalPanel(
            condition = "input.plot_type != ''", ns = ns,
            selectInput(ns("x_var"), "X axis", numeric_axis_choices, selected = "")
          ),
          conditionalPanel(
            condition = "input.plot_type == 'scatter'", ns = ns,
            selectInput(ns("y_var"), "Y axis", numeric_axis_choices, selected = "")
          ),
          conditionalPanel(
            condition = "input.plot_type == 'histogram'", ns = ns,
            numericInput(ns("bins"), "Number of bins", value = 50, min = 1, max = 500, step = 1)
          ),
          conditionalPanel(
            condition = "input.plot_type != ''", ns = ns,
            div(
              class = "inline-checks",
              checkboxInput(ns("log_x"), "Log10 X axis", FALSE),
              checkboxInput(ns("log_y"), "Log10 Y axis", FALSE)
            )
          )
        ),
        div(
          class = "option-section aesthetics-options",
          h5("Aesthetics"),
          conditionalPanel(
            condition = "input.plot_type == ''", ns = ns,
            p(class = "section-prompt", "Choose a plot type to see its aesthetics.")
          ),
          conditionalPanel(
            condition = "input.plot_type != ''", ns = ns,
            selectInput(
              ns("colour_mode"), "Outline colour",
              c("Constant" = "constant", "Pass / fail" = "status", "A variable" = "variable"),
              selected = "constant"
            ),
            conditionalPanel(
              condition = "input.colour_mode == 'variable'", ns = ns,
              selectInput(ns("colour_var"), "Colour variable", aesthetic_choices, selected = "")
            ),
            selectInput(
              ns("fill_mode"), "Fill",
              c("Constant" = "constant", "Pass / fail" = "status", "A variable" = "variable"),
              selected = "constant"
            ),
            conditionalPanel(
              condition = "input.fill_mode == 'variable'", ns = ns,
              selectInput(ns("fill_var"), "Fill variable", aesthetic_choices, selected = "")
            ),
            conditionalPanel(
              condition = "input.plot_type == 'scatter'", ns = ns,
              selectInput(
                ns("size_mode"), "Point size",
                c("Constant" = "constant", "A numeric variable" = "variable"),
                selected = "constant"
              ),
              conditionalPanel(
                condition = "input.size_mode == 'variable'", ns = ns,
                selectInput(
                  ns("size_var"), "Size variable",
                  c("Choose a numerical column" = "", named_choices(numeric_columns)), selected = ""
                )
              ),
              conditionalPanel(
                condition = "input.size_mode == 'constant'", ns = ns,
                numericInput(ns("point_size"), "Constant point size", 4, min = 0.5, max = 20, step = 0.5)
              )
            )
          )
        ),
        div(
          class = "option-section filter-options",
          h5("Filtering"),
          conditionalPanel(
            condition = "input.plot_type == ''", ns = ns,
            p(class = "section-prompt", "Choose a plot type before adding filters.")
          ),
          conditionalPanel(
            condition = "input.plot_type != ''", ns = ns,
            radioButtons(
              ns("filter_mode"), "Axis thresholds",
              c("Visualise pass / fail" = "visualise", "Apply and remove failed rows" = "apply"),
              selected = "visualise", inline = TRUE
            ),
            div(
              class = "axis-filter-grid",
              div(
                class = "axis-threshold-control",
                checkboxInput(ns("x_threshold_on"), "Add X threshold", FALSE),
                conditionalPanel(
                  condition = "input.x_threshold_on", ns = ns,
                  div(
                    class = "threshold-row",
                    selectInput(ns("x_threshold_op"), NULL, numeric_operator_choices, selected = ">=", width = "68%"),
                    numericInput(ns("x_threshold"), NULL, value = NA, width = "30%")
                  )
                )
              ),
              conditionalPanel(
                condition = "input.plot_type == 'scatter'", ns = ns,
                div(
                  class = "axis-threshold-control",
                  checkboxInput(ns("y_threshold_on"), "Add Y threshold", FALSE),
                  conditionalPanel(
                    condition = "input.y_threshold_on", ns = ns,
                    div(
                      class = "threshold-row",
                      selectInput(ns("y_threshold_op"), NULL, numeric_operator_choices, selected = ">=", width = "68%"),
                      numericInput(ns("y_threshold"), NULL, value = NA, width = "30%")
                    )
                  )
                )
              )
            ),
            hr(),
            div(
              class = "filter-toolbar",
              div(
                tags$strong("Additional filters"),
                p(class = "small text-muted", "These always remove matching or non-matching rows before plotting.")
              ),
              actionButton(ns("add_filter"), "Add another filter", icon = icon("plus"), class = "btn-sm btn-outline-primary")
            ),
            div(
              class = "filter-output",
              tagList(lapply(seq_len(8L), function(filter_id) {
                div(
                  id = ns(paste0("filter_slot_", filter_id)),
                  style = "display:none;",
                  filter_rule_ui(ns(paste0("filter_", filter_id)))
                )
              }))
            )
          )
        )
      ),
      div(
        class = "plot-results",
        card(
          class = "resizable-panel main-plot-panel",
          card_header(
            tags$strong(textOutput(ns("plot_heading"), inline = TRUE)),
            downloadButton(ns("download_plot"), "PNG", class = "btn-sm btn-outline-secondary")
          ),
          plotlyOutput(ns("main_plot"), height = "540px"),
          panel_handles()
        ),
        conditionalPanel(
          condition = "input.plot_type == 'scatter'", ns = ns,
          layout_columns(
            col_widths = c(8, 4),
            card(
              class = "resizable-panel damage-panel",
              card_header(tags$strong("Damage profile (scatter selection)")),
              tags$p(class = "panel-description", "Observed and fitted terminal substitutions for the clicked point."),
              plotlyOutput(ns("damage_plot"), height = "360px"),
              panel_handles()
            ),
            card(
              class = "selected-panel",
              card_header(tags$strong("Selected point")),
              DTOutput(ns("selected_table"))
            )
          )
        )
      )
    )
  )
}

plot_panel_server <- function(id, number, remove_callback) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    filter_ids <- reactiveVal(integer())
    selected_key <- reactiveVal(NULL)

    observeEvent(input$remove_plot, remove_callback(), ignoreInit = TRUE)

    remove_filter <- function(filter_id) {
      force(filter_id)
      function() {
        filter_ids(setdiff(filter_ids(), filter_id))
        session$sendCustomMessage(
          "toggle-plot-panel",
          list(id = ns(paste0("filter_slot_", filter_id)), show = FALSE)
        )
      }
    }

    filter_modules <- lapply(seq_len(8L), function(filter_id) {
      filter_rule_server(paste0("filter_", filter_id), remove_filter(filter_id))
    })

    observeEvent(input$add_filter, {
      available <- setdiff(seq_len(8L), filter_ids())
      if (!length(available)) {
        showNotification("Each plot can have up to eight additional filters.", type = "warning")
        return()
      }
      filter_id <- available[[1L]]
      filter_modules[[filter_id]]$reset()
      filter_ids(c(filter_ids(), filter_id))
      session$sendCustomMessage(
        "toggle-plot-panel",
        list(id = ns(paste0("filter_slot_", filter_id)), show = TRUE)
      )
    })

    read_additional_filters <- reactive({
      lapply(filter_ids(), function(filter_id) filter_modules[[filter_id]]$rule())
    })

    base_filtered_data <- reactive({
      data <- copy(course_data)
      for (rule in read_additional_filters()) {
        if (!nzchar(rule$variable %||% "") || !rule$variable %in% names(data)) next
        if (is.null(rule$value) || !length(rule$value) || all(is.na(rule$value))) next
        x <- data[[rule$variable]]
        if (rule$variable %in% numeric_columns) {
          match <- compare_numeric(x, rule$operator, suppressWarnings(as.numeric(rule$value[[1L]])))
        } else {
          match <- x %in% as.character(rule$value)
          if (identical(rule$operator, "not_in")) match <- !match
          match[is.na(match)] <- FALSE
        }
        keep <- if (identical(rule$action, "remove")) !match else match
        data <- data[keep]
      }
      data
    })

    thresholded_data <- reactive({
      req(input$plot_type %in% c("scatter", "histogram"))
      data <- base_filtered_data()
      x <- input$x_var %||% ""
      req(x %in% numeric_columns)
      pass <- rep(TRUE, nrow(data))
      if (isTRUE(input$x_threshold_on)) {
        pass <- pass & compare_numeric(data[[x]], input$x_threshold_op %||% ">=", suppressWarnings(as.numeric(input$x_threshold)))
      }
      if (identical(input$plot_type, "scatter")) {
        y <- input$y_var %||% ""
        req(y %in% numeric_columns)
        if (isTRUE(input$y_threshold_on)) {
          pass <- pass & compare_numeric(data[[y]], input$y_threshold_op %||% ">=", suppressWarnings(as.numeric(input$y_threshold)))
        }
      }
      data[, .filter_status := factor(ifelse(pass, "Pass", "Fail"), levels = c("Pass", "Fail"))]
      if (identical(input$filter_mode, "apply")) data <- data[.filter_status == "Pass"]
      data
    })

    hover_columns <- reactive({
      filter_variables <- vapply(read_additional_filters(), function(rule) rule$variable %||% "", character(1))
      unique(c(
        "sample_id", "depth_cm", "name_unicorn", "rank_unicorn", "kingdom", "Zfit_new", "A",
        input$x_var,
        if (identical(input$plot_type, "scatter")) input$y_var else NULL,
        filter_variables[nzchar(filter_variables)]
      ))
    })

    make_hover <- function(data) {
      columns <- intersect(hover_columns(), names(data))
      labels <- vapply(columns, pretty_label, character(1))
      vapply(seq_len(nrow(data)), function(row) {
        pieces <- vapply(seq_along(columns), function(index) {
          value <- htmltools::htmlEscape(format_value(data[[columns[[index]]]][[row]]))
          paste0("<b>", labels[[index]], ":</b> ", value)
        }, character(1))
        paste(pieces, collapse = "<br>")
      }, character(1))
    }

    main_ggplot <- reactive({
      data <- thresholded_data()
      x <- input$x_var
      req(x %in% numeric_columns)
      if (isTRUE(input$log_x)) data <- data[is.finite(get(x)) & get(x) > 0]
      if (identical(input$plot_type, "scatter") && isTRUE(input$log_y)) {
        y_for_log <- input$y_var
        req(y_for_log %in% numeric_columns)
        data <- data[is.finite(get(y_for_log)) & get(y_for_log) > 0]
      }
      validate(need(nrow(data) > 0L, "No positive finite rows remain for the selected log axes."))
      colour_mode <- input$colour_mode %||% "constant"
      fill_mode <- input$fill_mode %||% "constant"
      if (colour_mode == "variable") {
        validate(need((input$colour_var %||% "") %in% filter_columns, "Choose an outline-colour variable."))
      }
      if (fill_mode == "variable") {
        validate(need((input$fill_var %||% "") %in% filter_columns, "Choose a fill variable."))
      }

      if (identical(input$plot_type, "histogram")) {
        mapping <- aes(x = .data[[x]], y = after_stat(count))
        if (colour_mode == "status") mapping$colour <- rlang::quo(.data$.filter_status)
        if (colour_mode == "variable" && input$colour_var %in% names(data)) mapping$colour <- rlang::quo(.data[[!!input$colour_var]])
        if (fill_mode == "status") mapping$fill <- rlang::quo(.data$.filter_status)
        if (fill_mode == "variable" && input$fill_var %in% names(data)) mapping$fill <- rlang::quo(.data[[!!input$fill_var]])

        histogram_args <- list(bins = max(1L, as.integer(input$bins %||% 50L)), alpha = 0.88, na.rm = TRUE)
        if (colour_mode == "constant") histogram_args$colour <- "white"
        if (fill_mode == "constant") histogram_args$fill <- "#3f8f8a"

        plot <- ggplot(data, mapping) +
          do.call(geom_histogram, histogram_args) +
          labs(x = pretty_label(x), y = "Count", colour = NULL, fill = NULL) +
          theme_bw() +
          theme(panel.grid.minor = element_blank(), legend.position = "bottom")
        if (identical(input$filter_mode, "visualise") && isTRUE(input$x_threshold_on) && is.finite(as.numeric(input$x_threshold))) {
          plot <- plot + geom_vline(xintercept = as.numeric(input$x_threshold), linetype = "dotted", linewidth = 0.8)
        }
      } else {
        y <- input$y_var
        req(y %in% numeric_columns)
        data[, .hover := make_hover(data)]
        mapping <- aes(x = .data[[x]], y = .data[[y]], key = .row_key, text = .hover)
        size_mode <- input$size_mode %||% "constant"
        if (size_mode == "variable") {
          validate(need((input$size_var %||% "") %in% numeric_columns, "Choose a point-size variable."))
        }
        if (colour_mode == "status") mapping$colour <- rlang::quo(.data$.filter_status)
        if (colour_mode == "variable" && input$colour_var %in% names(data)) mapping$colour <- rlang::quo(.data[[!!input$colour_var]])
        if (fill_mode == "status") mapping$fill <- rlang::quo(.data$.filter_status)
        if (fill_mode == "variable" && input$fill_var %in% names(data)) mapping$fill <- rlang::quo(.data[[!!input$fill_var]])
        if (size_mode == "variable" && input$size_var %in% names(data)) mapping$size <- rlang::quo(.data[[!!input$size_var]])

        point_args <- list(shape = 21, stroke = 0.1, alpha = 0.82, na.rm = TRUE)
        if (colour_mode == "constant") point_args$colour <- "black"
        if (fill_mode == "constant") point_args$fill <- "NA"
        if (size_mode == "constant") point_args$size <- input$point_size %||% 4

        plot <- ggplot(data, mapping) +
          do.call(geom_point, point_args) +
          labs(x = pretty_label(x), y = pretty_label(y), colour = NULL, fill = NULL, size = NULL) +
          theme_bw() +
          theme(panel.grid.minor = element_line(colour = "grey92"), legend.position = "bottom")
        if (identical(input$filter_mode, "visualise") && isTRUE(input$x_threshold_on) && is.finite(as.numeric(input$x_threshold))) {
          plot <- plot + geom_vline(xintercept = as.numeric(input$x_threshold), linetype = "dotted", linewidth = 0.8)
        }
        if (identical(input$filter_mode, "visualise") && isTRUE(input$y_threshold_on) && is.finite(as.numeric(input$y_threshold))) {
          plot <- plot + geom_hline(yintercept = as.numeric(input$y_threshold), linetype = "dotted", linewidth = 0.8)
        }
      }

      if (isTRUE(input$log_x)) plot <- plot + scale_x_log10(labels = log10_axis_labels)
      if (isTRUE(input$log_y)) plot <- plot + scale_y_log10(labels = log10_axis_labels)
      if (colour_mode == "status") {
        plot <- plot + scale_colour_manual(values = c("Pass" = "#2E8B57", "Fail" = "#C44A3A"), drop = FALSE)
      }
      if (fill_mode == "status") {
        plot <- plot + scale_fill_manual(values = c("Pass" = "#61B878", "Fail" = "#E36A5D"), drop = FALSE)
      } else if(fill_mode == "variable") {
        plot <- plot + scale_fill_viridis_c(option = "viridis")
      }
      plot
    })

    output$main_plot <- renderPlotly({
      validate(need(input$plot_type %in% c("scatter", "histogram"), "Choose a plot type."))
      validate(need((input$x_var %||% "") %in% numeric_columns, "Choose an X-axis variable."))
      if (identical(input$plot_type, "scatter")) {
        validate(need((input$y_var %||% "") %in% numeric_columns, "Choose a Y-axis variable."))
      }
      plot <- main_ggplot()
      if (identical(input$plot_type, "scatter")) {
        ggplotly(plot, tooltip = "text", source = paste0("teaching_scatter_", number)) |>
          config(displaylogo = FALSE, responsive = TRUE, modeBarButtonsToRemove = c("lasso2d", "select2d")) |>
          event_register("plotly_click")
      } else {
        ggplotly(plot, source = paste0("teaching_histogram_", number)) |>
          config(displaylogo = FALSE, responsive = TRUE)
      }
    })

    click_event <- reactive({
      event_id <- paste0("plotly_click-teaching_scatter_", number)
      raw_event <- session$rootScope()$input[[event_id]]
      if (is.null(raw_event)) return(NULL)
      jsonlite::parse_json(raw_event, simplifyVector = TRUE)
    })

    observeEvent(click_event(), {
      click <- click_event()
      if (!is.null(click$key)) selected_key(as.integer(click$key))
    }, ignoreInit = TRUE)

    observeEvent(list(input$plot_type, input$x_var, input$y_var), selected_key(NULL), ignoreInit = TRUE)

    selected_record <- reactive({
      req(identical(input$plot_type, "scatter"), !is.null(selected_key()))
      record <- course_data[.row_key == selected_key()]
      req(nrow(record) == 1L)
      record
    })

    output$damage_plot <- renderPlotly({
      record <- selected_record()
      profile <- build_damage_profile(record)
      plot <- make_damage_plot(profile, 1)
      validate(need(!is.null(plot), "No per-position damage profile is available for this point."))
      plot
    })

    output$selected_table <- renderDT({
      record <- selected_record()
      columns <- intersect(hover_columns(), names(record))
      details <- data.frame(
        Field = vapply(columns, pretty_label, character(1)),
        Value = vapply(columns, function(column) format_value(record[[column]][[1L]]), character(1)),
        check.names = FALSE
      )
      datatable(details, rownames = FALSE, options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollY = "340px"), escape = TRUE)
    }, server = FALSE)

    output$plot_heading <- renderText({
      if (!input$plot_type %in% c("scatter", "histogram")) return("Choose a plot type")
      if (!(input$x_var %||% "") %in% numeric_columns) return("Choose an X-axis variable")
      if (identical(input$plot_type, "scatter") && !(input$y_var %||% "") %in% numeric_columns) return("Choose a Y-axis variable")
      data <- thresholded_data()
      sprintf("%s — %s rows", if (identical(input$plot_type, "histogram")) "Histogram" else "Scatter plot", scales::comma(nrow(data)))
    })

    output$download_plot <- downloadHandler(
      filename = function() sprintf("aeGenomics_plot_%02d.png", number),
      content = function(file) ggsave(file, plot = main_ggplot(), width = 11, height = 7, dpi = 300, bg = "white")
    )
  })
}

app_css <- HTML("
:root { --ae-teal:#176b68; --ae-deep:#153c43; --ae-bg:#f4f7f7; --ae-muted:#65757a; }
body { background:var(--ae-bg); }
.navbar { box-shadow:0 2px 10px rgba(16,53,58,.14); }
.app-intro { padding:1rem 1.2rem .1rem; }
.app-intro h2 { color:var(--ae-deep); margin-bottom:.3rem; }
.toolbar-card { margin:0 1.2rem 1rem; min-height:64px; background:white; border:1px solid #d9e2e1; box-shadow:0 2px 8px rgba(26,65,68,.06); }
.toolbar-card .card-body { display:flex; align-items:center; padding:.7rem 1rem; }
.plots-stack { display:flex; flex-direction:column; gap:1rem; margin:0 1.2rem 1rem; min-width:0; }
.builder-card { border:1px solid #d9e2e1; box-shadow:0 3px 12px rgba(26,65,68,.08); }
.builder-card > .card-header { background:#eaf2f1; display:flex; align-items:center; justify-content:space-between; }
.plot-card-title { font-size:1.05rem; font-weight:700; color:var(--ae-deep); }
.plot-panel-body { display:flex; flex-direction:column; }
.plot-options-grid { display:grid; grid-template-columns:minmax(0,1fr) minmax(0,1fr); gap:1rem; padding:1rem; background:#f7faf9; border-bottom:1px solid #dce5e3; }
.option-section { background:white; border:1px solid #d9e4e2; border-radius:.55rem; padding:1rem 1rem .65rem; box-shadow:0 1px 4px rgba(26,65,68,.04); }
.option-section h5 { color:var(--ae-deep); font-weight:700; margin:0 0 .8rem; padding-bottom:.5rem; border-bottom:2px solid #dfecea; }
.filter-options { grid-column:1 / -1; }
.section-prompt, .filter-prompt { color:var(--ae-muted); font-size:.9rem; font-style:italic; }
.inline-checks { display:flex; gap:1.5rem; flex-wrap:wrap; }
.inline-checks .form-check { margin-bottom:.4rem; }
.axis-filter-grid { display:grid; grid-template-columns:repeat(2,minmax(240px,1fr)); gap:1rem; }
.filter-toolbar { display:flex; align-items:center; justify-content:space-between; gap:1rem; }
.filter-toolbar p { margin:0; }
.filter-output { margin-top:.5rem; }
.plot-results { min-width:0; display:flex; flex-direction:column; gap:1rem; }
.resizable-panel { position:relative; min-width:420px; min-height:340px; overflow:auto; }
.main-plot-panel { height:625px; }
.damage-panel { height:470px; }
.selected-panel { height:470px; overflow:hidden; }
.panel-description { color:var(--ae-muted); font-size:.88rem; margin:.55rem 1rem 0; }
.threshold-row { display:flex; gap:2%; align-items:flex-start; }
.axis-threshold-control { min-width:0; }
.axis-threshold-control .form-check { margin-bottom:.35rem; }
.threshold-row .form-group { margin-bottom:.4rem; }
.filter-sentence { background:white; border:1px solid #dce5e3; border-radius:.45rem; padding:.55rem; margin:.6rem 0; display:flex; flex-wrap:wrap; gap:.3rem; align-items:flex-start; }
.filter-sentence > span { padding-top:.42rem; white-space:nowrap; }
.filter-sentence .form-group { margin:0; }
.filter-sentence .selectize-control { margin-bottom:0; min-width:220px; }
.filter-rule-controls { display:flex; flex:1 1 360px; min-width:300px; gap:.35rem; align-items:flex-start; }
.filter-rule-controls > .form-group:last-child { flex:1 1 auto; }
.filter-rule-controls > .filter-prompt { padding-top:.42rem; }
.remove-filter { margin-top:.15rem; }
.data-card { margin:0 1.2rem 1.5rem; }
.data-card .card-body { overflow-x:auto; }
.resize-handle { position:absolute; z-index:1000; width:18px; height:18px; }
.resize-nw { top:0; left:0; cursor:nwse-resize; }
.resize-ne { top:0; right:0; cursor:nesw-resize; }
.resize-sw { bottom:0; left:0; cursor:nesw-resize; }
.resize-se { bottom:0; right:0; cursor:nwse-resize; }
.resize-handle::after { content:''; position:absolute; width:8px; height:8px; border-color:rgba(23,107,104,.6); border-style:solid; }
.resize-nw::after { top:4px; left:4px; border-width:2px 0 0 2px; }
.resize-ne::after { top:4px; right:4px; border-width:2px 2px 0 0; }
.resize-sw::after { bottom:4px; left:4px; border-width:0 0 2px 2px; }
.resize-se::after { bottom:4px; right:4px; border-width:0 2px 2px 0; }
@media (max-width:900px) {
  .plot-options-grid { grid-template-columns:1fr; }
  .filter-options { grid-column:auto; }
  .axis-filter-grid { grid-template-columns:1fr; }
  .filter-toolbar { align-items:flex-start; flex-direction:column; }
  .resizable-panel { min-width:320px; }
}
")

resize_js <- HTML("
(function(){
  'use strict';
  var drag=null;
  $(document).on('shiny:connected',function(){
    Shiny.addCustomMessageHandler('toggle-plot-panel',function(message){
      var slot=document.getElementById(message.id);
      if(!slot)return;
      slot.style.display=message.show?'block':'none';
      if(message.show){
        window.setTimeout(function(){
          var graphs=slot.querySelectorAll('.plotly.html-widget');
          graphs.forEach(function(graph){if(window.Plotly)Plotly.Plots.resize(graph);});
        },80);
      }
    });
  });
  document.addEventListener('mousedown',function(event){
    var handle=event.target.closest('.resize-handle');
    if(!handle)return;
    event.preventDefault();
    var panel=handle.closest('.resizable-panel'), rect=panel.getBoundingClientRect();
    drag={panel:panel,handle:handle,x:event.clientX,y:event.clientY,width:rect.width,height:rect.height,tx:Number(panel.dataset.tx||0),ty:Number(panel.dataset.ty||0)};
    document.body.style.userSelect='none';
  });
  document.addEventListener('mousemove',function(event){
    if(!drag)return;
    var dx=event.clientX-drag.x,dy=event.clientY-drag.y;
    var left=drag.handle.classList.contains('resize-nw')||drag.handle.classList.contains('resize-sw');
    var top=drag.handle.classList.contains('resize-nw')||drag.handle.classList.contains('resize-ne');
    var width=Math.max(620,drag.width+(left?-dx:dx));
    var height=Math.max(340,drag.height+(top?-dy:dy));
    var tx=drag.tx+(left?drag.width-width:0),ty=drag.ty+(top?drag.height-height:0);
    drag.panel.style.width=width+'px';drag.panel.style.height=height+'px';
    drag.panel.style.transform='translate('+tx+'px,'+ty+'px)';drag.panel.dataset.tx=tx;drag.panel.dataset.ty=ty;
    var graph=drag.panel.querySelector('.plotly.html-widget');if(graph&&window.Plotly)Plotly.Plots.resize(graph);
  });
  document.addEventListener('mouseup',function(){drag=null;document.body.style.userSelect='';});
})();
")

ui <- page_navbar(
  title = "aeGenomics2026",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#176b68"),
  header = tagList(tags$style(app_css), tags$script(resize_js)),
  nav_panel(
    "Explore",
    div(
      class = "app-intro",
      h2("Explore merged Unicorn and metaDMG results"),
      p(
        "Build up to ten plots. Time to explore the data!"
      )
    ),
    div(
      class = "plots-stack",
      tagList(lapply(seq_len(10L), function(number) {
        div(
          id = paste0("plot_slot_", number),
          class = "plot-slot",
          style = "display:none;",
          plot_panel_ui(paste0("plot_", number), number)
        )
      }))
    ),
    card(
      class = "toolbar-card",
      card_body(
        div(
          class = "d-flex align-items-center gap-3",
          actionButton("add_plot", "Add plot panel", icon = icon("plus"), class = "btn-primary"),
          textOutput("plot_count", inline = TRUE)
        )
      )
    ),
    card(
      class = "data-card",
      card_header(tags$strong("Merged data table")),
      card_body(
        DTOutput("data_table")
      )
    )
  )
)

server <- function(input, output, session) {
  plot_ids <- reactiveVal(integer())

  remove_plot <- function(number) {
    force(number)
    function() {
      plot_ids(setdiff(plot_ids(), number))
      session$sendCustomMessage(
        "toggle-plot-panel",
        list(id = paste0("plot_slot_", number), show = FALSE)
      )
    }
  }

  lapply(seq_len(10L), function(number) {
    plot_panel_server(paste0("plot_", number), number, remove_plot(number))
  })

  observeEvent(input$add_plot, {
    if (length(plot_ids()) >= 10L) {
      showNotification("The dashboard supports a maximum of ten plots.", type = "warning")
      return()
    }
    number <- setdiff(seq_len(10L), plot_ids())[[1L]]
    plot_ids(c(plot_ids(), number))
    session$sendCustomMessage(
      "toggle-plot-panel",
      list(id = paste0("plot_slot_", number), show = TRUE)
    )
  })

  output$plot_count <- renderText({
    sprintf("%d of 10 plot panels", length(plot_ids()))
  })

  output$data_table <- renderDT({
    shown_columns <- setdiff(names(course_data), c(".row_key", damage_columns, "lca", "taxa_path"))
    datatable(
      course_data[, ..shown_columns],
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 15, scrollX = TRUE, deferRender = TRUE, autoWidth = FALSE)
    )
  }, server = TRUE)
}

shinyApp(ui, server)
