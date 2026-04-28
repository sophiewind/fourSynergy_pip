#devtools::install_github('https://github.com/sophiewind/fourSynergy')
library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyBS)
library(shinyvalidate)
library(dplyr)
library(ggplot2)
library(reshape2)
library(DT)
library(stringr)
library(tidyr)
library(fourSynergy)
library(bslib)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(cowplot)
library(UpSetR)
library(gVenn)
library(shinycustomloader)
library(DESeq2)  #
library(bamsignals)  #
library(karyoploteR)  #
library(tibble)  #
library(plotly)  #
library(yaml)
library(clusterProfiler)
source('.///utils.R')
library(fresh)

# ui --------------------------------------------------------------------------
# Create the theme
mytheme <- create_theme(
    adminlte_color(
        light_blue = "#2C2B2B",  # top right
        olive = "#2C2B2B",
        green = "#2C2B2B",
        lime = "#E6E6E6"
    ),    adminlte_sidebar(
        width = "400px",
        dark_bg = "#2C2B2B",
        dark_hover_bg = "#1E1E1E",
        dark_color = "#FFF",
        dark_hover_color = "#A4FFB3",
    ),
    adminlte_global(
        content_bg = "#FFF",
        box_bg = "#E6E6E6",
        info_box_bg = "#2C2B2B"
    )
)

# Definition der UI
ui <- shinyUI(
    dashboardPage(
        title = "fourSynergy",
        dashboardHeader(
            title = tags$div(
                tags$img(src = "logo_wide_2.png",
                         height = "50px", style = "margin-right: 10px;")
            )
        ),

        dashboardSidebar(
            sidebarMenu(
                tags$style(HTML(".sidebar-menu li a { font-size: 22px; }")),
                menuItem("Upload", tabName = "upload", icon = icon("upload")),
                menuItem("QC", tabName = "qc", icon = icon("chart-simple")),
                menuItem("Base Tools", tabName = "base_tools", icon = icon("dice-one")),
                menuItem("Ensemble", tabName = "ensemble", icon = icon("dice-four")),
                menuItem("Differential Analysis", tabName = "diff_analysis",
                         icon = icon("dna")),
                menuItem("Pipeline and Support", tabName = "pipeline",
                         icon = icon("circle-info")))),

        dashboardBody(
            use_theme(mytheme),
            useShinyjs(),
            tags$style(HTML("
                 .value-box-title {
                    font-size: 16px!important;
                    line-height: 1.2!important;
                  }
                 .irs-bar {background: #2D973F!important;}
                 .irs-handle {background: #2D973F; border: 1px solid #2D973F}
                 .irs-bar-edge {background: #2D973F}
                 .irs-single {background: #2D973F}
                 .irs-min,.irs-max {background: #2D973F}
                 .irs-tooltip {background-color: #007bff; color: white}
                 .irs-grid-text { font-size: 10px!important; }
                  [class*='js-irs-'].irs-single {background: #2D973F}
                  .progress-bar{background-color:#2D973F;}
                ")),
            tags$style(HTML(".js-irs-0 .irs-single {background: #2D973F}")),
            tags$style(HTML(".js-irs-1 .irs-single {background: #2D973F}")),
            tags$style(HTML(".js-irs-2 .irs-single {background: #2D973F}")),

            tabItems(
                ## Upload ----------------------------------
                tabItem(tabName = "upload",
                        h1(icon("upload"), "Upload files"),
                        hr(),
                        fluidRow(
                            column(12,
                                   h2("Upload config",
                                      bsButton(inputId = "pa",
                                               label = "",
                                               icon = icon("info"),
                                               size = "small")),
                                   bsPopover(id = "pa",
                                             title = "Config file",
                                             content = paste0(
                                                 "You can find more information to the config file info.yaml here: https://github.com/sophiewind/fourSynergy_pip"
                                             ),
                                             placement = "right",
                                             trigger = "hover",
                                             options = list(container = "body")
                                   ),
                                   fileInput("config", "Upload config (Datasets/[projectname]/info.yaml):", accept = ".yaml", width = 800),
                                   textOutput("config_valid_val"),
                                   conditionalPanel(
                                       condition = "output.config_status == 1",
                                       h2("fourSynergy pipeline results upload",
                                          bsButton(inputId = "pip",
                                                   label = "",
                                                   icon = icon("info"),
                                                   size = "small")),
                                       bsPopover(id = "pip",
                                                 title = "Pipeline results",
                                                 content = paste0(
                                                     "snakemake generates the shiny_in directory as long as the `--void_shiny` flag is not active. If any files are missing, it will be reported here. While missing quality control files can be ignored, the following files are required: sia files, BAM files, and virtual fragment library-related files."
                                                 ),
                                                 placement = "right",
                                                 trigger = "hover",
                                                 options = list(container = "body")
                                       ),
                                       fileInput("shiny_in", "Upload all files from (results/[projectname]/shiny_in/):", multiple = TRUE, width = 800),
                                       uiOutput("shiny_in_info"),
                                       uiOutput("multiqc_info"),
                                       uiOutput("sia_info"),
                                       uiOutput("track_info"),
                                       uiOutput("basic_info"),
                                   ),
                                   conditionalPanel(
                                       condition = "output.config_status == -1",
                                       fluidRow(
                                           style = "text-align: center; margin-top: 50px;",
                                           p(
                                               style = "color: black;",
                                               "Config file not valid."
                                           ),
                                           uiOutput("error_message_config")
                                       )),
                                   conditionalPanel(
                                       condition = "output.config_status == 0",
                                       fluidRow(
                                           style = "text-align: center; margin-top: 50px;",
                                           p(
                                               style = "color: black;",
                                               "Please upload a config file."
                                           )
                                       ))),

                        )
                ),

                tabItem(tabName = "qc",
                        conditionalPanel(
                            condition = "output.file_ready",
                            h1(icon("chart-simple"), "Quality control"),
                            hr(),
                            fluidRow(
                                column(12,
                                       h2(icon("database"), "Metadata"),
                                       h3(textOutput('sample_name')),
                                       fluidRow(column(6,
                                                       value_box(
                                                           title = "Viewpoint",
                                                           value = textOutput("meta_vp"),
                                                           showcase = icon("eye"),
                                                           theme = 'green',
                                                           showcase_layout = "top right",
                                                           full_screen = TRUE,
                                                       )),
                                                column(6,
                                                       value_box(
                                                           title = "First and 2nd restiction enzyme",
                                                           value = textOutput("meta_re"),
                                                           showcase = icon("scissors"),
                                                           theme = 'green',
                                                           showcase_layout = "top right",
                                                           margin = margin(t = 10, b = 10, l = 10,
                                                                           r = 10),
                                                       ))),
                                       div(style = "height: 20px;"),
                                       fluidRow(column(6,
                                                       value_box(
                                                           title = "Organism",
                                                           value = textOutput("meta_orga"),
                                                           showcase = uiOutput("orga_icon"),
                                                           theme = "green",
                                                           showcase_layout = "top right",
                                                           full_screen = TRUE
                                                       )),
                                                column(6,
                                                       value_box(
                                                           title = "Read length",
                                                           value = textOutput("meta_rl"),
                                                           showcase = icon("ruler"),
                                                           theme = 'green',
                                                           showcase_layout = "top right",
                                                           full_screen = TRUE
                                                       ),),
                                       ),
                                       div(style = "height: 20px;"),
                                       fluidRow(column(6,
                                                       value_box(
                                                           title = "Condition",
                                                           subtitle = textOutput("meta_rep_co"),
                                                           value = textOutput("meta_cond_1"),
                                                           showcase = icon("capsules"),
                                                           theme = "green",
                                                           showcase_layout = "top right",
                                                           full_screen = TRUE
                                                       ),),
                                                column(6,
                                                       value_box(
                                                           title = "Control",
                                                           subtitle = textOutput("meta_rep_ct"),
                                                           value = textOutput("meta_ctrl_2"),
                                                           showcase = icon("x"),
                                                           theme = "green",
                                                           showcase_layout = "top right",
                                                           full_screen = TRUE
                                                       ))),
                                       hr(),
                                       h2(icon ("chart-simple"), "FASTQC statistics"),
                                       h3("Sequencing statistics"),
                                       dataTableOutput("fastqc_tab"),
                                       column(6,
                                              h3("Sequence counts"),
                                              withLoader(plotOutput("seq_plot"), loader = "dnaspin"),
                                              loaderOptions = list(color = "#2D973F")),
                                       column(6,
                                              h3("Per sequence quality scores"),
                                              withLoader(plotOutput("seq_qual"), loader = "dnaspin")),
                                       hr(),
                                       h2(icon("bars-staggered"), "Alignment statistics"),
                                       dataTableOutput("align_tab"),
                                       hr(),
                                       h2(icon("arrow-trend-up"), "4C-seq statistics"),
                                       dataTableOutput("basic_tab")
                                ),

                            ),
                            hr(),
                            "More quality control information can be found ",
                            "under results/[project_name]/qc and results/[project_name]/mutiqc_data. The ",
                            "full MultiQC report should be located at ",
                            "Datasets/[project_name]/multiqc_report.html."
                        ),
                        conditionalPanel(
                            condition = "!output.file_ready",
                            fluidRow(
                                style = "text-align: center; margin-top: 50px;",
                                p(
                                    style = "color: black;",
                                    "Please upload the necessary files to perform the analyses."
                                )
                            ))),

                ## Base Tools ---------------------------
                tabItem(tabName = "base_tools",
                        conditionalPanel(
                            condition = "output.file_ready",
                            h1(icon("dice-one"), "Base tool results",
                               bsButton(inputId = "bt",
                                        label = "",
                                        icon = icon("info"),
                                        size = "small")),
                            bsPopover(id = "bt",
                                      title = "Base tools",
                                      content = paste0("Our ensemble algorithm relies on four base tools: PeakC, r3Cseq, FourSig, and 4C-ker. Each tool is utilized with varying window size parameters, unless otherwise specified, in which case the default window size is employed. To visualize replicate data, we overlay plots of individual replicates on top of one another."),
                                      placement = "right",
                                      trigger = "hover",
                                      options = list(container = "body")
                            ),
                            hr(),
                            h2("Base tools karyoplot", bsButton(inputId = "kpbase",
                                                                label = "",
                                                                icon = icon("info"),
                                                                size = "small")),
                            bsPopover(id = "kpbase",
                                      title = "Karyoplot",
                                      content = paste0(
                                          "In the karyoplot, the tracks depict the smoothed read counts surrounding the viewpoint area. The rectangles below the karyoplot indicate whether a particular tool identifies the corresponding genomic region as an interacting region. Note that the karyoplot and the track plot below may exhibit differences due to variations in resolution: whereas the upper karyoplot displays smoothed read counts, the lower track plot presents unsmoothed data at a higher resolution."
                                      ),
                                      placement = "right",
                                      trigger = "hover",
                                      options = list(container = "body")
                            ),
                            fluidRow(
                                column(6, withLoader(plotOutput("karyo_base", height = 900),
                                                     loader = "dnaspin")),
                                conditionalPanel('output.karyo_base_ready',
                                                 column(6,
                                                        box(title = 'Karyoplot settings',
                                                            solidHeader = TRUE,
                                                            status = "primary",
                                                            width = 12,
                                                            h3("Select genes of interest"),
                                                            p("Please be aware that plotting genes can take some time..."),
                                                            withLoader(uiOutput("genes_select_base"),
                                                                       loader = "dnaspin"),
                                                            hr(),
                                                            fluidRow(
                                                                column(1),
                                                                column(10, uiOutput("in_regions_ui")),
                                                                column(1)
                                                            ),
                                                            br(),
                                                            fluidRow(
                                                                column(12,
                                                                       actionButton("button_genes_base", "Update karyoplot", icon = icon("reload"))
                                                                )))),
                                )),
                            h2("Base tool track plot"),
                            sliderInput("ymax_base", "Enter ymax.",
                                        min = 0,
                                        max = 40000,
                                        value = 20000,
                                        step = 1000,
                                        ticks = TRUE,
                                        width = "50%"),

                            withLoader(plotOutput("plot_bt", height = 800), loader = "dnaspin"),
                            h2('Interaction calls per tool'),
                            h3('Condition'),
                            dataTableOutput("tab_base_cond"),
                            h3('Control'),
                            dataTableOutput("tab_base_ctrl"),
                            h2("Compare condition and control"),
                            box(title = "UpsetR",
                                solidHeader = TRUE,
                                width = 12,
                                status = 'primary',
                                uiOutput("rb_upset_ui"),
                                hr(),
                                checkboxGroupInput("cb_tools", label = "Tools to compare",
                                                   choices = c("foursig_1",
                                                               "foursig_3",
                                                               "foursig_5",
                                                               "foursig_11",
                                                               "r3c_2000",
                                                               "r3c_5000",
                                                               "r3c_10000",
                                                               "peakc_11",
                                                               "peakc_21",
                                                               "peakc_31",
                                                               "peakc_51",
                                                               "r4cker_nearbait"),
                                                   inline = T),
                            ),
                            plotOutput("upset_base", height = 1200),
                            hr(),
                            #),
                            conditionalPanel(
                                condition = "!output.file_ready",
                                fluidRow(
                                    style = "text-align: center; margin-top: 50px;",
                                    p(
                                        style = "color: black;",
                                        "Please upload the necessary files to perform the analyses."
                                    )
                                )))),

                ## Ensemble -----------------------------
                tabItem(tabName = "ensemble",
                        conditionalPanel(
                            condition = "output.file_ready",
                            h1(icon("dice-four"), "Ensemble interaction analysis"),
                            hr(),
                            column(6,
                                   box(title = "Ensemble interaction calling settings",
                                       solidHeader = TRUE,
                                       status = "primary",
                                       width = 12,
                                       radioButtons("rb_model",
                                                    "Choose F1 or AUPRC optimized model.",
                                                    choices = c("F1", "AUPRC")),
                                       bsButton(inputId = "pu",
                                                label = "",
                                                icon = icon("info"),
                                                size = "small"),
                                       bsPopover(id = "pu",
                                                 title = "Model",
                                                 content = paste0(
                                                     "Weighted voting is either optimized on F1 score or Area Under Precision-Recall Curve."
                                                 ),
                                                 placement = "right",
                                                 trigger = "hover",
                                                 options = list(container = "body")
                                       ),
                                       actionButton("run_ens", "Start interaction calling.", icon = icon("play"))),
                            ),
                            column(6,
                                   conditionalPanel("output.ens_clicked",
                                                    box(title = 'Karyoplot settings',
                                                        solidHeader = TRUE,
                                                        status = "primary",
                                                        width = 12,
                                                        conditionalPanel("output.karyo_ens_ready",
                                                                         h3("Select genes of interest"), width = 12,
                                                                         p("Please be aware that plotting genes can take some time..."),
                                                                         withLoader(uiOutput("genes_select_ens"), loader = "dnaspin"),
                                                                         checkboxInput("spider_ens", "Show spiderplot.", FALSE),

                                                                         fluidRow(
                                                                             column(1, ''),
                                                                             column(10, uiOutput("in_regions_ens_ui")),
                                                                             column(1)
                                                                         ),
                                                                         br(),
                                                                         fluidRow(
                                                                             column(12,
                                                                                    actionButton("button_genes_ens", "Update karyoplot", icon = icon("reload"))
                                                                             )

                                                                             #textInput("in_regions_ens", "Please enter genomic regions to highlight.", "chr:start-end"),

                                                                         ),
                                                        )
                                                    )
                                   )),
                            conditionalPanel("output.ens_clicked",
                                             h2("Karyoplot fourSynergy", bsButton(inputId = "kpens",
                                                                                  label = "",
                                                                                  icon = icon("info"),
                                                                                  size = "small")),
                                             bsPopover(id = "kpens",
                                                       title = "Karyoplot",
                                                       content = paste0(
                                                           "In the karyoplot, the tracks depict the smoothed read counts surrounding the viewpoint area. The rectangles below the karyoplot indicate whether a fourSynergy identifies the corresponding genomic region as an interacting region. Note that the karyoplot and the track plot below may exhibit differences due to variations in resolution: whereas the upper karyoplot displays smoothed read counts, the lower track plot presents unsmoothed data at a higher resolution."
                                                       ),
                                                       placement = "right",
                                                       trigger = "hover",
                                                       options = list(container = "body")
                                             ),
                                             withLoader(plotOutput("karyo_ens", height = 1200), loader = "dnaspin"),
                                             h2("Track plots ensemble"),
                                             sliderInput("ymax_ens", "Enter ymax.",
                                                         min = 0,
                                                         max = 40000,
                                                         value = 20000,
                                                         step = 5000,     # Smaller steps for flexibility
                                                         ticks = TRUE,    # Auto-shows nice ticks
                                                         width = "50%",
                                                         pre = " "),

                                             withLoader(plotOutput("plot_bt_ens", height = 800), loader = "dnaspin"),
                                             h2('Interaction calls fourSynergy'),
                                             fluidRow(
                                                 column(6, h3('Condition'),
                                                        dataTableOutput("tab_ens_cond")),
                                                 column(6,
                                                        h3('Control'),
                                                        dataTableOutput("tab_ens_ctrl"))),
                                             hr(),
                                             h2(icon("code-compare"), "Compare conditions"),
                                             plotOutput("upset_ens"),
                                             selectInput("area_venn", " from Venn diagram to see which regions are included:",
                                                         choices = c("condition only",
                                                                     "overlap",
                                                                     "control only")),
                                             dataTableOutput("tab_ens")),
                            #),
                            conditionalPanel(
                                condition = "!output.file_ready",
                                fluidRow(
                                    style = "text-align: center; margin-top: 50px;",
                                    p(
                                        style = "color: black;",
                                        "Please upload the necessary files to perform the analyses."
                                    )
                                )))),

                ## Differential Analysis ----
                tabItem(tabName = "diff_analysis",
                        conditionalPanel(
                            condition = "output.file_ready",
                            h3(icon("dna"), "Differential interaction analysis"),
                            hr(),
                            fluidRow(
                                column(6,
                                       value_box(
                                           title = "Condition",
                                           subtitle = textOutput("meta_rep_co"),
                                           value = textOutput("meta_cond_2"),
                                           showcase = icon("capsules"),
                                           theme = "green",
                                           showcase_layout = "top right",
                                           full_screen = TRUE
                                       ),
                                       br(),

                                       value_box(
                                           title = "Control",
                                           subtitle = textOutput("meta_rep_ct"),
                                           value = textOutput("meta_ctrl_1"),
                                           showcase = icon("x"),
                                           theme = "green",
                                           showcase_layout = "top right",
                                           full_screen = TRUE
                                       )),
                                column(6,
                                       box(title = 'Karyoplot settings',
                                           width = 12,
                                           solidHeader = TRUE,
                                           status = "primary",
                                           conditionalPanel("output.karyo_diff_ready",
                                                            p("Please be aware that plotting genes can take some time..."),
                                                            uiOutput('genes_select_diff'),
                                                            column(1.5,),
                                                            column(11, uiOutput("in_regions_diff_ui")),
                                                            column(1.5,),
                                                            checkboxInput("spider_diff", "Show spiderplot.", FALSE),
                                                            actionButton("button_genes_diff", "Update karyoplot", icon = icon("reload")))))),
                            h2("Differential interaction karyoplot", bsButton(inputId = "kpdiff",
                                                                              label = "",
                                                                              icon = icon("info"),
                                                                              size = "small")),
                            bsPopover(id = "kpdiff",
                                      title = "Karyoplot",
                                      content = paste0(
                                          "In the karyoplot, the tracks depict the smoothed read counts surrounding the viewpoint area. The rectangles below the karyoplot indicate whether fourSynergy identifies the corresponding genomic region as an interacting region."
                                      ),
                                      placement = "right",
                                      trigger = "hover",
                                      options = list(container = "body")
                            ),
                            withLoader(plotOutput('karyo_diff', height = 800), loader = "dnaspin"),
                            fluidRow(
                                column(6,
                                       h2("Heatmap"),
                                       plotOutput("hm_diff", height = 800)),
                                column(6,
                                       h2("MA"),
                                       plotOutput("plot_ma", height = 800))),
                            h2("Results differential analysis (DESeq2)"),
                            dataTableOutput("tab_diff"),
                        ),
                        conditionalPanel(
                            condition = "!output.file_ready",
                            fluidRow(
                                style = "text-align: center; margin-top: 50px;",
                                p(
                                    style = "color: black;",
                                    "Please upload the required files to perform the analyses."
                                )
                            ))),
                tabItem("pipeline",
                        #fluidRow(
                         h2(icon("timeline"), "Pipeline"),
                        h3("Base tools used:"),
                        "r3Cseq (https://doi.org/10.1093/nar/gkt373), fourSig (https://doi.org/10.1093/nar/gku156), peakC (https://doi.org/10.1093/nar/gky443) and R.4Cker (https://doi.org/10.1371/journal.pcbi.1004780)",
                        "The pipeline is written in snakemake. You can find comprehensive documentation here: https://snakemake.github.io/",
                        hr(),
                        h2(icon("comments"), "Support"),
                        "If you have questions, feature requests or encounter any bugs, please contact Sophie Wind (sophie.wind@uni-muenster.de).\n",
                        hr(),
                        h2(icon("code"), "Code availibilty"),
                        "You can find the Bioconductor R package fourSynergy here: https://bioconductor.org/packages/fourSynergy <br/>",
                        "The fourSynergy R package is available here: https://github.com/sophiewind/fourSynergy <br/>",
                        "The fourSynergy snakemake pipeline can be found here: https://github.com/sophiewind/fourSynergy_pip",
                        hr(),
                        h2(icon("file-lines"), "Documentation"),
                        "A comprehensive description of fourSynergy is in the Vignette (https://bioconductor.org/packages/devel/bioc/vignettes/fourSynergy/inst/doc/fourSynergy_vignette.html).<br/>",
                        "The fourSynergy publication can be found here: https://github.com/sophiewind/fourSynergy_pip",
                        hr(),
                        h2(icon("bug"), "Bugs"),
                        "Bugs can be reported here: https://github.com/sophiewind/fourSynergy_pip/issues",
                        hr(),
                        h2(icon("balance-scale"), "Legal notice"),
                        p(tags$a(href = "https://www.medizin.uni-muenster.de/imi/impressum.html",
                                "Imprint / Legal notice (University of Münster)",
                                target = "_blank", rel = "noopener"))
                )
            )
        )
    )
)
# Definition des Servers
server <- function(input, output, session) {
    options(shiny.maxRequestSize = 2500 * 1024^2)

    # Upload logic -------------------------------------------------
    destfile_paths <- reactiveVal()
    trackfile_paths <- reactiveVal()
    config_path <- reactiveVal()
    multiqc_paths <- reactiveVal()
    config <- reactiveVal()
    karyo_base_ready <- reactiveVal()
    karyo_ens_ready <- reactiveVal()
    karyo_diff_ready <- reactiveVal()
    ens_clicked <- reactiveValues(clicked = FALSE)
    basic_paths <- reactiveVal()
    files_ready <- reactiveVal(0)
    config_valid <- reactiveVal(0)
    highlight_regions <- reactiveVal(NULL)
    error_values <- reactiveValues(error_message = NULL)

    ## Region io ####
    output$in_regions_ui <- renderUI({
        req(sia())
        tagList(
            tags$style(type = 'text/css',
                       ".irs-bar,.irs-handle,.irs-bar-edge,.irs-single,.irs-min,.irs-max {background: #2D973F;}"),
            tags$style(type = 'text/css', ".js-irs-0.irs-bar,.js-irs-0.irs-handle,.js-irs-0.irs-bar-edge,.js-irs-0.irs-single {background: #2D973F;}"),
            tags$style(HTML(".js-irs-2 .irs-single {background: #2D973F}")),
            div(id = 'big_slider',
                sliderInput("in_regions",
                            label = 'Select area of interest',
                            min = start(sia()@vfl[1]),
                            max = end(sia()@vfl[length(sia()@vfl)]),
                            value = c(start(sia()@vp), end(sia()@vp))
                )
            )
        )
    })

    output$in_regions <- renderUI({
        req(sia())
        tagList(
            tags$style(type = 'text/css',
                       ".irs-bar,.irs-handle,.irs-bar-edge,.irs-single {background: #2D973F;}"),
            tags$style(type = 'text/css', ".js-irs-0.irs-bar,.js-irs-0.irs-handle,.js-irs-0.irs-bar-edge,.js-irs-0.irs-single {background: #2D973F;}"),
            div(id = 'big_slider',
                sliderInput("in_region", min = start(sia()@vfl[1]), max = end(sia()@vfl[length(sia()@vfl)]), 100)
            )
        )
    })

    output$in_regions_ens_ui <- renderUI({
        req(sia())
        tagList(
            tags$style(type = 'text/css',
                       ".irs-bar,.irs-handle,.irs-bar-edge,.irs-single,.irs-min,.irs-max {background: #2D973F;}"),
            tags$style(type = 'text/css', ".js-irs-0.irs-bar,.js-irs-0.irs-handle,.js-irs-0.irs-bar-edge,.js-irs-0.irs-single {background: #2D973F;}"),
            div(id = 'big_slider',
                sliderInput("in_regions_ens",
                            label = 'Select area of interest',
                            min = start(sia()@vfl[1]),
                            max = end(sia()@vfl[length(sia()@vfl)]),
                            value = c(start(sia()@vp), end(sia()@vp))
                )
            )
        )
    })

    output$in_regions_ens <- renderUI({
        req(sia())
        tagList(
            tags$style(type = 'text/css',
                       ".irs-bar,.irs-handle,.irs-bar-edge,.irs-single {background: #2D973F;}"),
            tags$style(type = 'text/css', ".js-irs-0.irs-bar,.js-irs-0.irs-handle,.js-irs-0.irs-bar-edge,.js-irs-0.irs-single {background: #2D973F;}"),
            div(id = 'big_slider',
                sliderInput("in_regions_ens", min = start(sia()@vfl[1]), max = end(sia()@vfl[length(sia()@vfl)]), 100)
            )
        )
    })

    output$in_regions_diff_ui <- renderUI({
        req(sia())
        tagList(
            tags$style(type = 'text/css',
                       ".irs-bar,.irs-handle,.irs-bar-edge,.irs-single,.irs-min,.irs-max {background: #2D973F;}"),
            tags$style(type = 'text/css', ".js-irs-0.irs-bar,.js-irs-0.irs-handle,.js-irs-0.irs-bar-edge,.js-irs-0.irs-single {background: #2D973F;}"),
            div(id = 'big_slider',
                sliderInput("in_regions_diff",
                            label = 'Select area of interest',
                            min = start(sia()@vfl[1]),
                            max = end(sia()@vfl[length(sia()@vfl)]),
                            value = c(start(sia()@vp), end(sia()@vp))
                )
            )
        )
    })

    output$in_regions_diff <- renderUI({
        req(sia())
        tagList(
            tags$style(type = 'text/css',
                       ".irs-bar,.irs-handle,.irs-bar-edge,.irs-single {background: #2D973F;}"),
            tags$style(type = 'text/css', ".js-irs-0.irs-bar,.js-irs-0.irs-handle,.js-irs-0.irs-bar-edge,.js-irs-0.irs-single {background: #2D973F;}"),
            div(id = 'big_slider',
                sliderInput("in_regions_diff", min = start(sia()@vfl[1]), max = end(sia()@vfl[length(sia()@vfl)]), 100)
            )
        )
    })


    ## Config check ####
    observe({
        req(input$config)
        cg <- input$config
        cg <- cg %>% mutate(dir = paste0(dirname(datapath), '/', name))

        # Increment files_ready
        isolate({
            files_ready(files_ready() + 1)
        })

        tryCatch(
            expr = {
                if (checkConfig(cg$datapath)) {
                    message("Check ok")
                    config_valid(1)
                    cat("config_valid after:", config_valid())
                    config(read_yaml(cg$datapath))
                    uploads_dir <- file.path(getwd(), paste0("/Datasets/", config()$author, "/"))
                    dir.create(uploads_dir, recursive = T)
                    print(paste("dir for uploaded config:", uploads_dir))
                    destfile <- file.path(uploads_dir, cg$name)
                    config_path(destfile)
                    file.copy(cg$datapath, destfile, overwrite = TRUE)
                } else {
                    message("Check failed")
                    config_valid(-1)
                }
            },
            error = function(e) {
                config_valid(-1)
                error_values$error_message <- e$message
            }
        )
    })

    output$error_message_config <- renderUI({
        req(error_values$error_message)
        tagList(
            tags$p(style = "color: red;",
                   paste("Error:", error_values$error_message)),
        )

    })

    ## Sia definition ####
    output$sia_info <- renderUI({
        req(input$shiny_in)
        sias <- input$shiny_in
        sias <- sias[endsWith(sias$name, ".bed"),]
        sias <- sias %>% mutate(dir = paste0(dirname(datapath), '/', name))
        req_sias <- paste0(config()$author, "_", rep(c("foursig_1",
                                                       "foursig_3",
                                                       "foursig_5",
                                                       "foursig_11",
                                                       "r3c_2000",
                                                       "r3c_5000",
                                                       "r3c_10000",
                                                       "peakcSig_11",
                                                       "peakcSig_21",
                                                       "peakcSig_31",
                                                       "peakcSig_51",
                                                       "r4cker_nearbait")), '_',
                           rep(c('condition'), each = 12), '_nearbait.bed')

        if (!is.null(config()$author)){
            req_sias <- c(req_sias, paste0(config()$author, "_", rep(c("foursig_1",
                                                                       "foursig_3",
                                                                       "foursig_5",
                                                                       "foursig_11",
                                                                       "r3c_2000",
                                                                       "r3c_5000",
                                                                       "r3c_10000",
                                                                       "peakcSig_11",
                                                                       "peakcSig_21",
                                                                       "peakcSig_31",
                                                                       "peakcSig_51",
                                                                       "r4cker_nearbait")), '_',
                                           rep(c('control'), each = 12), '_nearbait.bed'))
        }

        missing_sia <- req_sias[!req_sias %in% sias$name]

        if (length(missing_sia) > 0) {
            tagList(
                tags$p(style = "color: red;",
                       "The following interaction files are missing:"),
                tags$ul(
                    lapply(missing_sia, function(x) tags$li(x))
                ),
                tags$p(style = "color: black;",
                       "You can look for those in ./results/[Dataset]/sias.",
                       "The files are critical important for the analysis.",
                       "If they are not produced ",
                       "please rerun postprocessing and Basic4Cseq."),
            )
        } else {
            tags$p(style = "color: green;", "All paths are available.")

            # Increment files_ready
            isolate({
                files_ready(files_ready() + 1)
            })

            uploads_dir <- file.path(getwd(),  paste0("/results/", config()$author, "/sia/"))
            dir.create(uploads_dir, recursive = T)
            print(paste("Verzeichnis für hochgeladene Dateien:", uploads_dir))

            destfile <- file.path(uploads_dir, sias$name)

            destfile_paths(paste0(getwd(), "/results/", config()$author, "/"))

            for (i in seq_along(sias$name)) {
                print(sias$name[i])
                if (endsWith(sias$name[i], 'nearbait_area.bed')){
                    out <- paste0(getwd(), paste0("/results/", config()$author, "/nearbait_area.bed"))
                    file.copy(sias$datapath[i], out, overwrite = TRUE)
                    print(gsub(pattern = 'sia', '', destfile[i]))
                } else {
                    file.copy(sias$datapath[i], destfile[i], overwrite = TRUE)
                }
            }
        }
    })

    # Upload files ####
    output$track_info <- renderUI({
        req(input$shiny_in)

        # Tracks in
        tp <- input$shiny_in
        tp <-  tp[grepl("\\.(bam|bai|bedGraph)$", tp$name), ]
        tp <- tp %>% mutate(dir = paste0(dirname(datapath), '/', name))

        # Check if all files there
        if (!is.null(config()$control)){
            req_align <- paste0(rep(paste0(c(config()$condition, config()$control), "_", 
                                           rep(seq(1, max(config()$conditionRep)), each = 2)), each = 3),
                                c("_sorted.bam","_sorted.bam.bai", "_sorted.bedGraph"))
        } else {
            req_align <- paste0(rep(paste0(paste0(config()$condition, "_"), seq(1,max(config()$conditionRep))), each = 3),  
                                c("_sorted.bam","_sorted.bam.bai", "_sorted.bedGraph"))
        }

        missing_a <- req_align[!req_align %in% tp$name]

        if (length(missing_a) > 0) {
            tagList(
                tags$p(style = "color: red;",
                       "The following files are missing:"),
                tags$ul(
                    lapply(missing_a, function(x) tags$li(x)),
                ),
                tags$p(style = "color: black;",
                       "You can look for those in ./results/[Dataset]/alignment."),
            )
        } else {
            tags$p(style = "color: green;", "All paths are available.")

            # Increment files_ready
            isolate({
                files_ready(files_ready() + 1)
            })
            uploads_dir <- file.path(paste0(getwd(), "/results/", config()$author, "/alignment/"))
            dir.create(uploads_dir, recursive = T)
            print(paste("Verzeichnis für upload files:", uploads_dir))
            destfile <- file.path(uploads_dir, tp$name)
            trackfile_paths(paste0(getwd(), "/results/", config()$author, "/alignment/"))
            for (i in seq_along(tp$name[endsWith(tp$name, c("bam", "bam.bai", "bedGraph"))])) {
                print(destfile[i])
                file.copy(tp$datapath[i], destfile[i], overwrite = TRUE)
            }
            message(paste0("### ", paste0(getwd(), "/results/", config()$author, "/")))
        }
    })

    output$basic_info <- renderUI({
        req(input$shiny_in)

        # Tracks in
        bp <- input$shiny_in
        bp <- bp[endsWith(bp$name, "_stats.txt"),]
        bp <- bp %>% mutate(dir = paste0(dirname(datapath), '/', name))

        # Check if all files there
        if (!is.null(config()$control)){
            req_align <- paste0(rep(paste0(c(config()$condition, config()$control), "_",
                                           rep(seq(1, max(config()$conditionRep)), each = 2))),
                                "_stats.txt")
        } else {
            req_align <- paste0(rep(paste0(paste0(config()$condition, "_"), seq(1,max(config()$conditionRep)))),
                                "_stats.txt")
        }

        missing_a <- req_align[!req_align %in% bp$name]

        if (length(missing_a) > 0) {
            tagList(
                tags$p(style = "color: red;",
                       "The following files are missing:"),
                tags$ul(
                    lapply(missing_a, function(x) tags$li(x))
                ),
                tags$p(style = "color: black;",
                       "You can look for those in ",
                       "./results/[Dataset]/basic4cseq. The files are not ",
                       "required for the analysis."),
            )
        } else {
            tags$p(style = "color: green;", "All paths are available.")

            # Increment files_ready
            isolate({
                files_ready(files_ready() + 1)
            })

            uploads_dir <- file.path(paste0(
                getwd(), "/results/", config()$author, "/basic4cseq/stats/"))
            dir.create(uploads_dir, recursive = T)
            print(paste("dir for uploaded files:", uploads_dir))
            destfile <- file.path(uploads_dir, bp$name)
            basic_paths(paste0(getwd(), "/results/", config()$author,
                               "/basic4cseq/stats/"))
            for (i in seq_along(bp$name[endsWith(bp$name, "stats.txt")])) {
                print(destfile[i])
                file.copy(bp$datapath[i], destfile[i], overwrite = TRUE)
            }
        }
    })

    output$multiqc_info <- renderUI({
        req(input$shiny_in)
        mq <- input$shiny_in
        mq <- mq %>% mutate(dir = paste0(dirname(datapath), '/', name))

        # Create dir
        uploads_dir <- file.path(getwd(), "results/", config()$author, "/multiqc_data")
        print(paste("Dir for upload files:", uploads_dir))
        dir.create(uploads_dir, recursive = T)

        # Check if all files there
        req_mq <- c('multiqc_fastqc.txt',
                    'fastqc_sequence_counts_plot.txt',
                    'fastqc_per_sequence_quality_scores_plot.txt',
                    'samtools-flagstat-pct-table.txt')

        mq <- mq[mq$name %in% req_mq,]
        missing_mq <- req_mq[!req_mq %in% mq$name]

        if (length(missing_mq) > 0) {
            tagList(tags$p(style = "color: red;",
                           "The following QC files are missing:"),
                    tags$ul(
                        lapply(missing_mq, function(x) tags$li(x))
                    ),
                    tags$p("You can look for those in ./results/[Dataset]/qc ",
                           "or ./results/[Dataset]/multiqc_data. Those files are not ",
                           "required for the analysis."),
            )

        } else {
            tags$p(style = "color: green;", "All paths are available.")
            print(paste("Name of uploaded file:", mq$name))
        }

        destfile <- file.path(uploads_dir, mq$name)
        multiqc_paths(paste0(getwd(), "/results/", config()$author, "/multiqc_data/"))
        for (i in seq_along(mq$name)) {
            print(destfile[i])
            file.copy(mq$datapath[i], destfile[i], overwrite = TRUE)
        }
        message(paste0("### ", getwd(),
                       "/results/", config()$author, "/multiqc_data/"))
    })

    # Display conditional panels if file is uploaded successfully
    output$file_ready <- reactive({
        if (files_ready() >= 4) TRUE
    })

    output$config_status <- reactive({
        config_valid()
    })
    outputOptions(output, "file_ready", suspendWhenHidden = FALSE)
    outputOptions(output, "config_status", suspendWhenHidden = FALSE)

    # General validation function
    is_valid_region <- function(region_input) {
        grepl("chr[1-9XYM]+:\\d+-\\d+|chr:start-end", region_input)
    }

    # Create sia ----
    sia <- reactive({
        x <- destfile_paths()
        y <- config_path()
        z <- trackfile_paths()
        sia <- createIa(destfile_paths(),
                        config_path(),
                        trackfile_paths())
        output$genes_select_base <- get_gene_selection("base", sia)
        output$genes_select_ens <- get_gene_selection("ens", sia)
        output$genes_select_diff <- get_gene_selection("diff", sia)
        sia
    })

    observeEvent(input$run_ens, {
        ens_clicked$clicked <- TRUE
    })

    # Keep the rest of your code the same
    sia_ens <- eventReactive(input$run_ens, {
        consensusIa(sia(), model = input$rb_model)
    })

    # Display conditional panels if karyo is plotted
    output$ens_clicked <- reactive({
        req(sia_ens)
        ens_clicked$clicked
    })

    outputOptions(output, "ens_clicked", suspendWhenHidden = FALSE)

    sia_diff <-  eventReactive(input$run_ens,{
        differentialAnalysis(sia_ens())})

    # Metadata ----
    output$orga_icon <- renderUI({
        if (sia()@metadata$organism == "mm10") {
            icon("paw")
        } else {
            icon("person")
        }
    })

    output$sample_name <- renderText({
        paste0("Sample: ", sia()@metadata$author)
    })

    output$meta_re <- renderText({
        paste(sia()@metadata$REEnz, collapse = ", ")
    })

    output$meta_vp <- renderText({
        paste0("chr:", sia()@metadata$VPchr, "\n", start(sia()@vp))
    })
    output$meta_rl <- renderText({
        sia()@metadata$readLength
    })
    output$meta_cond_1 <- renderText({
        sia()@metadata$condition
    })

    output$meta_cond_2 <- renderText({
        sia()@metadata$condition
    })

    output$meta_rep_co <- renderText({
        max(sia()@metadata$conditionRep)
    })

    output$meta_orga <- renderText({
        sia()@metadata$organism
    })
    output$meta_ctrl_1 <- renderText({
        sia()@metadata$control
    })

    output$meta_ctrl_2 <- renderText({
        sia()@metadata$control
    })
    output$meta_rep_ct <- renderText({
        max(sia()@metadata$controlRep)
    })


    # Upload ----
    output$upload_text <- renderText({
        paste("You selected:", input$datei)
    })

    output$upload_text <- renderText({
        paste("You selected:", input$datei)
    })

    output$metadata <- renderDataTable({
        file.out <- paste0('fourSyerngy_metadata_', sia()@metadata$author)
        sia()@metadata %>%
            unlist() %>%
            as.matrix(ncol = 2) %>%
            as.data.frame() %>%
            datatable(colnames = "",
                      extensions = "Buttons",
                      options = list(pageLength = 21,
                                     dom = "Bfrtip",
                                     buttons = list(
                                         list(extend = "csv", filename = file.out),
                                         list(extend = "excel", filename = file.out),
                                         list(extend = "copy", filename = file.out))))
    }, server = FALSE)

    # QC ----
    output$fastqc_tab <- renderDataTable({
        req(files_ready() >= 2)
        file.out <- paste0('fourSyerngy_fastqc_', sia()@metadata$author)

        # Quality table
        read.delim(paste0(multiqc_paths(), "/multiqc_fastqc.txt")) %>%
            DT::datatable(., extensions = "Buttons",
                          options = list(
                              dom = "Bfrtip",
                              scrollX = TRUE,
                              buttons = list(
                                  list(extend = "csv", filename = file.out),
                                  list(extend = "excel", filename = file.out),
                                  list(extend = "copy", filename = file.out)
                              )))
    }, server = FALSE)

    output$seq_plot <- renderPlot({
        #browser()
        p <- read.delim(paste0(multiqc_paths(), "/fastqc_sequence_counts_plot.txt")) %>%
            melt() %>%
            ggplot(., aes(x = value, y = Sample, fill = variable)) +
            geom_bar(stat = 'identity') +
            theme(panel.background = element_blank(),
                  # panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text = element_text(size = 11),
                  axis.title = element_text(size = 11))
        (p)
    })

    output$seq_qual <- renderPlot({
        mc <- multiqc_paths()
        p <- read.delim(paste0(multiqc_paths(), "fastqc_per_sequence_quality_scores_plot.txt")) %>%
            pivot_longer(cols = starts_with("X"), names_to = "time", values_to = "val") %>%
            mutate(
                val = gsub("\\(|\\)", "", val),
            ) %>%
            separate(val, c('x', 'y'), sep = ', ') %>%
            mutate(x = as.numeric(x),
                   y = as.numeric(y)) %>%
            dplyr::select(Sample, time, x, y) %>%
            ggplot(., aes(x = x, y = y, color = Sample)) +
            geom_rect(aes(xmin = -Inf, xmax = 20, ymin = -Inf, ymax = Inf),
                      fill = "#E7C2C2", alpha = 0.03, color = NA) +
            geom_rect(aes(xmin = 20, xmax = 28.5, ymin = -Inf, ymax = Inf),
                      fill = "#E9DEC2", alpha = 0.03, color = NA) +
            geom_rect(aes(xmin = 28.5, xmax = Inf, ymin = -Inf, ymax = Inf),
                      fill = "#C2E5C2", alpha = 0.03, color = NA) +
            geom_line() +
            xlim(c(0, 40)) +
            labs(x = "Mean sequence quality (Phred score)",
                 y = "Count") +
            theme(axis.text = element_text(size = 11),
                  axis.title = element_text(size = 11))
        (p)
    })


    output$align_tab <- renderDataTable({
        file.out <- paste0('fourSyerngy_flagstat_', sia()@metadata$author)

        read.delim(paste0(multiqc_paths(), '/samtools-flagstat-pct-table.txt')) %>%
            datatable(colnames = paste0(gsub('\\.', ' ', colnames(.)), " (%)"),
                      rownames = FALSE,
                      extensions = "Buttons",
                      options = list(
                          dom = "Bfrtip",
                          buttons = list(
                              list(extend = "csv", filename = file.out),
                              list(extend = "excel", filename = file.out),
                              list(extend = "copy", filename = file.out)
                          )))
    }, server = FALSE)

    output$basic_tab <- renderDataTable({
        sia_val <- sia()
        file.out <- paste0('fourSyerngy_basic4cseq_', sia_val@metadata$author)
        p <- basic_paths()
        files <- list.files(basic_paths(),
                            pattern = ".*stats.txt", full.names = TRUE)
        dfs <- lapply(files, read.delim, sep = ":", header = FALSE)

        coll <- Reduce(rbind, dfs)
        coll$Sample <- rep(gsub("_stats.txt", "", basename(files)), each = 4)
        coll <- coll %>%
            pivot_wider(id_cols = V1, names_from = Sample, values_from = V2)

        colnames(coll)[1] <- ''
        df <- coll %>%
            t() %>%
            as.data.frame()

        colnames(df) <- df[1, ]

        df <- df[-1, ]

        datatable(df, extensions = "Buttons",
                  options = list(
                      dom = "Bfrtip",
                      buttons = list(
                          list(extend = "csv", filename = file.out),
                          list(extend = "excel", filename = file.out),
                          list(extend = "copy", filename = file.out)
                      )))
    }, server = FALSE)

    regions <- reactiveVal(data.frame(start = numeric(), end = numeric()))
    observeEvent(input$add_region, {
        req(input$highlight_on)
        new <- data.frame(
            start = input$in_regions[1],
            end = input$in_regions[2]
        )
        regions(rbind(regions(), new))
    })
    observeEvent(input$clear_regions, {
        regions(data.frame(start = integer(), end = integer()))
    })

    highlight_regions <- reactive({
        req(input$highlight_on)
        if(nrow(regions()) == 0 || is.null(sia()@vp)) return(NULL)
        regions_clean <- na.omit(regions())
        paste(
            paste0(seqname(sia()@vp), ":", regions_clean$start, "-", regions_clean$end),
            collapse = ", "
        )
    })

    # Base Tools ####
    output$base_tools_text <- renderText({
        ""
    })

    # Display conditional panels if karyo is plotted
    output$karyo_base_ready <- reactive({
        TRUE
    })

    outputOptions(output, "karyo_base_ready", suspendWhenHidden = FALSE)

    # Create a reactive value to track button clicks
    clicked <- reactiveValues(button_clicked = FALSE)

    # Define a reactive trigger
    plot_trigger <- eventReactive(input$button_genes_base, {
        if (is.null(input$in_regions)){
            reg <- paste0('chr', sia()@metadata$VPchr,':', sia()@metadata$VPpos,'-', sia()@metadata$VPpos)
        } else {
            reg <- paste0('chr', sia()@metadata$VPchr,':', input$in_regions[1],'-', input$in_regions[2])
        }
        # Capture current state
        list(
            sia_val = sia(),
            genes = if(is.null(input$genes_select_base)) NULL else input$genes_select_base,
            regions = reg#input$in_regions
        )
    }, ignoreNULL = FALSE)

    output$karyo_base <- renderPlot({
        message('Update karyo')
        req(plot_trigger())  # Wait for trigger
        current <- plot_trigger()
        plotIaIndiviualTools(
            current$sia_val,
            cex.chr = 2, cex.ideo = 1, cex.leg = 1,
            cex.y.track = 1, cex.y.lab = 1, cex.vp = 1.5,
            genes_of_interest = current$genes, highlight_regions = current$regions
        )
    })


    output$plot_bt <- renderPlot({
        req(sia())
        req(plot_trigger())
        req(input$ymax_base)

        sia_val <- sia()
        current <- plot_trigger()
        if (is.null(input$ymax_base)){
            max <- 20000
        } else {
            max <- input$ymax_base
        }
        p <- plotBaseTracks(sia_val, highlight_regions = current$regions, max_range = max)  #
        plot_grid(p[[1]] +
                      theme(axis.text.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            axis.title.x = element_blank(),
                            text = element_text(size = 16)),
                  p[[2]] +
                      theme(text = element_text(size = 16)), ncol = 1)

    })

    # Condition
    output$tab_base_cond <- renderDataTable({
        sia_val <- sia()
        file.out <- paste0('fourSyerngy_base_algorithm_condition_',
                           sia_val@metadata$author)

        vfl <- sia_val@vfl
        for (i in seq(1, length(sia_val@expInteractions))){
            ov <- findOverlaps(vfl, sia_val@expInteractions[[i]])
            if (!is.na(sia_val@expInteractions[[i]]$tool[1] )){
                mcols(vfl)[[sia_val@expInteractions[[i]]$tool[1]]] <- '-'
                mcols(vfl[queryHits(ov),])[sia_val@expInteractions[[i]]$tool[1]] <- 'X'
            } else {
                mcols(vfl)[[names(sia_val@expInteractions[i])]] <- '-'
            }

        }

        vfl.df <- vfl %>% as.data.frame(keep.extra.columns = TRUE)# %>%
        colnames(vfl.df) <- gsub('000', 'k', gsub('Sig', '', gsub('_condition', '', gsub('rep.', '', colnames(vfl.df)))))

        vfl.df %>%
            dplyr::select(-width, -strand) %>%
            filter(if_any(everything(.), ~. == "X")) %>%
            datatable(rownames = FALSE,
                      extensions = "Buttons",
                      options = list(
                          dom = "Bfrtip",
                          buttons = list(
                              list(extend = "csv", filename = file.out),
                              list(extend = "excel", filename = file.out),
                              list(extend = "copy", filename = file.out)
                          )))
    }, server = FALSE)

    output$tab_base_ctrl <- renderDataTable({
        sia_val <- sia()
        file.out <- paste0('fourSyerngy_base_algorithm_control_', sia_val@metadata$author)
        vfl <- sia_val@vfl

        if (!is.null(sia()@metadata$control)){
            for (i in seq(1, length(sia_val@ctrlInteractions))){
                ov <- findOverlaps(vfl, sia_val@ctrlInteractions[[i]])
                if (!is.na(sia_val@ctrlInteractions[[i]]$tool[1] )){
                    mcols(vfl)[[sia_val@ctrlInteractions[[i]]$tool[1]]] <- '-'
                    mcols(vfl[queryHits(ov),])[sia_val@ctrlInteractions[[i]]$tool[1]] <- 'X'
                } else {
                    mcols(vfl)[[names(sia_val@ctrlInteractions[i])]] <- '-'
                }
            }

            vfl.df <- vfl %>% as.data.frame(keep.extra.columns = TRUE)
            colnames(vfl.df) <- gsub('000', 'k',
                                     gsub('Sig', '',
                                          gsub('_control', '',
                                               gsub('rep.', '', colnames(vfl.df)))))


            vfl.df %>%
                dplyr::select(-width, -strand) %>%
                filter(if_any(everything(.), ~. == "X")) %>%
                datatable(rownames = FALSE,
                          extensions = "Buttons",
                          options = list(
                              dom = "Bfrtip",
                              buttons = list(
                                  list(extend = "csv", filename = file.out),
                                  list(extend = "excel", filename = file.out),
                                  list(extend = "copy", filename = file.out)
                              )))
        } else {
            data.frame('Seqnames' = character(), 'Start' = integer(),
                       'End' = integer(), 'Rating' = factor()) %>%
                datatable(rownames = FALSE,
                          extensions = "Buttons",
                          options = list(
                              dom = "Bfrtip",
                              buttons = list(
                                  list(extend = "csv", filename = file.out),
                                  list(extend = "excel", filename = file.out),
                                  list(extend = "copy", filename = file.out)
                              )))
        }

    }, server = FALSE)

    ## Upset ---
    output$rb_upset_ui <- renderUI({
        choices <- if (!is.null(sia()@metadata$control)){
            c("Tools in condition" = "tico", "Tools in control" = "tict",
              "Condition and control" = "coco")
        } else {
            c("Tools in condition" = "tico")
        }

        radioButtons("rb_upset", "Select what you want to compare",
                     choices = choices, selected = "tico")
    })


    observeEvent(input$rb_upset, {
        updateCheckboxGroupInput(session, "cb_tools", selected = NULL)
    })

    output$upset_base <- renderPlot({
        sia_val <- sia()
        if (input$rb_upset == "tico"){
            sel <- input$cb_tools
            if (length(sel) <= 2){
                ggplot(data.frame(), aes()) +
                    theme_void() +
                    labs(title = "Please select at least three tools to create an upset plot.")
            }
            else{
                sia_val@expInteractions[names(sia_val@expInteractions) %in% paste0('rep.', sel, '.condition')] %>%
                    as.list() %>%
                    lapply(function(x) as.data.frame(x, stringsAsFactors = FALSE)) %>%
                    lapply(function(x) x %>%
                               dplyr::mutate(region = paste0(seqnames, ":", start, "-", end))) %>%
                    lapply(function(x) x %>%
                               dplyr::select(region) %>%
                               unlist) %>%
                    fromList(.) %>%
                    upset(.,order.by = "freq", text.scale = 2, nsets = length(sel))
            }
        } else if (input$rb_upset == "tict"){
            sel <- input$cb_tools
            if (length(sel) <= 2){
                ggplot(data.frame(), aes()) +
                    theme_void() +
                    labs(title =  "Please select at least three tools to create an upset plot.")
            }
            else{
                sia_val@ctrlInteractions[names(sia_val@ctrlInteractions) %in% paste0('rep.', sel, '.control')] %>%
                    as.list() %>%
                    lapply(function(x) as.data.frame(x, stringsAsFactors = FALSE)) %>%
                    lapply(function(x) x %>%
                               dplyr::mutate(region = paste0(seqnames, ":", start, "-", end))) %>%
                    lapply(function(x) x %>%
                               dplyr::select(region) %>%
                               unlist) %>%
                    fromList(.) %>%
                    upset(.,order.by = "freq", text.scale = 2, nsets = length(sel))
            }
        } else {
            sel <- input$cb_tools

            if (length(sel) <= 1){
                ggplot(data.frame(), aes()) +
                    theme_void() +
                    labs(title =  "Please select at least two tools to create an upset plot.")
            }
            else{
                sel <- input$cb_tools
                comb <- c(sia_val@expInteractions, sia_val@ctrlInteractions)
                data <- comb[names(comb) %in% c(paste0("rep.", sel, ".condition"), paste0("rep.", sel, ".control"))] %>%
                    as.list() %>%
                    lapply(function(x) as.data.frame(x, stringsAsFactors = FALSE)) %>%
                    lapply(function(x) x %>%
                               dplyr::mutate(region = paste0(seqnames, ":", start, "-", end))) %>%
                    lapply(function(x) x %>%
                               dplyr::select(region) %>%
                               unlist) %>%
                    fromList(.)
                font <- ifelse(ncol(data) >= 8, 1, 2)
                upset(data, order.by = "freq", nsets = length(sel)*2, text.scale = font)
            }
        }
    })

    # Ensemble ####
    # Define a reactive trigger
    plot_trigger_ens <- eventReactive(input$button_genes_ens, {
        if (is.null(input$in_regions_ens)){
            reg <- paste0('chr', sia_ens()@metadata$VPchr,':', sia_ens()@metadata$VPpos,
                          '-', sia_ens()@metadata$VPpos)
        } else {
            reg <- paste0('chr', sia_ens()@metadata$VPchr,':', input$in_regions_ens[1],'-',
                          input$in_regions_ens[2])
        }
        # Capture current state
        list(
            sia_val = sia_ens(),
            genes = if(is.null(input$genes_select_ens)) NULL else input$genes_select_ens,
            regions = reg,
            spider = input$spider_ens
        )
    }, ignoreNULL = FALSE)

    output$karyo_ens <- renderPlot({
        message('Update karyo')
        req(plot_trigger_ens())
        current <- plot_trigger_ens()
        plotConsensusIa(current$sia_val, cex.chr = 2, cex.ideo = 1,
                        cex.leg = 1, cex.y.track = 1, cex.vp = 1.5, cex.y.lab = 1,
                        genes_of_interest = current$genes,
                        highlight_regions = current$regions,
                        plot_spider = current$spider
        )
    })
    output$karyo_ens_ready <- reactive({
        TRUE
    })

    outputOptions(output, "karyo_ens_ready", suspendWhenHidden = FALSE)
    output$plot_bt_ens <- renderPlot({
        req(sia_ens())
        req(plot_trigger_ens())
        req(input$ymax_ens)
        current <- plot_trigger_ens()
        if (is.null(input$ymax_ens)){
            max <- 20000
        } else {
            max <- input$ymax_ens
        }
        sia_val <- sia_ens()
        message('regions')
        p <- plotConsensusTracks(sia_val, highlight_regions = current$regions, max_range = max)
        plot_grid(p[[1]] +
                      theme(text = element_text(size = 16)),
                  p[[2]] +
                      theme(text = element_text(size = 16)), ncol = 2)
    })

    output$tab_ens_cond <- renderDataTable({
        file.out <- paste0('fourSyerngy_ensemble_condition_', sia_ens()@metadata$author)
        sia_ens()@expConsensus %>%
            as.data.frame(keep.extra.columns = TRUE) %>%
            filter(significance > 0) %>%
            dplyr::select(seqnames, start, end, rate_total_condition) %>%
            `colnames<-`(c('Seqnames', 'Start', 'End', 'Rating')) %>%
            datatable(rownames = FALSE, extensions = "Buttons",
                      options = list(
                          dom = "Bfrtip",
                          buttons = list(
                              list(extend = "csv", filename = file.out),
                              list(extend = "excel", filename = file.out),
                              list(extend = "copy", filename = file.out)
                          )))
    }, server = FALSE)

    output$tab_ens_ctrl <- renderDataTable({
        file.out <- paste0('fourSyerngy_ensemble_control_', sia_ens()@metadata$author)
        sia_ens()@ctrlConsensus %>%
            as.data.frame(keep.extra.columns = TRUE) %>%
            filter(significance > 0) %>%
            dplyr::select(seqnames, start, end, rate_total_control) %>%
            `colnames<-`(c('Seqnames', 'Start', 'End', 'Rating')) %>%
            datatable(rownames = FALSE, extensions = "Buttons",
                      options = list(
                          dom = "Bfrtip",
                          buttons = list(
                              list(extend = "csv", filename = file.out),
                              list(extend = "excel", filename = file.out),
                              list(extend = "copy", filename = file.out)
                          )))
    }, server = FALSE)

    output$upset_ens <- renderPlot({
        gVenn::computeOverlaps(GRangesList(Condition = sia_ens()@expConsensus[sia_ens()@expConsensus$significance > 0], Control = sia_ens()@ctrlConsensus[sia_ens()@ctrlConsensus$significance > 0])) %>%
            plotVenn(fontsize = 11, fills = list(fill = c('firebrick4', 'darkblue'), alpha = 0.7))
    })

    output$tab_ens <- renderDataTable({
        file.out <- paste0('fourSyerngy_ensemble_', input$area_venn,
                           '_', sia_ens()@metadata$author)
        ov <- gVenn::computeOverlaps(GRangesList(Condition = sia_ens()@expConsensus[sia_ens()@expConsensus$significance > 0], Control = sia_ens()@ctrlConsensus[sia_ens()@ctrlConsensus$significance > 0]))
        if (input$area_venn == "condition only"){
            set <- ov$reduced_regions[ov$reduced_regions$intersect_category == '10']
        } else if (input$area_venn == 'overlap'){
            set <- ov$reduced_regions[ov$reduced_regions$intersect_category == '11']
        } else {
            set <- ov$reduced_regions[ov$reduced_regions$intersect_category == '01']
        }
        set %>%
            as.data.frame() %>%
            dplyr::select(-width, -strand, -intersect_category) %>%
            `colnames<-`(c('Seqnames', 'Start', 'End')) %>%
            datatable(., extensions = "Buttons",
                      rownames = FALSE,
                      options = list(
                          dom = "Bfrtip",
                          buttons = list(
                              list(extend = "csv", filename = file.out),
                              list(extend = "excel", filename = file.out),
                              list(extend = "copy", filename = file.out)
                          )))
    }, server = FALSE)

    # Diff ----
    # Display conditional panels if karyo is plotted
    output$karyo_diff_ready <- reactive({
        TRUE
    })
    outputOptions(output, "karyo_diff_ready", suspendWhenHidden = FALSE)

    ## Karyo ensemble ####
    # Define a reactive trigger
    plot_trigger_diff <- eventReactive(input$button_genes_diff, {
        if (is.null(input$in_regions_diff)){
            reg <- paste0('chr', sia_diff()@metadata$VPchr,':', sia_diff()@metadata$VPpos,
                          '-', sia_diff()@metadata$VPpos)
        } else {
            reg <- paste0('chr', sia_diff()@metadata$VPchr,':', input$in_regions_diff[1],'-',
                          input$in_regions_diff[2])
        }
        # Capture current state
        list(
            sia_val = sia_diff(),
            genes = if(is.null(input$genes_select_diff)) NULL else input$genes_select_diff,
            regions = reg,
            spider = input$spider_diff
        )
    }, ignoreNULL = FALSE)

    output$karyo_diff <- renderPlot({
        message('Update karyo')
        req(plot_trigger_diff())
        current <- plot_trigger_diff()
        plotDiffIa(current$sia_val, cex.leg = 1, cex.y.lab = 1,
                   genes_of_interest = current$genes,
                   highlight_regions = current$regions,
                   plot_spider = current$spider, cex.y.track = 1)
    })

    output$tab_diff <- renderDataTable({
        file.out <- paste0('fourSyerngy_differential_', sia_diff()@metadata$author)
        sia_diff()@differential %>%
            as.data.frame() %>%
            datatable(., extensions = "Buttons",
                      options = list(
                          dom = "Bfrtip",
                          buttons = list(
                              list(extend = "csv", filename = file.out),
                              list(extend = "excel", filename = file.out),
                              list(extend = "copy", filename = file.out)
                          )))
    }, server = FALSE)

    output$plot_ma <- renderPlot({
        sia_diff()@differential %>%
            plotMA()
    })

    output$hm_diff <- renderPlot({
        sia_diff()@dds %>%
            counts() %>%
            heatmap
    })

}

# Starten der App
shinyApp(ui = ui, server = server)
