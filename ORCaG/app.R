#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
#Author:Pavlo_hrab
#Build as a part of HackBio internship
#version -> 0.5


#Import required libraries
library(shiny)
library(tidyverse)
library(plyr)
library(ggplot2)
library(hrbrthemes)
library(forcats)
library(ggtree) #Bioconductor library
library(viridis)
library(ape)
library(epiR)
library(tidytree)
library(gridExtra)
library(formattable)
library(shinycssloaders)
library(grid)
library(ggplotify)
library(knitr)
library(ggpubr)
library(ggsignif)
library(phylocanvas)
library(shinyjs)
library(getopt) #Installed via Bioconductor
library(optparse)
library(plotrix)
library(readODS) 
library(plotly)
library(epitools)
#If there are import problems, feel free to install packages via install.packages() or via Bioconductor


# Define UI for application
ui <- fluidPage(

    # Application title
    titlePanel("ORCaG"),

    # Sidebar with a several inputs 
    sidebarLayout(
        sidebarPanel(
          #File upload input
            fileInput("gisaid", 
                      "Upload Gisaid dataset:"),
            #File upload input
            fileInput("status", 
                      "Upload status assignment csv file:"),
            #File upload input
            fileInput("tree", 
                      "Upload newick tree:"),
            #Multiple select input for tree labels choices
            selectInput("tree_labels", 
                        "Tree labels are:", list("Virus.name column"=0,"Accession.ID column"=1),
                        selected = 0),
            #Multiple select input for itol node coloring
            selectInput("itol_background", 
                        "Choose iTOL column for clade coloring in iTOL:", list("Country"="Country", "Host"="Host", "Gender"="Gender", "Patient.status"="Patient.status", "Patient.age"="Patient.age" , "Specimen"="Specimen" ,"Lineage"="Lineage", "Clade"="Clade" ),
                        selected = "Patient.status"),
            #Button click to generate iTOL files
            actionButton("itol","Generate iTOL file and annotations" ),
            #Button to download generated iTOL files
            downloadButton("itol_download","Download iTOL files" ),
            #Link to iTOL website for tree uploading
            helpText(a("Go to iTOL upload", href = "https://itol.embl.de/upload.cgi", target="_blank")),
            #Link to our project tree
            helpText(a("Our example iTOL tree", href = "https://itol.embl.de/tree/1172201664490531598537092", target="_blank")),
            #Multiple select input for Hospitalized patient status
            selectInput("status_choice", 
                        "Hospitalized patient status:", list("Default"=0,"Low risk"=1,"High risk"=2, "Remove from data"=3),
                        selected = 0),
            #Check box input for what to do with alive/live/symptomatic patients
            selectInput("more_data", 
                        "Include alive/live/symptomatic patients in the analysis:", list("Default"=0,"Include/Low risk"=1,
                                                                                         "Remove from data"=2),
                        selected = 0),
            #Button to download currently used dataset
            downloadButton("current_dataset","Download currently used dataset"),
            #Slider to choose age cutoff
            sliderInput("age_cutoff",
                        "Age cutoff:",
                        min = 1,
                        max = 100,
                        value=50)
        ),

        # Show a plot of the generated distribution
        mainPanel(
          plotlyOutput("oddsPlot"),
          #Texts for OR for specific age cutoff view
          textOutput("age_cut"),
          textOutput("age_mes"),
          textOutput("age_min"),
          textOutput("age_max"),
          #Tree view
          plotOutput("tree"),
          plotOutput("ageDensity"),
          plotOutput("Countries"),
          plotlyOutput("country_count"),
          plotlyOutput("country_odds"),
          plotOutput("boxplot_age"),
          plotOutput("cladeAge"),
          #Interactive tree view with phylocabvas
          phylocanvasOutput("tree_interactive")
        )
    )
)

#Server logic to draw texts and charts
server <- function(input, output) {
  options(warn = -1)
  #Main function to generate iTOL files. From https://github.com/mgoeker/table2itol 
    create_itol_files <- function(infiles, identifier = "ID", label = "Label",
                                background = "", identifier2 = "", directory = ".", colour.file = "",
                                gradient.file = "", separator = "\t", na.strings = paste0(c("", "(null)",
                                                                                            "NA"), collapse = separator), quote = "\"", abort = FALSE,
                                conversion = "none", double.to.bars = FALSE, emblems = "", template = "%s",
                                max.size = 20L, favour = 1, width = 0.5, precision = 1L, restrict = "",
                                opacity = 1) {
    
    
    OLDOPT <- options(warn = 1L)
    on.exit(options(OLDOPT))
    
    
    # EL  ellipse
    # RE  rectangle
    # TL  left pointing triangle
    # TR  right pointing triangle
    # DI  rhombus (diamond)
    # HH  horizontal hexagon
    # HV  vertical hexagon
    # PL  left pointing pentagram
    # PR  right pointing pentagram
    # PU  up pointing pentagram
    # PD  down pointing pentagram
    # OC  octagon
    # GP  rectangle (gap; black filled rectangle with 1/3 normal height)
    #
    SYMBOLS <- c("EL", "RE", "TL", "TR", "DI", "HH", "HV",
                 "PL", "PR", "PU", "PD", "OC", "GP")
    
    
    # Restricted set: 1 = square, 2 = circle, 3 = star, 4 = right triangle,
    # 5 = left triangle.
    #
    BRANCH_SYMBOLS <- seq_len(5L)
    
    
    # Like branch symbols. We omit #6 (check mark) because it does not display
    # nicely.
    #
    BINARY_SYMBOLS <- seq_len(5L)
    
    
    BLACK <- "#000000"
    
    
    LIGHTGREY <- "#E5E5E5" # R's grey90
    
    
    WHITE <- "#FFFFFF"
    
    
    OUTPUT_SEPARATOR <- "\t"
    
    
    # Used as end points of colour gradients, with white as other end point,
    # and for colouring binary data (logical vectors)
    #
    SPECIAL_COLORS <- c("#1f78b4", "#e31a1c", "#33a02c", "#b15928",
                        "#6a3d9a", "#ff7f00")
    
    
    # Colour vectors collected by Jan P. Meier-Kolthoff.
    #
    COLOURS <- list(
      # Dark2; colour-blind-safe
      JMK01 = "#1b9e77",
      # Dark2; colour-blind-safe
      JMK02 = c("#d95f02", "#1b9e77"),
      # Dark2; colour-blind-safe
      JMK03 = c("#1b9e77", "#d95f02", "#7570b3"),
      # 4-class Paired; colour-blind-safe
      JMK04 = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"),
      # 5-class Accent; print-friendly
      JMK05 = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                "#fb9a99"),
      # 6-class Paired; print-friendly
      JMK06 = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                "#fb9a99", "#e31a1c"),
      # 7-class Paired; print-friendly
      JMK07 = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                "#fb9a99", "#e31a1c", "#fdbf6f"),
      # Dark2; print-friendly
      JMK08 = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a",
                "#66a61e", "#e6ab02", "#a6761d", "#666666"),
      # 9-class Set1; print-friendly
      JMK09 = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
                "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999"),
      # 10-class Paired
      JMK10 = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"),
      # 11-class Paired
      JMK11 = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
                "#ffff99"),
      # 12-class Paired
      JMK12 = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
                "#ffff99", "#b15928"),
      ## from here on: iwanthue (all colours, hard)
      JMK13 = c("#8393c7", "#8ad256", "#6a49c5", "#d2b351",
                "#cb55c3", "#4d4040", "#c4527c", "#57743d", "#d85439", "#7accb1",
                "#925136", "#ceb2ab", "#512f67"),
      JMK14 = c("#a2d1cd", "#5d39a8", "#71d14c", "#cb56c7",
                "#7ed094", "#4d4040", "#7077b8", "#c28b4c", "#cd9dae", "#c64a34",
                "#55868c", "#cccb51", "#b2436e", "#567137"),
      JMK15 = c("#92d4ad", "#6842c1", "#6ecf58", "#cb4ec2",
                "#55733d", "#4d4040", "#c99447", "#9083cb", "#c9d14f", "#4d2c63",
                "#cea4a2", "#d54f38", "#71a6bd", "#ca507f", "#823f33"),
      JMK16 = c("#76a5bd", "#bfdf44", "#cf4bab", "#66c95b",
                "#7c42c5", "#4d4040", "#7279ca", "#c27837", "#4b2a62", "#c7b956",
                "#cc8cb5", "#536e3b", "#d74746", "#84d3ae", "#893b42", "#cdb19a"),
      JMK17 = c("#823f35", "#77d952", "#6d44c4", "#78d5a1",
                "#cf4a70", "#4d4040", "#ca53bd", "#69923c", "#6d7fc4", "#d1d04e",
                "#532b63", "#d64d31", "#4b623d", "#ca96b7", "#78b5c2", "#ccbf9b",
                "#c58741"),
      JMK18 = c("#697bc5", "#5e9742", "#6641c0", "#7bdc57",
                "#c954c9", "#4d4040", "#4d2b62", "#73d6ac", "#d6493d", "#75adbe",
                "#c54883", "#526339", "#caca9b", "#7b332e", "#cfcf49", "#c89dc8",
                "#c58738", "#c78980"),
      JMK19 = c("#9e693f", "#9147d5", "#c9d747", "#9482d3",
                "#61913d", "#4d4040", "#6dd85e", "#d049a4", "#76d0b6", "#d5493c",
                "#6897bb", "#d7993d", "#553291", "#c7cb8a", "#472f5b", "#cd7993",
                "#496340", "#ccb8bc", "#7f2c3a"),
      JMK20 = c("#7295c1", "#d44b38", "#6ad14f", "#6a3bc0",
                "#cedb44", "#4d4040", "#77d192", "#cb4fc3", "#b1b85f", "#7772cc",
                "#d9973b", "#4f2b62", "#79d1cf", "#cc497b", "#4a6c2e", "#c990b5",
                "#752e30", "#d1c5ac", "#a26f47", "#537e71"),
      JMK21 = c("#90b5d9", "#d6532d", "#c84ccc", "#74d147",
                "#512d79", "#4d4040", "#6740c8", "#cace49", "#6b79d1", "#6ccc84",
                "#c8478c", "#74c4b8", "#cc4458", "#4f6781", "#cb9142", "#552443",
                "#c6cb97", "#82442d", "#c489c5", "#546d37", "#cb9896"),
      JMK22 = c("#392c51", "#4d4040", "#642c79", "#792d3b",
                "#6a3ec6", "#875b30", "#4f7231", "#547f72", "#d24637", "#6d71ce",
                "#d2497e", "#cd4fc8", "#6a8fbc", "#d88742", "#c78dc6", "#cc9795",
                "#c7af40", "#68cd55", "#72d4a6", "#9ecfd6", "#c9cb8f", "#c3de48"),
      JMK23 = c("#8ad93f", "#c749c4", "#5e8f3d", "#6639be",
                "#73d979", "#4d4040", "#d4ca4a", "#6c6ccc", "#d78c3b", "#6485b9",
                "#d24635", "#70d4ae", "#cc4279", "#cbcb99", "#4c295f", "#ce867e",
                "#793130", "#84cbd7", "#896c35", "#c27bbb", "#364e27", "#cab2cb",
                "#5b837b"),
      JMK24 = c("#ccc79a", "#6a42c7", "#d0a540", "#cc49c9",
                "#6dd755", "#4d4040", "#de5a26", "#7cc7d0", "#cc3f47", "#78d8a5",
                "#5e2d78", "#c9da51", "#6679d0", "#bf7348", "#c6b7d8", "#5f903c",
                "#c47ec5", "#6a5b29", "#ce4684", "#497359", "#772d38", "#c3858c",
                "#352444", "#5b7a9e"),
      JMK25 = c("#6ba43c", "#c74ace", "#cbe14b", "#6847cd",
                "#6ede53", "#4d4040", "#cbb248", "#592e82", "#d6842f", "#5e78c1",
                "#76dd99", "#c6438e", "#4b8047", "#cf4c67", "#7acdc4", "#d2472f",
                "#7ba5c4", "#79322f", "#c388cf", "#78662f", "#45294d", "#c8cd9d",
                "#3e5d4a", "#d08c6c", "#c698a9"),
      JMK26 = c("#73d991", "#b44adb", "#71d94d", "#cf4cb4",
                "#ccde4d", "#4d4040", "#ceae44", "#5a41c2", "#cdd09c", "#652e7a",
                "#83d7ce", "#dc4338", "#536e83", "#d34a79", "#5d9073", "#c68dc7",
                "#619339", "#85b1d7", "#da8340", "#6978cb", "#9d4533", "#34284e",
                "#d09e9e", "#732d41", "#364e25", "#866a38"),
      JMK27 = c("#363258", "#6ed853", "#5b3fc7", "#c9de43",
                "#b54ad9", "#4d4040", "#5c2c7e", "#b7d17b", "#cf4a83", "#6ed9a4",
                "#cd4450", "#8fd3d5", "#d74527", "#769ac1", "#d27d3f", "#6d75cf",
                "#d4af42", "#4f8c3b", "#d14eba", "#568778", "#c692c8", "#344625",
                "#d4c7a6", "#722e4c", "#c88988", "#7a3a25", "#86783a"),
      JMK28 = c("#7f3a27", "#71da53", "#c14bd4", "#55933d",
                "#626ad0", "#4d4040", "#623ac4", "#cbd943", "#542c79", "#c1d483",
                "#bc7fd0", "#6ad7a3", "#d84330", "#71bec7", "#ce7537", "#6f99d8",
                "#d5aa43", "#546586", "#7c7233", "#ce429f", "#3e6344", "#ce7d9f",
                "#2d1d38", "#c6b3ce", "#793151", "#bfcbae", "#d24566", "#c8927d"),
      JMK29 = c("#cdc2c2", "#663dc8", "#76dd51", "#c64ece",
                "#cfda49", "#4d4040", "#549e3f", "#7577da", "#d3522e", "#7cd6ce",
                "#d4425b", "#77de9a", "#542a7e", "#d1d395", "#321e3d", "#d74a98",
                "#95963d", "#586095", "#db9a3e", "#77abd9", "#8b3c67", "#639575",
                "#d08982", "#456129", "#ca92cc", "#896134", "#597984", "#742c28",
                "#283a28"),
      JMK30 = c("#31223c", "#bbe141", "#c94edb", "#65d559",
                "#8b3899", "#4d4040", "#613ec8", "#df9b36", "#6e75d5", "#c16c39",
                "#402a74", "#cfc248", "#da47a4", "#63d6ad", "#d94330", "#6abccd",
                "#c58181", "#617fae", "#7f2f2c", "#b5cfb8", "#833b65", "#b5d888",
                "#cc88cb", "#4e8a3b", "#d6466a", "#476d58", "#d2b284", "#544320",
                "#c9b6d0", "#867c36"),
      JMK31 = c("#913d83", "#ced242", "#6643d0", "#79d949",
                "#c249d4", "#4d4040", "#db45a4", "#68dc88", "#3a1f4f", "#c3d483",
                "#532e8e", "#da983e", "#6d79d5", "#9b4b29", "#d085d5", "#8b7d3b",
                "#c9a0c0", "#54913d", "#dc4b32", "#72d4b1", "#8f3e58", "#90d0d8",
                "#592720", "#d2c7a9", "#21262c", "#d64769", "#3b4f25", "#6ea2cf",
                "#cd887a", "#5c6089", "#568477"),
      JMK32 = c("#8f8b38", "#663cc8", "#6bd546", "#c74cce",
                "#b1d773", "#4d4040", "#c6e03a", "#59287c", "#5edb86", "#d14592",
                "#7ad9b1", "#da4627", "#719cd8", "#dc973a", "#6e71d7", "#dbc348",
                "#ca84c8", "#4c8b3a", "#d5445a", "#84ccd6", "#7f3353", "#d3c99f",
                "#2e1c38", "#ca7442", "#5a558b", "#803325", "#537286", "#cc8585",
                "#314826", "#cab3cc", "#7e6136", "#618d75"),
      JMK33 = c("#d64e9e", "#6cd54c", "#dd49d1", "#c8dd41",
                "#a152dd", "#4d4040", "#5139c2", "#ceaa3b", "#432d7c", "#c6d179",
                "#8f379a", "#70d68c", "#d9432f", "#6ad5be", "#d5416a", "#76c2d7",
                "#d87a71", "#6a75d5", "#836834", "#c988d1", "#598939", "#7a3260",
                "#bed3b3", "#8f372e", "#6082b3", "#d47c35", "#312749", "#d4ac8b",
                "#314825", "#cab9d7", "#4b211f", "#ad788b", "#568275"),
      JMK34 = c("#d8436c", "#653cc7", "#b4dc41", "#d143d0",
                "#5fd857", "#4d4040", "#a4db84", "#c64496", "#6adcad", "#de4830",
                "#6aa3d9", "#d98731", "#6271d1", "#dec841", "#b062cd", "#528e36",
                "#c28acd", "#675b2c", "#cbb7d3", "#a53332", "#528089", "#532878",
                "#d9d393", "#2a1e3c", "#8ed4d3", "#834629", "#5e5e8a", "#a08e3c",
                "#2b482a", "#d78763", "#619470", "#c87b8d", "#702944", "#c3a994"),
      JMK35 = c("#72d4cf", "#ccdf3e", "#5533c1", "#70d951",
                "#ac42d6", "#4d4040", "#6d66dc", "#b9c866", "#562a84", "#71da99",
                "#db43c7", "#518f39", "#d04497", "#314826", "#bc6cc9", "#5d8b74",
                "#d2416d", "#72abd3", "#dd461f", "#6078c6", "#d7ab3b", "#c49ad6",
                "#7d6b2f", "#cab8c4", "#3c1a20", "#c8ddb6", "#312652", "#cfb182",
                "#7c3463", "#c98271", "#576782", "#d24243", "#cb7a99", "#82372d",
                "#cf7734"),
      JMK36 = c("#6ade4b", "#6344d3", "#7bdc86", "#b746d4",
                "#65a234", "#4d4040", "#dbc941", "#552c93", "#bee148", "#dc3fb4",
                "#62d7b4", "#903a7e", "#4a8245", "#cf74d0", "#da993a", "#3e255f",
                "#c0d3b2", "#291d2d", "#cdce7e", "#752c41", "#7dcbd6", "#c43c44",
                "#669bcf", "#de4e28", "#5b5e83", "#c97449", "#bd92d0", "#847933",
                "#d7417a", "#558279", "#d07d92", "#364525", "#ceb9d0", "#763d23",
                "#6872d2", "#be9880"),
      JMK37 = c("#645b8e", "#80dc40", "#4f2ea4", "#69dc7b",
                "#d848cd", "#4d4040", "#8548da", "#c7d84e", "#96368e", "#afd995",
                "#d54227", "#61d9b9", "#db4187", "#4a9339", "#cd83d6", "#7a8431",
                "#6870d5", "#e3bc3b", "#6b9bd7", "#d87935", "#6fbfcf", "#cd3e50",
                "#c3d8c8", "#772e29", "#dbc38b", "#3f2267", "#bf9340", "#cab1d6",
                "#304726", "#b2918d", "#2a1f35", "#d5816f", "#5e8c6b", "#c77192",
                "#497080", "#7d592d", "#732d52"),
      JMK38 = c("#cf8ad0", "#74e042", "#b946da", "#5be080",
                "#5834c1", "#4d4040", "#d248bb", "#59a434", "#8064d4", "#b4dc4e",
                "#893876", "#96db99", "#d9478a", "#499052", "#627bcf", "#dfd238",
                "#47277a", "#908f39", "#79a2d8", "#d79234", "#4c7788", "#df502c",
                "#625984", "#d7d27b", "#2e1d3b", "#6bdac4", "#d34557", "#6a8b73",
                "#9e4427", "#cfb5cd", "#78562e", "#7cc6d5", "#26392b", "#cdcfb2",
                "#702735", "#bd7984", "#405924", "#d59571"),
      JMK39 = c("#8b308f", "#74dd41", "#6939ca", "#cce346",
                "#d545d2", "#4d4040", "#b271dd", "#e39b39", "#5050bc", "#cabc46",
                "#3a1f64", "#5cde7e", "#d9428e", "#57a56d", "#d63949", "#76dfc2",
                "#7e3052", "#b7e28f", "#d286c6", "#66a234", "#6d83d8", "#d65629",
                "#76c3d2", "#843326", "#6aa0d5", "#9c762c", "#5f5488", "#d48e70",
                "#4a6a81", "#d36778", "#466b2c", "#b28491", "#273825", "#c1b47a",
                "#301b31", "#d0d2bd", "#6c552d", "#c9b8d8", "#5f8675"),
      JMK40 = c("#3c2b5d", "#dee032", "#ab48d5", "#5bd749",
                "#db49c6", "#4d4040", "#5c42d0", "#a4e040", "#462687", "#d8b136",
                "#8d3989", "#60d076", "#d7468f", "#63d8b5", "#de4528", "#77c7d6",
                "#d13a55", "#5f8c7b", "#ce88d5", "#759b31", "#696ecd", "#de8739",
                "#6f9ad6", "#b75738", "#aadc90", "#946d89", "#d0dc6a", "#2c1a25",
                "#c6d8bc", "#782849", "#ceb977", "#283f27", "#d9798c", "#447c3d",
                "#ceb8d4", "#635b2d", "#c79783", "#733426", "#476682", "#98762e")
    )
    
    
    GENERATED_FILES <- new.env(TRUE, emptyenv())
    
    
    # Note that type.convert() does not recognize some spellings of false/true.
    # The value for TRUE must be the first one of each row.
    #
    BINARY_VALUE_SPELLINGS <- matrix(c(
      "y", "n",
      "t", "f",
      "true", "false",
      "yes", "no",
      "on", "off"
    ), 5L, 2L, TRUE)
    
    
    ## Helper functions
    
    
    # E.g. anyNA() is only available from 3.1.0 on.
    #
    check_R_version <- function(wanted = numeric_version("3.2.0")) {
      if (getRversion() < wanted)
        stop(sprintf("need a newer version of R, %s or higher", wanted))
      invisible(TRUE)
    }
    
    
    # Checking makes sense because colour vectors can be user-defined.
    #
    check_colour_vectors <- function(x, lengthcheck) {
      if (lengthcheck && !identical(seq_along(x), lengths(x, FALSE)))
        stop("incorrectly arranged colour vectors")
      bad <- !vapply(x, is.character, NA)
      if (any(bad))
        stop("wrong data type of colour vector(s) no. ",
             paste0(seq_along(x)[bad], collapse = ", "))
      bad <- vapply(x, anyDuplicated.default, 0L)
      if (any(bad))
        stop("duplicated value in colour vector(s) no. ",
             paste0(seq_along(x)[bad > 0L], collapse = ", "))
      invisible(TRUE)
    }
    
    
    # This also converts trivial names for colours into RGB codes.
    #
    standardize_colour <- function(x, opacity) {
      x <- col2rgb(x, TRUE)
      x["alpha", ] <- as.integer(x["alpha", ] * opacity)
      tolower(rgb(x["red", ], x["green", ], x["blue", ],
                  if (all(x["alpha", ] == 255L)) NULL else x["alpha", ], NULL, 255L))
    }
    
    
    # For input of user-defined colour vectors.
    #
    read_colour_vectors <- function(file, upto) {
      if (!nzchar(file, FALSE))
        return(NULL)
      x <- yaml::yaml.load_file(file)
      if (!is.list(x))
        x <- list(x)
      n <- lengths(x)
      x[n > 0L & n <= upto]
    }
    
    
    # Input method dispatch is based on file extension. Depends on extra library
    # for Excel and Libreoffice/Openoffice files, respectively.  Must ensure
    # character vectors are converted to factors.
    #
    read_file <- function(file, sep, na, quote) {
      
      read_xl <- function(sheet, path, na) {
        # for some reason read_excel() yields a 'tibble' instead of a data frame,
        # which does not display the same behaviour of `[`; hence we convert
        tryCatch(expr = as.data.frame(readxl::read_excel(path = path, na = na,
                                                         sheet = sheet, col_names = TRUE, col_types = NULL, skip = 0L,
                                                         trim_ws = FALSE)), error = function(e) {
                                                           warning(e) # a typical error is to encounter an empty sheet
                                                           data.frame() # now we can treat this later on ourselves
                                                         })
      }
      
      rescue_integers <- function(x) { # necessary for 'tibble' input
        is_whole_number <- function(x, tolerance = .Machine$double.eps ^ 0.5) {
          all(is.na(x) | abs(x - round(x)) < tolerance) # see ?is.integer
        }
        for (i in which(vapply(x, is.double, NA)))
          if (is_whole_number(x[, i]))
            storage.mode(x[, i]) <- "integer"
          x
      }
      
      rescue_factors <- function(x) { # not necessary for CSV input
        for (i in which(vapply(x, is.character, NA)))
          x[, i] <- factor(x[, i])
        x
      }
      
      switch(
        EXPR = tolower(tools::file_ext(file)),
        ods = lapply(lapply(lapply(X = readODS::ods_sheets(file), path = file,
                                   FUN = readODS::read_ods, na = na[[1L]], col_names = TRUE,
                                   col_types = NULL, formula_as_formula = FALSE, skip = 0L, range = NULL),
                            rescue_integers), rescue_factors),
        xls =,
        xlsx = lapply(lapply(lapply(readxl::excel_sheets(file),
                                    read_xl, file, na), rescue_integers), rescue_factors),
        list(read.table(file = file, header = TRUE, sep = sep, quote = quote,
                        na.strings = na, fill = FALSE, stringsAsFactors = TRUE, dec = ".",
                        check.names = FALSE, comment.char = ""))
      )
      
    }
    
    
    # Checking makes sense if 'chr' is used to join strings together.
    #
    assert_no_forbidden_character <- function(chr, ...) {
      x <- list(...)
      for (str in lapply(x[vapply(x, is.factor, NA)], levels.default)) {
        bad <- grepl(chr, str, FALSE, FALSE, TRUE)
        if (any(bad))
          stop(sprintf("string '%s' contains forbidden character '%s'",
                       str[bad][[1L]], chr))
      }
      invisible(TRUE)
    }
    
    
    # We add white at the end, assuming this represents NA, when NA values occur.
    #
    select_colours <- function(size, hasna) {
      if (hasna) {
        message(sprintf("Fetching %i colour(s) ...", size - 1L))
        if (size < 2L)
          WHITE
        else
          c(COLOURS[[size - 1L]], WHITE)
      } else {
        message(sprintf("Fetching %i colour(s) ...", size))
        if (size < 1L)
          character()
        else
          COLOURS[[size]]
      }
    }
    
    
    # Used for generating legend titles.
    #
    pretty_str <- function(x) {
      chartr("_.", "  ", x)
    }
    
    
    # We assume NA values have already been removed.
    #
    legend_range <- function(x, precision) {
      if (length(precision))
        sprintf(sprintf("%%s (%%.%if)", precision), c("Min.", "Max."), range(x))
      else
        sprintf("%s (%i)", c("Min.", "Max."), range(x))
    }
    
    
    # Used to not display branch symbols associated with certain values.
    #
    mask_if_requested <- function(x, cutoff, restriction) {
      outliers <- function(x, n) {
        me <- median(x, na.rm = TRUE)
        ma <- mad(x, na.rm = TRUE)
        x > me + ma * n | x < me - ma * n
      }
      if (is.na(cutoff))
        return(logical(length(x)))
      switch(
        EXPR = restriction,
        atleast = x < cutoff,
        beyond = !outliers(x, cutoff),
        larger = x <= cutoff,
        smaller = x >= cutoff,
        upto = x > cutoff,
        within = outliers(x, cutoff),
        stop(sprintf("unkown 'restriction' value '%s'", restriction))
      )
    }
    
    
    # Used to modify several vectors at once. Needs at least one argument.
    #
    coordinated_na_removal <- function(...) {
      ok <- !is.na(..1)
      if (all(ok))
        return(FALSE)
      args <- list(...)
      names(args) <- all.names(match.call(), FALSE, -1L, FALSE)
      parentframe <- parent.frame()
      for (name in unique.default(names(args)))
        assign(name, args[[name]][ok], parentframe)
      TRUE
    }
    
    
    # Helper function to convert output vectors or matrices.
    #
    convert_to <- function(x, unavailable)  {
      storage.mode(x) <- typeof(unavailable)
      x[is.na(x)] <- unavailable
      x
    }
    
    
    # Used for generating the output filename.
    #
    itol_filename <- function(colname, kind, directory) {
      result <- file.path(directory, sprintf("iTOL_%s-%s.txt", kind,
                                             gsub("\\W", "_", colname, FALSE, TRUE)))
      if (exists(result, GENERATED_FILES))
        stop(sprintf("name clash: file '%s' has already been generated", result))
      GENERATED_FILES[[result]] <- TRUE
      message(sprintf("Generating %s file for column '%s' ...", kind, colname))
      result
    }
    
    
    # Here '...' contains the data part.
    #
    print_itol <- function(outdir, title, annotation, ...) {
      
      join <- function(x) {
        if (!length(x))
          return(NULL)
        if (is.null(names(x)))
          stop("non-empty annotation lists must have names")
        if (!all(vapply(x, is.atomic, NA)))
          stop("non-empty annotation lists must contain only atomic values")
        x <- x[sort.list(names(x))]
        sizes <- lengths(x, FALSE)
        for (i in which(sizes > 1L))
          x[[i]] <- paste0(x[[i]], collapse = OUTPUT_SEPARATOR)
        for (i in which(!sizes))
          x[[i]] <- ""
        paste(names(x), unlist(x, FALSE, FALSE), sep = OUTPUT_SEPARATOR)
      }
      
      if (is.character(annotation)) {
        colname <- annotation
        annotation <- NULL
      } else {
        colname <- get("DATASET_LABEL", annotation)
      }
      
      kind <- switch(
        EXPR = title,
        branchsymbols = "DATASET_SYMBOL",
        collapse = "COLLAPSE",
        labels = "LABELS",
        treecolors = "TREE_COLORS",
        binary =,
        colorstrip =,
        domains =,
        gradient =,
        heatmap =,
        simplebar =,
        text = sprintf("DATASET_%s", toupper(title)),
        stop(sprintf("unknown title '%s'", title))
      )
      
      separator <- sprintf("SEPARATOR %s", switch(
        EXPR = OUTPUT_SEPARATOR,
        `\t` = "TAB",
        stop(sprintf("output separator '%s' not yet supported", OUTPUT_SEPARATOR))
      ))
      
      file <- itol_filename(colname, title, outdir)
      
      cat(kind, separator, join(annotation), "DATA", file = file,
          labels = NULL, sep = "\n", fill = FALSE, append = FALSE)
      cat(paste(..., sep = OUTPUT_SEPARATOR, collapse = NULL), file = file,
          labels = NULL, sep = "\n", fill = FALSE, append = TRUE)
      
    }
    
    
    ## Functions for special columns
    
    
    # For labelling the leaves.
    #
    emit_itol_labeltexts <- function(x, ids, name, outdir, ...) {
      coordinated_na_removal(x, ids)
      print_itol(outdir, "labels", name, ids, x)
    }
    
    
    # For colouring the leaves. 'x' is a factor or coerced to a factor, hence NAs
    # do not get removed.
    #
    emit_itol_labelcolors <- function(x, ids, name, outdir, ...) {
      if (!is.factor(x))
        x <- factor(x)
      x <- addNA(x, TRUE)
      size <- length(levels.default(x))
      if (size > length(COLOURS)) {
        warning(sprintf("skipping column '%s', which yields > %i levels",
                        name, length(COLOURS)))
        return()
      }
      annotation <- list(
        COLOR = "#a6cee3",
        DATASET_LABEL = name,
        LEGEND_COLORS = select_colours(size, anyNA(levels.default(x))),
        LEGEND_LABELS = levels.default(x),
        LEGEND_SHAPES = rep.int(1L, size),
        LEGEND_TITLE = pretty_str(name)
      )
      print_itol(outdir, "treecolors", annotation,
                 ids, "range", annotation$LEGEND_COLORS[x], x)
    }
    
    
    ## Functions for columns according to data type (class)
    
    
    # Output varies depending on the number of colours and symbols chosen and/or
    # available. 'x' is a factor, hence NAs do not get removed.
    #
    emit_itol_factor <- function(x, ids, name, outdir, symbols, maxclrs,
                                 favour, borwid, ...) {
      
      product <- function(x, y) {
        cbind(rep(x = x, each = length(y)), rep.int(y, length(x)))
      }
      
      x <- addNA(x, TRUE)
      annot1 <- list(DATASET_LABEL = name, MARGIN = 5, COLOR = "#bebada")
      size <- length(levels.default(x))
      
      if (size > maxclrs * length(SYMBOLS)) {
        
        # additional columns: position, colour, style, size_factor, rotation
        print_itol(outdir, "text", annot1, ids, x,
                   -1, BLACK, "normal", 0.75, 0)
        
      } else if (length(symbols) || size > maxclrs) {
        
        if (length(symbols)) {
          symbols <- vapply(split.default(symbols, x), `[[`, "", 1L)
          clrs <- select_colours(size, anyNA(levels.default(x)))
        } else {
          nsym <- ncls <- ceiling(sqrt(size))
          nsym <- round(nsym / favour, 0L)
          ncls <- round(ncls * favour, 0L)
          if (nsym > length(SYMBOLS) || ncls > maxclrs) {
            msg <- sprintf(
              "Column '%s': # symbols (%i) or # colours (%i) inacceptable",
              name, nsym, ncls)
            if (favour >= 1) {
              ncls <- maxclrs
              nsym <- ceiling(size / ncls)
            } else {
              nsym <- length(SYMBOLS)
              ncls <- ceiling(size / nsym)
            }
            message(msg, sprintf(", trying %i/%i instead.", nsym, ncls))
          }
          clrs <- select_colours(ncls, FALSE) # NA treated below
          symbols <- product(SYMBOLS, clrs)
          clrs <- symbols[, 2L]
          if (anyNA(levels.default(x))) # ensure last selected position is white
            clrs[[size]] <- WHITE
          symbols <- symbols[, 1L]
        }
        
        annotation <- c(annot1, list(
          BACKBONE_HEIGHT = 0, # controls the height of the midline
          BACKBONE_COLOR = WHITE, # controls the colour of the midline
          # we are hiding it by drawing it white
          BORDER_WIDTH = borwid,
          HEIGHT_FACTOR = 1,
          LEGEND_COLORS = clrs[seq_len(size)],
          LEGEND_LABELS = levels.default(x),
          LEGEND_SHAPES = symbols[seq_len(size)],
          LEGEND_TITLE = pretty_str(name),
          SHOW_DOMAIN_LABELS = 0,
          WIDTH = 25
        ))
        assert_no_forbidden_character("|", x)
        joint <- paste(symbols[x], 0L, 10L, clrs[x], as.character(x), sep = "|")
        print_itol(outdir, "domains", annotation, ids, 10L, joint)
        
      } else {
        
        annotation <- c(annot1, list(
          BORDER_WIDTH = borwid,
          LEGEND_COLORS = select_colours(size, anyNA(levels.default(x))),
          LEGEND_LABELS = levels.default(x),
          LEGEND_SHAPES = rep.int(1L, size),
          LEGEND_TITLE = pretty_str(name),
          STRIP_WIDTH = 25
        ))
        print_itol(outdir, "colorstrip", annotation,
                   ids, annotation$LEGEND_COLORS[x], x)
        
      }
    }
    
    
    # Integer vectors yield a bar chart.
    #
    emit_itol_integer <- function(x, ids, name, outdir, precision, ...) {
      if (!is.object(x)) # then we can assume it is really an integer vector
        precision <- NULL
      coordinated_na_removal(x, ids)
      annotation <- list(
        COLOR = BLACK,
        DATASET_LABEL = name,
        LEGEND_COLORS = BLACK,
        LEGEND_LABELS = paste0(legend_range(x, precision), collapse = ", "),
        LEGEND_SHAPES = 1L,
        LEGEND_TITLE = pretty_str(name),
        MARGIN = 5,
        WIDTH = 200
      )
      print_itol(outdir, "simplebar", annotation, ids, x)
    }
    
    
    # For treating vectors of mode 'double' like integers.
    #
    emit_itol_pseudointeger <- function(...) {
      emit_itol_integer(...)
    }
    
    
    # Should not normally occur in input.
    #
    emit_itol_list <- function(x, ids, name, outdir, ...) {
      message(sprintf("Skipping column '%s' of mode 'list' ...", name))
    }
    
    
    # For logical vectors. NAs do not get removed.
    #
    emit_itol_logical <- function(x, ids, name, outdir, bincolor, binsymbol,
                                  borwid, ...) {
      annotation <- list(
        BORDER_WIDTH = borwid,
        COLOR = "#4daf4a",
        DATASET_LABEL = name,
        FIELD_COLORS = bincolor,
        FIELD_LABELS = pretty_str(name),
        FIELD_SHAPES = binsymbol,
        LEGEND_COLORS = bincolor,
        LEGEND_LABELS = pretty_str(name),
        LEGEND_SHAPES = 1L,
        LEGEND_TITLE = pretty_str(name),
        MARGIN = 5,
        WIDTH = 20
      )
      print_itol(outdir, "binary", annotation, ids, convert_to(x, -1L))
    }
    
    
    # Vectors of mode 'double' (of class 'numeric' in R) yield a colour gradient.
    #
    emit_itol_numeric <- function(x, ids, name, outdir, endcolor,
                                  precision, borwid, ...) {
      coordinated_na_removal(x, ids)
      annotation <- list(
        BORDER_WIDTH = borwid,
        COLOR = "#fb9a99",
        COLOR_MAX = endcolor,
        COLOR_MIN = LIGHTGREY,
        DATASET_LABEL = name,
        LEGEND_COLORS = c(LIGHTGREY, endcolor),
        LEGEND_LABELS = legend_range(x, precision),
        LEGEND_SHAPES = c(1L, 1L),
        LEGEND_TITLE = pretty_str(name),
        MARGIN = 5,
        STRIP_WIDTH = 50
      )
      print_itol(outdir, "gradient", annotation, ids, x)
    }
    
    
    emit_itol_matrix <- function(x, ids, name, outdir, endcolor,
                                 precision, borwid, ...) {
      if (is.integer(x))
        precision <- 0L
      annotation <- list(
        BORDER_WIDTH = borwid,
        COLOR = "#fb9a99",
        COLOR_MAX = endcolor,
        COLOR_MIN = LIGHTGREY,
        COLOR_NAN = WHITE,
        DATASET_LABEL = name,
        FIELD_LABELS = colnames(x),
        LEGEND_COLORS = c(LIGHTGREY, endcolor),
        LEGEND_LABELS = legend_range(x, precision),
        LEGEND_SHAPES = c(1L, 1L),
        LEGEND_TITLE = pretty_str(name),
        MARGIN = 5,
        STRIP_WIDTH = 50
      )
      x <- apply(X = convert_to(x, "X"), MARGIN = 1L,
                 FUN = paste0, collapse = OUTPUT_SEPARATOR)
      print_itol(outdir, "heatmap", annotation, ids, x)
    }
    
    
    # Vectors of class 'factor' yield instructions for naming subtrees.
    #
    emit_branch_annotation_factor <- function(x, ids, name, outdir, ...) {
      coordinated_na_removal(x, ids)
      print_itol(outdir, "labels", name, ids, as.character(x))
    }
    
    
    # Input is supposed to represent support values to be shown as text.
    #
    emit_branch_annotation_integer <- function(x, ids, name, outdir, branchpos,
                                               cutoff, restriction, ...) {
      coordinated_na_removal(x, ids)
      annotation <- list(DATASET_LABEL = name, MARGIN = 5, COLOR = "#bebada")
      mask <- mask_if_requested(x, cutoff, tolower(restriction))
      if (any(mask)) {
        x[mask] <- NA_real_
        coordinated_na_removal(x, ids)
      }
      # additional columns: position, colour, style, size_factor, rotation
      print_itol(outdir, "text", annotation, ids, x,
                 branchpos, BLACK, "normal", 0.75, 0)
    }
    
    
    # Should not normally occur in input.
    #
    emit_branch_annotation_list <- function(x, ids, name, outdir, ...) {
      message(sprintf("Skipping column '%s' of mode 'list' ...", name))
    }
    
    
    # Vectors of class 'logical' select rows for the generation of instructions
    # for collapsing subtrees.
    #
    emit_branch_annotation_logical <- function(x, ids, name, outdir, ...) {
      x[is.na(x)] <- FALSE
      print_itol(outdir, "collapse", name, ids[x])
    }
    
    
    # Vectors of mode 'double' (of class 'numeric' in R) yield a colour gradient
    # within the branch symbols.
    #
    emit_branch_annotation_numeric <- function(x, ids, name, outdir, branchpos,
                                               cutoff, restriction, symbol, endcolor, maxsize, precision, ...) {
      coordinated_na_removal(x, ids)
      legendcolors <- c(LIGHTGREY, endcolor)
      if (grepl("[A-Z]", restriction, FALSE, TRUE)) {
        mask <- mask_if_requested(x, cutoff, tolower(restriction))
        if (any(mask)) {
          x[mask] <- NA_real_
          xclrs <- plotrix::color.scale(x = x[!is.na(x)], extremes = legendcolors)
          coordinated_na_removal(x, xclrs, ids)
        } else {
          xclrs <- plotrix::color.scale(x = x, extremes = legendcolors)
        }
      } else {
        xclrs <- plotrix::color.scale(x = x, extremes = legendcolors)
        mask <- mask_if_requested(x, cutoff, restriction)
        if (any(mask)) {
          x[mask] <- NA_real_
          coordinated_na_removal(x, xclrs, ids)
        }
      }
      annotation <- list(
        COLOR = endcolor,
        DATASET_LABEL = name,
        LEGEND_TITLE = pretty_str(name),
        LEGEND_SHAPES = c(symbol, symbol),
        LEGEND_COLORS = legendcolors,
        LEGEND_LABELS = legend_range(x, precision),
        MAXIMUM_SIZE = maxsize
      )
      print_itol(outdir, "branchsymbols", annotation,
                 # columns: ID, symbol, size, colour, fill, position
                 ids, symbol, maxsize, xclrs, 1L, branchpos)
    }
    
    
    # Helper function for fix_column_types. Only accepts dates in the canonical
    # format.
    #
    character2timediff <- function(x) {
      present <- !is.na(x) & nzchar(x, TRUE)
      result <- as.Date(x, "%Y-%m-%d")
      if (sum(is.na(result[present])) * 2L > length(result[present]))
        return(NULL)
      as.double(result - min(result, na.rm = TRUE) + 1L)
    }
    
    
    # Useful when R does not get the type right because of special notations;
    # also for user-defined type modifications.
    #
    fix_column_types <- function(x, convint, convdbl) {
      
      # convert binary integer vectors to logical vectors
      if (!is.element(convint, c("keep", "double")))
        for (i in which(vapply(x, is.integer, NA)))
          if (all(is.element(x[, i], c(0L, 1L, NA_integer_))))
            storage.mode(x[, i]) <- "logical"
          
          # convert factors to logical vectors if values look like boolean values
          if (convint != "keep")
            for (i in which(vapply(x, is.factor, NA))) {
              values <- tolower(levels.default(x[, i]))
              truevalue <- NA_character_
              for (j in seq_len(nrow(BINARY_VALUE_SPELLINGS)))
                if (all(is.element(values, BINARY_VALUE_SPELLINGS[j, ]))) {
                  truevalue <- BINARY_VALUE_SPELLINGS[j, 1L]
                  break
                }
              if (!is.na(truevalue))
                x[, i] <- tolower(x[, i]) == truevalue
            }
          
          if (convdbl)
            for (i in which(vapply(x, is.double, NA)))
              class(x[, i]) <- "pseudointeger"
            
            # convert integers and logical vectors to other data types if requested;
            # convert factors that look like dates to to double vectors if requested
            switch(
              EXPR = convint,
              keep = NULL,
              none = for (i in which(vapply(x, is.logical, NA)))
                if (anyNA(x[, i]))
                  x[, i] <- factor(x[, i]),
              factor = {
                for (i in which(vapply(x, is.integer, NA)))
                  x[, i] <- factor(x[, i])
                for (i in which(vapply(x, is.logical, NA)))
                  x[, i] <- factor(x[, i])
              },
              double = {
                for (i in which(vapply(x, is.integer, NA)))
                  storage.mode(x[, i]) <- "double"
                for (i in which(vapply(x, is.logical, NA)))
                  x[is.na(x[, i]), i] <- FALSE
                for (i in which(vapply(x, is.factor, NA))) {
                  timediff <- character2timediff(as.character(x[, i]))
                  if (length(timediff))
                    x[, i] <- timediff
                }
              },
              stop(sprintf("invalid integer/logical vector conversion indicator '%s'",
                           convint))
            )
            
            x
    }
    
    
    # Helper function for itol_files().
    #
    assort <- function(f, x) {
      idx <- split.default(seq_along(f), f)
      if (length(x)) {
        result <- vector(typeof(x), length(f))
        for (i in idx)
          result[i] <- rep_len(x, length(i))
      } else {
        result <- vector("double", length(f))
        for (i in idx)
          result[i] <- seq_along(i) / (length(i) + 1L)
      }
      result
    }
    
    
    # Helper function for itol_files().
    #
    get_col <- function(name, x, strict) {
      if (length(name) != 1L)
        stop("need a single column name for identifying special column")
      result <- match(name, names(x), 0L)
      if (!result)
        if (strict)
          stop(sprintf(
            "selected column '%s' does not exist -- must select one of %s",
            name, paste0(sprintf("'%s'", names(x)), collapse = ", ")))
      else
        warning(sprintf("cannot find column '%s', skipping data", name))
      result
    }
    
    
    # Helper function for branch_annotation_files().
    #
    force_logical <- function(x) {
      result <- switch(
        EXPR = class(x),
        logical = x,
        factor = nzchar(as.character(x), TRUE),
        character = nzchar(x, TRUE),
        list = lengths(x) > 0L,
        as.logical(x)
      )
      result[is.na(result)] <- FALSE
      result
    }
    
    # Called by itol_files() instead of jumping to the main branch.
    #
    branch_annotation_files <- function(x, icol, jcol, scol, idpat, precision,
                                        outdir, strict, maxsize, restrict) {
      
      if (nzchar(restrict, FALSE)) {
        cutoff <- as.double(basename(restrict))
        restriction <- dirname(restrict)
        if (restriction == ".") # when only a number was provided
          restriction <- "upto"
      } else {
        cutoff <- NA_real_
        restriction <- "upto"
      }
      
      idpos <- get_col(icol, x, TRUE)
      jpos <- get_col(jcol, x, TRUE)
      
      if (length(scol) && all(nzchar(scol, FALSE))) {
        spos <- get_col(scol, x, strict)
        if (spos)
          x <- x[force_logical(x[, spos]), , drop = FALSE]
      }
      
      assert_no_forbidden_character("|", x[, idpos], x[, jpos])
      icol <- ifelse(is.na(x[, jpos]), sprintf(idpat, x[, idpos]),
                     paste(sprintf(idpat, x[, idpos]), sprintf(idpat, x[, jpos]), sep = "|"))
      x <- x[, -c(idpos, jpos), drop = FALSE]
      
      klass <- vapply(x, class, "")
      
      # normal columns, dispatch done according to data type (class)
      mapply(FUN = function(fun, ...) fun(...), x = x, name = names(x),
             fun = lapply(sprintf("emit_branch_annotation_%s", klass), match.fun),
             branchpos = assort(klass, NULL), SIMPLIFY = FALSE,
             symbol = assort(klass, BRANCH_SYMBOLS), USE.NAMES = FALSE,
             endcolor = assort(klass, SPECIAL_COLORS),
             MoreArgs = list(ids = icol, precision = precision, outdir = outdir,
                             maxsize = maxsize, cutoff = cutoff, restriction = restriction))
      
      invisible(TRUE)
    }
    
    
    ## Main
    
    
    # The main function, taking care of all columns of data frame 'x'.
    #
    itol_files <- function(x, lcol, bcol, icol, jcol, scol, idpat, precision,
                           maxsize, favour, strict, convint, convdbl, outdir, borwid, restrict) {
      
      # identifier column (mandatory in strict mode), step 1
      idpos <- get_col(icol, x, strict)
      if (!idpos)
        return(invisible(FALSE))
      names(idpos) <- names(x)[[idpos]]
      
      if (anyNA(x[, idpos]))
        x <- x[!is.na(x[, idpos]), , drop = FALSE]
      
      if (!all(dim(x))) {
        if (strict)
          stop("encountered empty data frame")
        else
          warning("skipping empty data frame")
        return(invisible(FALSE))
      }
      
      x <- fix_column_types(x, convint, convdbl)
      
      # must be done before the first use of 'outdir'
      if (!dir.exists(outdir))
        dir.create(outdir)
      
      # generate branch symbols, skip normal run
      if (length(jcol) && all(nzchar(jcol, FALSE)))
        return(branch_annotation_files(x = x, icol = icol, jcol = jcol,
                                       scol = scol, idpat = idpat, precision = precision, outdir = outdir,
                                       strict = strict, maxsize = maxsize, restrict = restrict))
      
      # identifier column (mandatory in strict mode), step 2
      icol <- x[, idpos]
      if (is.factor(icol))
        icol <- as.character(icol)
      icol <- sprintf(idpat, icol)
      
      # label column (mandatory in strict mode)
      lpos <- get_col(lcol, x, strict)
      if (lpos)
        emit_itol_labeltexts(x = x[, lpos], ids = icol,
                             name = names(x)[[lpos]], outdir = outdir)
      
      # background colour column (optional in strict mode)
      if (length(bcol) && all(nzchar(bcol, FALSE))) {
        cpos <- get_col(bcol, x, strict)
        if (cpos)
          emit_itol_labelcolors(x = x[, cpos], ids = icol,
                                name = names(x)[[cpos]], outdir = outdir)
      } else {
        cpos <- 0L
      }
      
      # symbol-defining column (optional in strict mode)
      symbols <- NULL
      if (length(scol) && all(nzchar(scol, FALSE))) {
        spos <- get_col(scol, x, strict)
        if (spos) {
          symbols <- x[, spos]
          if (!is.factor(symbols) || anyNA(symbols) ||
              length(levels.default(symbols)) > length(SYMBOLS)) {
            warning("column '", scol, "' is either not a factor or ",
                    "has too many levels to be used for deriving symbols")
            symbols <- NULL
          } else {
            symbols <- SYMBOLS[symbols]
          }
        }
      }
      
      # normal columns, dispatch done according to data type (class)
      x <- x[, -c(idpos, lpos, cpos), drop = FALSE]
      klass <- vapply(x, class, "")
      
      # join all remaining columns in a matrix when applicable
      if (length(klass) > 1L && all(duplicated.default(klass)[-1L]))
        switch(
          EXPR = klass[[1L]],
          integer =,
          numeric = { # 'pseudointeger' should not be accepted here
            klass <- "matrix"
            x <- list(as.matrix(x))
            names(x) <- names(idpos)
          },
          NULL
        )
      
      clrs <- assort(klass, SPECIAL_COLORS)
      mapply(FUN = function(fun, ...) fun(...), x = x, name = names(x),
             fun = lapply(sprintf("emit_itol_%s", klass), match.fun), endcolor = clrs,
             bincolor = clrs, binsymbol = assort(klass, seq_along(BINARY_SYMBOLS)),
             MoreArgs = list(ids = icol, precision = precision, outdir = outdir,
                             symbols = symbols, maxclrs = maxsize, favour = favour,
                             borwid = borwid), SIMPLIFY = FALSE, USE.NAMES = FALSE)
      
      invisible(TRUE)
      
    }
    
    check_R_version()
    
    # assignment of input colour vectors is solely by vector length
    for (clrs in read_colour_vectors(colour.file, length(COLOURS)))
      COLOURS[[length(clrs)]] <- clrs
    COLOURS[] <- lapply(COLOURS, standardize_colour, opacity)
    check_colour_vectors(COLOURS, TRUE)
    
    # any length allowed, last one wins
    for (clrs in read_colour_vectors(gradient.file, Inf))
      SPECIAL_COLORS <- clrs
    SPECIAL_COLORS <- standardize_colour(SPECIAL_COLORS, opacity)
    check_colour_vectors(list(SPECIAL_COLORS), FALSE)
    
    LIGHTGREY <- standardize_colour(LIGHTGREY, opacity)
    
    na.strings <- unlist(strsplit(na.strings, separator, TRUE), FALSE, FALSE)
    if (!length(na.strings))
      na.strings <- ""
    
    for (infile in infiles)
      # note that read_file() is supposed to return a list of data frames
      lapply(X = read_file(infile, separator, na.strings, quote),
             FUN = itol_files, bcol = background, precision = precision, lcol = label,
             icol = identifier, scol = emblems, idpat = template, maxsize = max.size,
             favour = favour, strict = abort, jcol = identifier2, borwid = width,
             outdir = if (nzchar(directory, FALSE)) directory else dirname(infile),
             restrict = restrict, convint = conversion, convdbl = double.to.bars)
    
    invisible(NULL)
    
  }
  
    #Define values that would be changed while app is working. Datasets are also stores as reactive values
    values <- reactiveValues(test_data = NULL, primary_data=NULL, assignment_data=NULL,
                             age_cuttoff = NULL, assignment_data_backup = NULL, age_min=NULL, 
                             age_max=NULL)
    
    #This observes event of gisaid dataset upload. If uploaded cleaned the age and stires as reactive value
    observeEvent(input$gisaid, {
        #Read the dataset
        values$primary_data <- read.delim(input$gisaid$datapath, header = TRUE, sep = '\t')
        #Rename columns to those we are going to use
        colnames(values$primary_data) <- c( "Virus.name", "Accession.ID", "Collection.date" , "Location" , "Host", "Additional.location.information",
                                     "Gender" , "Patient.age", "Patient.status","Passage", "Specimen","Additional.host.information",    
                                     "Lineage", "Clade")
        #Patient.age cleaning
        values$primary_data$Patient.status <- as.factor(values$primary_data$Patient.status)
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, ">|<", "")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "20 - 30", "25")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "18-49", "33.5")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "4 months", "0.4")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "20-30", "25")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "33 years 5 months", "33.5")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "36, 11 months", "36.9")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "59 years 1 months", "59.1")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "2020-1975", "NA")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "31 years 6 months", "31.5")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "45-49", "47")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "55-59", "57")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "75-79", "77")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "2 months", "0.2")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "2020-1976", "NA")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "50-54", "52")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "51 years 3 months", "51.3")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "75, 0.4", "75.4")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "2020-1986", "NA")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "25-29", "27")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "30-34", "32")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "50-64", "57")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "61, 11 months", "61.9")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "65-69", "67")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "7 months", "0.6")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "years", "")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "s|'s", "")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "1 month", "0.1")
        values$primary_data$Patient.age <- str_replace_all(values$primary_data$Patient.age, "10-20", "15")
        values$primary_data$Patient.age <- as.numeric(values$primary_data$Patient.age) 
        #Drop the values for age less than 1 year
        values$primary_data <- values$primary_data %>%
          filter(Patient.age >= 1)
        #Write file
        write.csv(values$primary_data, file="primary_data.csv", row.names = F)
    })
    #Observe atatus_assignment.csv file upload
    observeEvent(input$status,{
      #Read the cleaned primary_data file
      values$primary_data <- read.csv("primary_data.csv")
      #Read the uploaded data
        values$assignment_data <- read.csv(input$status$datapath, header = F)
        #Writes csv file
        write.csv(values$assignment_data, file = "assignment_data.csv", row.names = F)
        values$assignment_data_backup <- data.frame(values$assignment_data)
        #Select first two rows
        values$assignment_data <- values$assignment_data %>%
            select(V1, V2)
        #Some cleaning of Patient.status column in primary_data
        values$primary_data$Patient.status <- as.character(values$primary_data$Patient.status)
        values$primary_data$Patient.status <- trimws(values$primary_data$Patient.status, which=c("both"))
        values$primary_data$Patient.age <- as.numeric(values$primary_data$Patient.age)
        #Map values from status_assignment to primary_data
        values$primary_data$Patient.status<- values$primary_data$Patient.status %>%
            mapvalues( from = c(values$assignment_data$V1), to = c(values$assignment_data$V2) )
        #Clean all the NAs
        values$primary_data$Patient.status[values$primary_data$Patient.status == "-"] <- ""
        values$primary_data$Patient.status <- as.factor(values$primary_data$Patient.status)
        values$primary_data$Patient.status <- droplevels(values$primary_data$Patient.status)
        values$primary_data$Patient.age <- as.numeric(values$primary_data$Patient.age)
        
        #copy dataset to another variable
        values$test_data <- values$primary_data
        #Create Coutry column
        values$test_data <- values$test_data %>%
            separate(Location,c("Continent", "Country", "Region", "City", "Small1", "Small2"), " / ") %>%
            select(-c(Continent, Region, City, Small1, Small2)) %>%
            drop_na()
        values$test_data$Country <- as.factor(values$test_data$Country)
        values$test_data$Patient.age <- as.numeric(values$test_data$Patient.age)
        values$test_data$Patient.status <- as.factor(values$test_data$Patient.status)
        values$test_data$Country <- as.factor(values$test_data$Country)
        values$test_data$Patient.status[values$test_data$Patient.status==""]<-NA
        #Drop NAs from the dataset. REady for calculations
        values$test_data <- drop_na(values$test_data)
    })
    # If Hospitalised status is changed or alive/live/.. box is checked reads the local datasets, 
    # change status_assignment.csv values  and remap the values to Patient.status category
    observeEvent({input$status_choice
      input$more_data},ignoreInit = TRUE, {
      #Read stored files
      values$primary_data <- read.csv("primary_data.csv")
        values$assignment_data=read.csv("assignment_data.csv")
        values$assignment_data <- values$assignment_data %>%
            select(V1, V2)
        #Logic how to rename categories in status_assignment dataset
        if (as.numeric(input$status_choice) == 1) {
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hopsitalized', 'Low risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalised', 'Low risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitaized', 'Low risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized', 'Low risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized', 'Low risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized or to be hospitalized', 'Low risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized, oxygenotherapy, diarrhea', 'Low risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized, released', 'Low risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized/Released', 'Low risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized; Stable', 'Low risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='In-hospital', 'Low risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Initially hospitalized, but now improved and discharged', 'Low risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Released', 'Low risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Released', 'Low risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Released, Live', 'Low risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Still hospitalized', 'Low risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='hospitalized', 'Low risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='hospitalized or to be hospitalized', 'Low risk', V2))
        }
        if (as.numeric(input$status_choice)  == 2) {
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hopsitalized', 'High risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitaized', 'High risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalised', 'High risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized', 'High risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized', 'High risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized or to be hospitalized', 'High risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized, oxygenotherapy, diarrhea', 'High risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized, released', 'High risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized/Released', 'High risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized; Stable', 'High risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='In-hospital', 'High risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Initially hospitalized, but now improved and discharged', 'High risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Released', 'High risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Released', 'High risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Released, Live', 'High risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Still hospitalized', 'High risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='hospitalized', 'High risk', V2))
            values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='hospitalized or to be hospitalized', 'High risk', V2))
        }
        if (as.numeric(input$status_choice)  == 3) {
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hopsitalized', '-', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitaized', '-', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalised', '-', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized', '-', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized', '-', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized or to be hospitalized', '-', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized, oxygenotherapy, diarrhea', '-', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized, released', '-', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized/Released', '-', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Hospitalized; Stable', '-', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='In-hospital', '-', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Initially hospitalized, but now improved and discharged', '-', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Released', '-', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Released', '-', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Released, Live', '-', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Still hospitalized','-', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='hospitalized', '-', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='hospitalized or to be hospitalized', '-', V2))
        }
        #Logic how to treat alive/live/symptomatic patients
        if (as.numeric(input$more_data) == 1) {
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Live', 'Low risk', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Alive', 'Low risk', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Symptomatic', 'Low risk', V2))
          
        } 
        if (as.numeric(input$more_data) == 2){
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Live', '-', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Alive', '-', V2))
          values$assignment_data$V2 <- with(values$assignment_data, ifelse(V1=='Symptomatic', '-', V2))
        } 
        # Final values mapping and dataset cleaning
        values$primary_data$Patient.status <- as.character(values$primary_data$Patient.status)
        values$primary_data$Patient.status <- trimws(values$primary_data$Patient.status, which=c("both"))
        values$primary_data$Patient.age <- as.numeric(values$primary_data$Patient.age)
        values$primary_data$Patient.status<- values$primary_data$Patient.status %>%
            mapvalues( from = c(values$assignment_data$V1), to = c(values$assignment_data$V2) )
        values$primary_data$Patient.status[values$primary_data$Patient.status == "-"] <- ""
        values$primary_data$Patient.status <- as.factor(values$primary_data$Patient.status)
        values$primary_data$Patient.status <- droplevels(values$primary_data$Patient.status)
        values$primary_data$Patient.age <- as.numeric(values$primary_data$Patient.age)
        #Copy dataset to test_data. Will work on this dataset
        values$test_data <- values$primary_data
        values$test_data <- values$test_data %>%
            separate(Location,c("Continent", "Country", "Region", "City", "Small1", "Small2"), " / ") %>%
            select(-c(Continent, Region, City, Small1, Small2)) %>%
            drop_na()
        values$test_data$Country <- as.factor(values$test_data$Country)
        values$test_data$Patient.age <- as.numeric(values$test_data$Patient.age)
        values$test_data$Patient.status <- as.factor(values$test_data$Patient.status)
        values$test_data$Country <- as.factor(values$test_data$Country)
        values$test_data$Patient.status[values$test_data$Patient.status==""]<-NA
        values$test_data <- drop_na(values$test_data)
        write.csv(values$test_data, file = "current_dataset.csv", row.names = F)

    })
    
    #If "Generate ITOL annotations" button is clicked and tree, gisaid data and status_assignment are uploaded do:
    observeEvent(input$itol,{
      #Check if all data is uploaded
      #Do nothing if not
      req(input$gisaid)
      req(input$status)
      req(input$tree)
      
      #Logic which labels to use for annotation files
      if (as.numeric(input$tree_labels) == 1){
        
        #Make local dataset
        itol_data <- values$test_data
        #Make age categories
        itol_data$Patient.age[itol_data$Patient.age <= 10] <- "0-10"
        itol_data$Patient.age[itol_data$Patient.age > 10 & itol_data$Patient.age <= 20] <- "11-20"
        itol_data$Patient.age[itol_data$Patient.age > 20 & itol_data$Patient.age <= 30] <- "21-30"
        itol_data$Patient.age[itol_data$Patient.age > 30 & itol_data$Patient.age <= 40] <- "31-40"
        itol_data$Patient.age[itol_data$Patient.age > 40 & itol_data$Patient.age <= 50] <- "41-50"
        itol_data$Patient.age[itol_data$Patient.age > 50 & itol_data$Patient.age <= 60] <- "51-60"
        itol_data$Patient.age[itol_data$Patient.age > 60 & itol_data$Patient.age <= 70] <- "61-70"
        itol_data$Patient.age[itol_data$Patient.age > 70 & itol_data$Patient.age <= 80] <- "71-80"
        itol_data$Patient.age[itol_data$Patient.age > 80 & itol_data$Patient.age <= 90] <- "81-90"
        itol_data$Patient.age[itol_data$Patient.age > 90] <- ">90"
        itol_data$Patient.age <- as.factor(itol_data$Patient.age)
        #Select data subset
        itol_data <- itol_data %>%
          select(Accession.ID, Patient.status, Clade, Country, Gender, Patient.age, Lineage) 
        itol_data$Country <- as.factor(itol_data$Country)
        itol_data$Clade <- as.factor(itol_data$Clade)
        itol_data$Lineage <- as.factor(itol_data$Lineage)
        #Write this data
        write.table(itol_data, file = "itol.files.tsv", sep = '\t', row.names = F )
        #Use written data to generate ITOL annotation files
        create_itol_files(infiles = c("itol.files.tsv"), identifier = "Accession.ID", label = "Accession.ID", na.strings = "X",background = input$itol_background )
      }
      if (as.numeric(input$tree_labels) == 0){
        #Make local dataset
        itol_data <- values$test_data
        #Make age categorie
        itol_data$Patient.age[itol_data$Patient.age <= 10] <- "0-10"
        itol_data$Patient.age[itol_data$Patient.age > 10 & itol_data$Patient.age <= 20] <- "11-20"
        itol_data$Patient.age[itol_data$Patient.age > 20 & itol_data$Patient.age <= 30] <- "21-30"
        itol_data$Patient.age[itol_data$Patient.age > 30 & itol_data$Patient.age <= 40] <- "31-40"
        itol_data$Patient.age[itol_data$Patient.age > 40 & itol_data$Patient.age <= 50] <- "41-50"
        itol_data$Patient.age[itol_data$Patient.age > 50 & itol_data$Patient.age <= 60] <- "51-60"
        itol_data$Patient.age[itol_data$Patient.age > 60 & itol_data$Patient.age <= 70] <- "61-70"
        itol_data$Patient.age[itol_data$Patient.age > 70 & itol_data$Patient.age <= 80] <- "71-80"
        itol_data$Patient.age[itol_data$Patient.age > 80 & itol_data$Patient.age <= 90] <- "81-90"
        itol_data$Patient.age[itol_data$Patient.age > 90] <- ">90"
        itol_data$Patient.age <- as.factor(itol_data$Patient.age)
        #Select data subset
        itol_data <- itol_data %>%
          select(Virus.name, Patient.status, Clade, Country, Gender, Lineage)
        itol_data$Country <- as.factor(itol_data$Country)
        itol_data$Clade <- as.factor(itol_data$Clade)
        itol_data$Lineage <- as.factor(itol_data$Lineage)
        #Write this data
        write.table(itol_data, file = "itol.files.tsv", sep = '\t', row.names = F )
        #Use written data to generate ITOL annotation files
        create_itol_files(infiles = c("itol.files.tsv"), identifier = "Virus.name", label = "Virus.name", na.strings = "X",background = input$itol_background )
        
        
      }
      
      #Find iTOL annotation files + uploaded nwk tree, zip them  
      flst <- c()
      files_in_dir <- list.files()
      for (file_names in files_in_dir) {
        if (grepl('iTOL', file_names, fixed = TRUE)) {
          flst <- c(flst, file_names)
        }
      } 
      
      zip(zipfile = "itol_annotations.zip", files = flst)
      
    })
    
    
    
    #Output plots
    output$oddsPlot <- renderPlotly({
      #Do the staff only if gisaid dataset and status_assignment file are uploaded
        req(input$gisaid)
        req(input$status)
        
        #OR calculations
        ages = 1
        odds_old <- c()
        odds_min <- c()
        odds_max <- c()
        #Loop to generate vector of odds ratio calculations of every
        # age theshold from 1 to max(Patient.age)-1
        
        for (ages in seq(1:(max(values$test_data$Patient.age)-1))) {
          age <- values$test_data[,c(8,9)] %>%
            drop_na() %>%
            mutate(age=ifelse(Patient.age > ages, 1, 0)) %>%
            select(-Patient.age, status=Patient.status) %>%
            mutate(col_name=ifelse(age == 1, "old", "young")) %>%
            pivot_wider( names_from = col_name, values_from = age) %>%
            transmute( status=status, ">cutoff"=lengths(old), "<cutoff"=lengths(young)) %>%
            column_to_rownames('status')
          age_new <- epi.2by2(t(as.matrix(age))[c(">cutoff", "<cutoff"),c("High risk","Low risk")])
          age_new <- as.data.frame(age_new$res$OR.strata.wald)
          odds_old <- c(odds_old, age_new$est)
          odds_min <- c(odds_min,age_new$lower)
          odds_max <- c(odds_max,age_new$upper)
          ages <- ages +1
        }
        
        odds_old <- as.data.frame(odds_old) 
        odds_old["Age_threshold"] <- seq(1:(max(values$test_data$Patient.age)-1))
        odds_old["Min"] <- odds_min
        odds_old["Max"] <- odds_max
        odds_plot <- odds_old %>%
          ggplot(aes(x=Age_threshold, y=as.numeric(odds_old),
                     lower = Min, upper=Max)) + 
          geom_point(size = 3, alpha=.7, color="#69b3a2") +
          ggtitle("Odds ratio of being at high risk given the age threshold") +
          ylab("Odds ratio value") +
          xlab("Patient age threshold") 
        
        ggplotly(odds_plot)
        
    })
    #Print individual OR for given age cutoff
    output$age_cut <- renderText({
      #Only if required files are uploaded do:
      req(input$gisaid)
      req(input$status)
      values$age_cuttoff <- input$age_cutoff
      #OR calculation but with no loop with age_cutoff values
      age <- values$test_data[,c(8,9)] %>%
        drop_na() %>%
        mutate(age=ifelse(Patient.age > values$age_cuttoff, 1, 0)) %>%
        select(-Patient.age, status=Patient.status) %>%
        mutate(col_name=ifelse(age == 1, "old", "young")) %>%
        pivot_wider( names_from = col_name, values_from = age) %>%
        transmute( status=status, ">cutoff"=lengths(old), "<cutoff"=lengths(young)) %>%
        column_to_rownames('status')
      age_new <- epi.2by2(t(as.matrix(age))[c(">cutoff", "<cutoff"),c("High risk","Low risk")])
      age_new <- as.data.frame(age_new$res$OR.strata.wald)
      #Store CI values in global reactive variables
      values$age_min <- age_new$lower
      values$age_max <- age_new$upper
      paste("The odds ratio with the given age cutoff are",age_new$est)
  
      
      
    })
    #Render Text 
    output$age_mes <- renderText({
      req(input$gisaid)
      req(input$status)
      
      "With the given confidence interval (95%):"
      
      
      
    })
    
    #Render text from global variable (min CI)
    output$age_min <- renderText({
      req(input$gisaid)
      req(input$status)
      
      paste("Lower: ", values$age_min)
      
      
    })
    
    #Render text from global variable (max CI)
    output$age_max <- renderText({
      req(input$gisaid)
      req(input$status)
      
      paste("Upper: ", values$age_max)
      
      
    })
    
    #Generate tree with help of ggtree package
    output$tree <- renderPlot({
      #Do only if all the data is uploaded
        inFile_tree <- input$tree
        if(is.null(inFile_tree))     return(NULL)
        req(input$gisaid)
        req(input$status)
        req(input$tree)
        
        
        
        #Read tree
        tree <- read.tree(inFile_tree$datapath)%>% as.phylo()
        #Write tree
        write.tree(tree, file = "tree.nwk") 
        #Convert to treedata object
        tree <- tree %>% as.treedata()
        #Tree calculations
        if (as.numeric(input$tree_labels) == 0) {
          tree_subset <- values$test_data %>%
            select(Virus.name, Patient.status)%>%
            dplyr::rename(label=Virus.name)
        } 
        if (as.numeric(input$tree_labels) == 1){
        tree_subset <- values$test_data %>%
            select(Accession.ID, Patient.status) %>%
            dplyr::rename(label=Accession.ID)
        }
        #remove spaces in manes from Patient.status
        tree_subset <- apply(tree_subset,2,function(x)gsub('\\s+', '_',x))
        
      
        #Join tree and Patient.status data
        tree <- treeio::full_join(tree, tree_subset, by="label")
        
          ggtree(tree, aes(color=Patient.status), layout = "circular") +theme(plot.margin = unit(c(-5,-4,-6,-7), "cm"))
        
    })
    
    #Write currently used dataset to download hangler
    output$current_dataset <- downloadHandler(filename = function (){
      paste("current_dataset.csv")
    },
    content = function(file){
      write.csv(values$test_data, file = file, row.names = F)
    })
    
    #Download iTOL files in zip archive
    output$itol_download <- downloadHandler(filename = function(){
     paste("itol_annotations.zip")     
    },  
                                            content =  function(file){
                                              flst <- c()
                                              files_in_dir <- list.files()
                                              for (file_names in files_in_dir) {
                                                if (grepl('iTOL', file_names, fixed = TRUE)) {
                                                  flst <- c(flst, file_names)
                                                }
                                              
                                                if (grepl('nwk', file_names, fixed = TRUE)) {
                                                  flst <- c(flst, file_names)
                                                  print(file_names)
                                                }
                                                }
                                              #create the zip file
                                              zip(file,  flst) },
    contentType = "application/zip" )
    
    #Render interactive tree with phylocanvas
    output$tree_interactive <- renderPhylocanvas({
      #Check if all required daat is uploaded
      inFile_tree <- input$tree
      if(is.null(inFile_tree))     return(NULL)
      req(input$gisaid)
      req(input$status)
      req(input$tree)
      
      #Read tree
      tree <- read.tree(inFile_tree$datapath) %>% as.phylo()
      
      if (as.numeric(input$tree_labels) == 0) {
        tree_subset <- values$test_data %>%
          select(Virus.name, Patient.status) %>%
          dplyr::rename(label=Virus.name)
      } 
      if (as.numeric(input$tree_labels) == 1){
        tree_subset <- values$test_data %>%
          select(Accession.ID, Patient.status) %>%
          dplyr::rename(label=Accession.ID)
      }
      high_risk_nodes <- tree_subset %>%
        filter(Patient.status == "High risk")%>%
        select(label) %>%
        as.vector()
      
      low_risk_nodes <- tree_subset %>%
        filter(Patient.status == "Low risk")%>%
        select(label) %>%
        as.vector()
      
      high_risk_nodes <- apply(high_risk_nodes,2,function(x)gsub('\\s+', '_',x))
      low_risk_nodes <- apply(low_risk_nodes,2,function(x)gsub('\\s+', '_',x))
      
      phycanv <- phylocanvas(tree, treetype = "circular", nodesize = 1, textsize = 1)
      
      #Render tree
      phycanv
    })
    
    #Outputs age density plot with vertical red line
    output$ageDensity <- renderPlot({
      #check is required datasets are uploaded
        req(input$gisaid)
        req(input$status)
        #Read age slider input
        values$age_cuttoff <- input$age_cutoff
        #Plot the data with give age cutoff
        values$test_data %>% 
            ggplot(aes(x=Patient.age, fill=Patient.status, y=..count..)) +
            geom_density(adjust=1.5, alpha=.4) +
            geom_vline(xintercept = values$age_cuttoff, color='red')+
            ggtitle("Age density") +
            ylab("Number of patients") +
            xlab("Patient age")
    })
    
    #Plots Countries count for given age cutoff
    output$Countries <- renderPlot({
      #check if required packages are uploaded
      req(input$gisaid)
      req(input$status)
      #Read data from age cutoff slider
      values$age_cuttoff <- input$age_cutoff
      #Make data transformations to plot. Binary age coding
      age_cuttoff_data <- values$test_data %>%
        mutate(age=ifelse(Patient.age > values$age_cuttoff, 1, 0)) %>%
        mutate(col_name=ifelse(age == 1, "old", "young"))
      #Plot carplot for countries that are younger
      o_ctr <- age_cuttoff_data %>%
        filter(age ==1)%>%
        pivot_wider( names_from = col_name, values_from = age) %>%
        count('Country') %>%
        mutate(Country = fct_reorder(Country, freq)) %>%
        ggplot(aes(x=Country, y=freq )) +
        geom_bar(stat = 'identity',  fill="#f68060", alpha=.6, width=.4) +
        coord_flip() +
        ggtitle(paste0("Count of cases per country at age > ", values$age_cuttoff, " y")) +
        xlab("") +
        ylab('Count') +
        geom_text( aes(y=freq+4, label = freq))
      #Plot barplot for countries that are older
      y_ctr <- age_cuttoff_data %>%
        filter(age == 0) %>%
        pivot_wider( names_from = col_name, values_from = age) %>%
        count('Country') %>%
        mutate(Country = fct_reorder(Country, freq)) %>%
        ggplot(aes(x=Country, y=freq )) +
        geom_bar(stat = 'identity',  fill="#f68060", alpha=.6, width=.4) +
        coord_flip() +
        ggtitle(paste0("Count of cases per country at age < ", values$age_cuttoff, " y")) +
        xlab("") +
        ylab('Count') +
        geom_text( aes(y=freq+4, label = freq))
      
      #Make a list of plots
        ptlist <- list(o_ctr, y_ctr) 
        #Plot the plots one near another
        grid.arrange(grobs=ptlist,ncol=length(ptlist))

    })
    
    
    
    
    #Render boxplot chart for age
    output$boxplot_age <- renderPlot({
      #Check if all required data is uploded
        req(input$gisaid)
        req(input$status)
        #Plot age distribution boxplot to a variable
       pt1 <- values$test_data %>%
            ggplot( aes(x=Patient.status, y=Patient.age, fill=Patient.status)) +
            geom_boxplot() +
            scale_fill_viridis(discrete = TRUE, alpha=0.6) +
            geom_jitter(color="black", size=0.4, alpha=0.9) +
            theme_ipsum() +
            theme(
                legend.position="none",
                plot.title = element_text(size=11)
            ) +
            ggtitle("Distribution of age over High/Low risk categories") +
            xlab("")
       
       
       #Risk calculations. To accompany boxplot
       low_risk <- values$test_data %>%
         filter(Patient.status=="Low risk") %>%
         select(Patient.age) %>%
         summary()
       high_risk <- values$test_data %>%
         filter(Patient.status=="High risk") %>%
         select(Patient.age) %>%
         summary()
       high_risk <- as.data.frame(high_risk)  %>%
         separate(Freq, c("Name", "High risk"), ":") %>%
         select(Name,"High risk")
       low_risk <- as.data.frame(low_risk)  %>%
         separate(Freq, c("Name", "Low risk"), ":") %>%
         select("Low risk")
       #Actual plot
       risk_stats <- cbind( high_risk, low_risk)
      #Create ggplot object of a table
       pt2 <- as.ggplot(tableGrob(risk_stats))
       ptlist <- list(pt1, pt2)  
       #Plot boxplot and a table side-by-side
       grid.arrange(grobs=ptlist,ncol=length(ptlist))
    })
    
    #Ouput clade distribution plots
    output$cladeAge <- renderPlot({
      req(input$gisaid)
      req(input$status)
      
      #First plot with distribution of age over Clades
      pt1 <- values$test_data %>%
        ggplot( aes(x=Clade, y=Patient.age, fill=Clade)) +
        geom_boxplot() +
        scale_fill_viridis(discrete = TRUE, alpha=0.6) +
        geom_jitter(color="black", size=0.4, alpha=0.9) +
        theme_ipsum() +
        theme(
          legend.position="none",
          plot.title = element_text(size=11)
        ) +
        ggtitle("Distribution of age over Clades") +
        xlab("") +
        #Significance comparisons
        geom_signif(comparisons = list(c("G", "GH"), c("GH", "GR"),c("GR", "L"), c("L", "O"), c("O", "S"), c("S", "V")), map_signif_level = T)
      #Second plot with distribution of age over Clades with Patient.status assignment
      pt2 <- values$test_data %>%
        ggplot( aes(x=Clade, y=Patient.age, fill=Patient.status)) +
        geom_boxplot() +
        scale_fill_viridis(discrete = TRUE, alpha=0.6) +
        geom_jitter(color="black", size=0.4, alpha=0.9) +
        theme_ipsum() +
        theme(
          plot.title = element_text(size=11)
        ) +
        ggtitle("Distribution of age over Clades with Patient.status assignment") +
        xlab("") +
        #Statistical significance comparison between groups within 1 clade
        stat_compare_means (label = "p.signif")
      
      
      ptlist <- list(pt1, pt2)  
      #Plot side-by-side
      grid.arrange(grobs=ptlist,ncol=length(ptlist))
      
    })
    
    #Create Country per clade counts interactive plot
    output$country_count <- renderPlotly({
      req(input$gisaid)
      req(input$status)
      #Create empty dataframe to fill it in with country counts
      country_count <- data.frame()
      #Iterate over countries and count Clade frequency for each country
      for (country in as.vector(levels(values$test_data$Country))) {
        tmp <- values$test_data %>%
          select(Country, Clade) %>%
          filter(Country == country) %>%
          count()
        
        country_count <- rbind(country_count, tmp)
      }
      
      #Plot the dataframe into variable
      plot_1 <- country_count %>%
        ggplot(aes(x=Country, y=freq,  color=Clade, group=Clade)) +
        geom_point(alpha=0.7) +
        theme(axis.text.x = element_text(face = "bold", color = "#993333", 
                                         angle = 45,  hjust = 1)) +
        xlab("") + 
        ylab("Patient count")
      #Plot interactive plot
      ggplotly(plot_1)
      
    })
    
    #Create country odds ratio interactive plot
    output$country_odds <- renderPlotly({
      req(input$gisaid)
      req(input$status)
      #Country odds ratio calculation
      #Construct the overal dataframe for High/Low risks calculations
      status_freq <- values$test_data %>%
        select(Patient.status) %>%
        count() %>%
        pivot_wider(names_from = Patient.status, values_from=freq)
      #tmp dataframe to store clade count information (general)
      clade_freq_tmp <- values$test_data %>%
        select(Clade) %>%
        count()
      #Get the median value for patient count in single clade
      clade_freq_median <- median(clade_freq_tmp$freq)
      #Get the ratio of high/low categories. Use this ratio and median value for the clade to calculate reference clade
      status_freq <-  status_freq %>% 
        mutate(ratio=`High risk`/`Low risk`, part=1/(length(levels(as.factor(values$test_data$Clade)))+1)) %>%
        transmute(`High risk`=(clade_freq_median/(ratio+1)) * ratio, `Low risk`=median(clade_freq_tmp$freq)-`High risk`)
      #add clade name to the reference dummy clade
      status_freq$Clade <- "A_ref"
      #Calculate clade frequencies
      clade_freq <- values$test_data %>%
        select(Clade, Patient.status) %>%
        count() %>%
        pivot_wider(names_from = Patient.status, values_from=freq)
      #Join reference clade with the clade information
      whole_freq <- rbind(status_freq, clade_freq)
      #Store clade names to the value
      clade_names <- as.vector(whole_freq$Clade)
      #Select only risj categories for oddsratio calculations
      whole_freq <- whole_freq %>%
        select(`High risk`, `Low risk`)
      #oddsratio calculations
      whole_seq_odds <- suppressWarnings(oddsratio(as.matrix(whole_freq)))
      #Construct the new dataframe from Clade names and oddsratio calculation (columnbind the measure and p_values)
      country_odds <- cbind(Clade=clade_names,as.data.frame(whole_seq_odds$measure), as.data.frame(whole_seq_odds$p.value))
      country_odds$Clade <- as.factor(country_odds$Clade)
      #Plot the data into variable
      country_odds_plot <- country_odds %>%
        ggplot(aes(x=Clade, y=estimate, lower=lower, upper=upper, chi=chi.square, fisher=fisher.exact, mid=midp.exact)) + 
        geom_point(size = 3, alpha=.7, color="violet") +
        ggtitle("Odds ration of being at high risk given the Clade") +
        ylab("Odds ratio") +
        xlab("Clade") +
        geom_errorbar(aes(ymin=lower, ymax = upper))
      #Plot the interactive plot
      ggplotly(country_odds_plot)
      
    })
    

}  

# Run the application 
shinyApp(ui = ui, server = server)
