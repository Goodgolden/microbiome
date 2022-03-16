# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/


## packages ---------------------------------------------
# Sun Feb 14 22:42:58 2021 ------------------------------
library("shiny")
library("tidyverse")
library("here")


## options ---------------------------------------------
options(shiny.reactlog = TRUE)


## datasets ---------------------------------------------
# Sun Feb 14 22:43:17 2021 ------------------------------
merge_eligible <- 
  here::here(
    "seven_twenty_four",
    "new_data", 
    "724_merge_eligible.csv") %>%
  vroom::vroom() %>%
  select(Lib = 1, everything())

merge <- 
  here::here(
    "seven_twenty_four",
    "new_data",
    "724_merge_eligible.csv") %>%
  read.csv(row.names = 1,
           check.names = F,
           header = TRUE) 

alpha <- 
  here::here(
    "seven_twenty_four",
    "new_data", 
    "724_alpha_eligible.csv") %>%
  read.csv(row.names = 1,
           check.names = F,
           header = TRUE) %>%
  rownames_to_column("Lib")

meta <- 
  here::here(
    "seven_twenty_four",
    "data", 
    "724_meta.csv") %>%
  read_csv() %>%
  select(SYSID, PatNo, everything()) %>%
  as.data.frame()

merge_alpha <- 
  inner_join(meta, alpha, 
             by = "Lib") %>%
  filter(!is.na(SYSID))


## functions ---------------------------------------------
# Sun Feb 14 22:44:06 2021 ------------------------------

## a function from Kayla
## only keep the last level
## of the OTU_Names
name_split <-
  function(names) {
    names <- names
    save <- strsplit(names, "/")
    h <- 0
    for (i in 1:length(names)) {
      h[i] <- length(save[[i]])
    }
    i <- 0
    name.list <- NULL
    for (i in 1:length(save)) {
      name.list[i] <- save[[i]][h[i]]
    }
    return(name.list)
  }


## a function to subset
## easier to use in map()
select_eligible <-
  function(ID, dataset) {
    dataset %>%
      filter(SYSID == ID) %>%
      select(SYSID, VISITNUM, contains("/")) %>%
      rownames_to_column("Lib") %>%
      unite("SYS", c("SYSID", "VISITNUM", "Lib")) %>%
      janitor::adorn_totals("row") %>%
      remove_rownames() %>%
      column_to_rownames("SYS") %>%
      t() %>%
      as.data.frame() %>%
      arrange(desc(Total)) %>%
      head(12) %>%
      select(-Total) %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("SYS") %>%
      reshape2::melt(
        id.vars = "SYS",
        variable.name = "Taxa"
      ) %>%
      mutate(Taxa = name_split(as.character(Taxa))) %>%
      separate("SYS", into = c("SYSID", "VISITNUM", "Lib"))
  }

## easier to use in map()
select_ID <-
  function(ID, dataset) {
    dataset %>%
      filter(SYSID == ID)
  }


## shinyapp ---------------------------------------------
# Sun Feb 14 22:45:43 2021 ------------------------------

# ui for application that draws a histogram
ui <- fluidPage(
  ## title ---------------------------------------------
  titlePanel("Dashboard for 724 microbiome"),
    fluidRow(
      column(2,
             selectInput("individual",
                         "Patient",
                         choices = merge_eligible$SYSID,
                         width = "100%")),
      column(1,
             actionButton("run",
                          "Run!",
                          icon = icon("refresh"),
                          class = "btn-lg btn-success")),
      column(3,
             tableOutput("information")),
      column(3,
             tableOutput("biomarker")),
      column(1,
             downloadButton("download_pdf"))),

  ## table panels -----------------------------------------
   fluidRow(
     column(3,
            selectInput("items",
                        "Variables",
                        choices = colnames(merge_eligible[3:128]),
                        multiple = TRUE,
                        width = "100%")),
      column(9,
             tableOutput("immune"))),

    fluidRow(
      column(4,
             plotOutput("microbiome")),
      column(4,
             plotOutput("diversity")),
      column(4,
             plotOutput("overall"))))



server <- function(input, output) {
  selected <- reactive(input$individual)

  observeEvent(input$run, {  })

  output$information <- renderTable({
    req(input$individual, input$run)

    merge_eligible %>%
      filter(SYSID == input$individual) %>%
      select(SYSID, VISITNUM,
             fev = FEV_BEST,
             fvc = FVC_Best,
             frc = FRC_Best)})

  output$biomarker <- renderTable({
    req(input$individual, input$run)
    merge_eligible %>%
      filter(SYSID == input$individual) %>%
      select(il1= IL_1_FINALCONC,
             il17 = IL_17FINALCONC,
             il6 = IL_6FINALCONC,
             mmp2 = MMP2FINALCONC,
             mmp9 = MMP9FINALCONC)})

  output$immune <- renderTable({
    req(input$individual, input$run)
    merge_eligible %>%
      filter(SYSID == input$individual) %>%
      select(neutrophils = Neutrophils,
             macrophage = Monocyte__Macrophage,
             input$items)})

  summary1 <- reactive({
    req(input$individual, input$run)
    select_eligible(selected(), merge_eligible)
  })

  output$microbiome <- renderPlot({
    req(input$individual, input$run)

    summary1() %>%
      ggplot(aes(x = VISITNUM, y = value, fill = Taxa)) +
      ## should have saved this function
      ## probably add a ggsave() to pdf directly
      geom_bar(stat = "identity") +
      theme_classic() +
      ggthemes::scale_fill_tableau("Classic Cyclic") +
      ## so far the best color compositions for bar plot
      ## "Jewel Bright" only contains seven color
      ggthemes::scale_colour_tableau("Classic Cyclic") +
      ## the tableau is in the ggthemes
      labs(x = "Time of Visits") +
      labs(y = "Relative Abundance") +
      labs(title = unique(selected())) +
      ylim(0, 1)},
    res = 96)

    summary2 <- reactive({
      req(input$individual, input$run)
      select_ID(selected(), merge_alpha)
    })

    output$diversity <- renderPlot({
      req(input$individual, input$run)
      summary2() %>%
        ggplot(aes(x = VISITNUM, group = SYSID)) +
        geom_smooth(aes(y = `ShannonH Mean`, color = "ShannonH Mean"),  size = 1) +
        geom_smooth(aes(y = `ShannonE Mean`, color = "ShannonE Mean"), size = 1) +
        geom_smooth(aes(y = `Sobs Mean` / 20, color = "Sobs Mean"), size = 1) +
        geom_line(aes(y = `ShannonH Mean`, color = "ShannonH Mean"),  size = 1) +
        geom_line(aes(y = `ShannonE Mean`, color = "ShannonE Mean"), size = 1) +
        geom_line(aes(y = `Sobs Mean` / 20, color = "Sobs Mean"), size = 1) +
        geom_point(aes(y = `ShannonH Mean`, color = "ShannonH Mean"),  size = 2) +
        geom_point(aes(y = `ShannonE Mean`, color = "ShannonE Mean"), size = 2) +
        geom_point(aes(y = `Sobs Mean` / 20, color = "Sobs Mean"), size = 2) +
        theme_classic() +
        ylim(0, 5) +
        scale_y_continuous(name = "Shannon Mean", limits = c(0, 5),
                           sec.axis = sec_axis(trans = ~. * 20, name = "Sob Mean")) +
        labs(y = "Biodiversity") +
        labs(x = "Visit Times") +
        labs(title = selected()) +
        scale_color_discrete(name = "") +
        guides(color = guide_legend(override.aes=list(fill=NA))) +
        theme(legend.key = element_blank(),
              legend.background = element_blank(),
              legend.position = "bottom")},
    res = 96)



    output$overall <- renderPlot({
      req(input$individual, input$run)

      merge_alpha %>%
        mutate(VISITNUM = as.factor(VISITNUM)) %>%
        ggplot(aes(x = VISITNUM,
                   y = `ShannonE Mean`),
               group = VISITNUM) +
        geom_point(aes(color = VISITNUM)) +
        geom_boxplot(aes(fill = VISITNUM)) +
        theme_classic() +
        ggthemes::scale_fill_tableau("Jewel Bright") +
        ## so far the best color composition
        ggthemes::scale_colour_tableau("Jewel Bright") },
      res = 96)


    output$download_pdf <- downloadHandler(
      filename = "rendered_report.pdf",
      content = function(file) {
        res <- rmarkdown::render(
          "05_shiny.R",
          params = list(c(output$microbiome,
                          output$diversity,
                          output$overall)))
        file.rename(res, file)
      })

}

## run the app
# shinyApp(ui = ui, server = server)
