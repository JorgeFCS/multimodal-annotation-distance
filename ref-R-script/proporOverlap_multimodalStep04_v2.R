#### Select directories ####
#Separate gestures and transcription files in two different directories
choose.dir(caption = "Select gesture folder")
choose.dir(caption = "Select transcription folder")

#### Init ####
# Load packages
#lembrar de colocar install packages com if condition
library(dplyr)
library(ggplot2)
library(readr)
library(DescTools)
library(tidyverse)
library(rebus)

#Import data one file
transc_files <- list.files(path_transc, pattern = "*.txt", full.names = T)
utterances_raw <- sapply(transc_files, read_tsv, simplify=FALSE) %>% 
  bind_rows(.id = "id")

gestures_files <- list.files(path_gesto, pattern = "*.txt", full.names = T)
gestures_raw <- sapply(gestures_files, read_tsv, simplify=FALSE) %>% 
  bind_rows(.id = "id")


#### Data Wrangling ####
#Some data wrangling
#Change id to file name (You should change start parameter if filename is incorrect)
#Make sure that after wrangling the column file is correct. Otherwise, all other calculations may be incorrect
utterances_id<-utterances_raw$id
utterances_id<-sapply(utterances_id, function(x) {
  str_sub(x, start = nchar(path_transc)+2, end = -12)}, simplify = TRUE, USE.NAMES = FALSE)
utterances_raw$id<-utterances_id

gestures_id<-gestures_raw$id
gestures_id<-sapply(gestures_id, function(x) {
  str_sub(x, start = nchar(path_gesto)+2, end = -13)}, simplify = TRUE, USE.NAMES = FALSE)
gestures_raw$id<-gestures_id

#Cut off unused row
if("X14" %in% colnames(gestures_raw)){
  gestures_raw$X14 <- NULL
}

if("X5" %in% colnames(utterances_raw)){
  utterances_raw$X5 <- NULL
}

info <- utterances_raw$InfoStructure

utterances <- utterances_raw %>%
  select("id", "Begin Time - msec", "End Time - msec", "Transcription", "NTB") %>%
  group_by(Transcription) %>%
  mutate(IU_id = seq_along(Transcription)) %>%
  ungroup() %>%
  group_by(id, Transcription) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(utt_id = seq_len(n())) %>% 
  unnest(cols = data)

colnames(utterances) <- c("file", "utterance", "begin", "end", "iu", "iu_id", "utt_id")

gestures<-gestures_raw %>%
  mutate_all(~replace(., is.na(.), "?"))

colnames(gestures) <- c("file", "begin", "end", "ge_units", "ge_phrase", "ge_phase", "right_ge_position", "right_ge_handshape",
                        "right_ge_orientation", "right_ge_movement", "left_ge_position", "left_ge_handshape", "left_ge_orientation",
                        "left_ge_movement")
#, "right_handshape", "right_position", "right_orientation",      "right_movement")

gestures$begin<-as.numeric(gestures$begin)
gestures$end<-as.numeric(gestures$end)
gestures$ge_units<-as.numeric(gestures$ge_units)
gestures$ge_phrase<-as.numeric(gestures$ge_phrase)

#Check structure
str(utterances)
str(gestures)


#### Function Check Overlapping ####
#Function to check whether intervals overlap
check_overlapping <- function(begin1, end1, begin2, end2){
  x<-c(begin1, end1)
  y<-c(begin2, end2)
  overlapping_boolean <- as.logical(x %overlaps% y)
  overlapping_time <- Overlap(x, y)
  interval_time <- Interval(x, y)
  results<-c(overlapping_boolean, overlapping_time, interval_time)
  return(results)
}


#### Function Get Strokes ####
#Iterate over gestures to check least distant overlappings
#How to deal with "holds"? Should they be dealt with as part of the stroke?
get_strokes <- function(gestures) {
  gestures_wrangled <- data.frame(file=as.character(),
                                  ge_phrase=as.numeric(),
                                  stroke_begin=as.numeric(), 
                                  stroke_end=as.numeric(), 
                                  ge_phrase_begin=as.numeric(),
                                  ge_phrase_end=as.numeric(),
                                  stringsAsFactors=FALSE)
  i<-1
  for (i in 1:nrow(gestures)) {
    if (gestures$ge_phase[i] == "stroke") {
      row_index<-nrow(gestures_wrangled)+1
      gestures_wrangled[row_index,] <- NA
      gestures_wrangled$file[row_index]<-gestures$file[i]
      gestures_wrangled$ge_phrase[row_index]<-gestures$ge_phrase[i]
      gestures_wrangled$stroke_begin[row_index]<-gestures$begin[i]
      gestures_wrangled$stroke_end[row_index]<-gestures$end[i]
      stroke_ge_phrase<-gestures$ge_phrase[i]
      ge_phrase_backwards<-gestures$ge_phrase[i]
      pos_backwards <- i
      ge_phrase_forwards<-gestures$ge_phrase[i]
      pos_forwards <- i
      while (ge_phrase_backwards == stroke_ge_phrase) {
        gestures_wrangled$ge_phrase_begin[row_index]<-gestures$begin[pos_backwards]
        if (pos_backwards > 1) {
          pos_backwards <- pos_backwards-1
          ge_phrase_backwards<-gestures$ge_phrase[pos_backwards]
        } else {
          ge_phrase_backwards<-FALSE
        }
      }
      while (ge_phrase_forwards == stroke_ge_phrase) {
        gestures_wrangled$ge_phrase_end[row_index]<-gestures$end[pos_forwards]
        if (pos_forwards < nrow(gestures)) {
          pos_forwards<-pos_forwards+1
          ge_phrase_forwards<-gestures$ge_phrase[pos_forwards]
        } else {
          ge_phrase_forwards<-FALSE
        }
      }
    }
  }
  rm(row_index, stroke_ge_phrase, ge_phrase_backwards, pos_backwards, ge_phrase_forwards, pos_forwards, j, l)
  gestures_wrangled$file <- as.factor(gestures_wrangled$file)
  gestures_wrangled$ge_phrase <- as.factor(gestures_wrangled$ge_phrase)
  gestures_wrangled$stroke_begin <- as.numeric(gestures_wrangled$stroke_begin)
  gestures_wrangled$stroke_end <- as.numeric(gestures_wrangled$stroke_end)
  gestures_wrangled$ge_phrase_begin <- as.numeric(gestures_wrangled$ge_phrase_begin)
  gestures_wrangled$ge_phrase_end <- as.numeric(gestures_wrangled$ge_phrase_end)
  return(gestures_wrangled)
}


#### Define Function Get Closest ####
# Adcionar coluna com número do arquivo
#Create result dataframe called closest
get_closest <- function(gestures_wrangled_fun, utterances_fun) {
  closest <- gestures_wrangled_fun %>%
    add_column(overlap = F) %>%
    add_column(iu_begin = 0) %>%
    add_column(iu_end = 0) %>%
    add_column(overlapping_time = 0) %>%
    add_column(interval_time = 0) %>%
    add_column(utt_id = NA) %>%
    add_column(iu_id = NA) %>%
    add_column(iu_text = NA) %>%
    add_column(boundary_type = NA)
  
  # Find the IU that has the largest overlapping with each stroke from gestures_wrangled
  k <-1
  l <- 1
  for (k in 1:nrow(gestures_wrangled_fun)) {
    for (l in 1:nrow(utterances_fun)) {
      overlap <- check_overlapping(gestures_wrangled_fun$stroke_begin[k], gestures_wrangled_fun$stroke_end[k],
                                   utterances_fun$begin[l], utterances_fun$end[l])
      if (overlap[2] >= closest$overlapping_time[k]) {
        closest$overlap[k]<-overlap[1]
        closest$iu_begin[k]<-utterances_fun$begin[l]
        closest$iu_end[k]<-utterances_fun$end[l]
        closest$overlapping_time[k]<-overlap[2]
        closest$interval_time[k]<-overlap[3]
        closest$utt_id[k]<-utterances_fun$utt_id[l]
        closest$iu_id[k]<-utterances_fun$iu_id[l]
        closest$iu_text[k]<-utterances_fun$iu[l]
        # Boundary type only is correct when boundaries are transcribed
        closest$boundary_type[k]<-ifelse(test = str_detect(utterances_fun$iu[l], "//" %|% "\\+"),
                                         yes = "TB_or_interrupted",
                                         no = "NTB")
        
      }
    }
  }
  #Calculate boundaries differences
  closest_diff <- closest %>%
    mutate(diff_begin = ge_phrase_begin-iu_begin,
           diff_end = ge_phrase_end-iu_end,
           ge_phrase_dur = ge_phrase_end - ge_phrase_begin,
           overlapping_ratio = ifelse(test=overlapping_time>0, round(overlapping_time / ge_phrase_dur, 2), 0)) #*100 is to show the percent of overlap

  return(closest_diff)
}

#### Nearest IU by file ####
# This part of the code may take some time to run

#Split tibbles by file group
utterances_listed <- utterances %>%
  group_by(file) %>%
  group_split()

gestures_listed <- gestures %>%
  group_by(file) %>%
  group_split()

#Run the analyses on each group calling functions get_stroke and get_closest
#Results are saved in a list of the same length of file groups
if (length(utterances_listed) == length(gestures_listed)) {
  result_list <- list()
  for (i in 1:length(gestures_listed)) {
    gestures_unlisted <- gestures_listed[[i]]
    utterances_unlisted <- utterances_listed[[i]]
    strokes <- get_strokes(gestures_unlisted)
    closest <- get_closest(strokes, utterances_unlisted)
    result_list[[i]] <- closest
    rm(gestures_unlisted, utterances_unlisted, strokes, closest)
  }
} else {
  print("Gesture and Utterance tables don't have the same number of files!")
}

#Bind result list and clean workspace
results <- result_list %>%
  bind_rows()

results <- results %>%
  filter(overlap == 1)

#### Graphs ####
#List to store ggplot objects
plots<-list()

#Centered difference value histogram
plots[[1]]<-results %>%
  gather("differences", "time", 16:17) %>%
  filter(time <20000)%>%
  ggplot(aes(scale(time), fill=as.factor(differences))) +
  geom_histogram(alpha=.6, position = "identity", bins=20) +
  labs(x="Centered values", y="Frequency", title="Differences of Boundary Alignment between GE-Phrase and IU",
       fill="Boundary")+
  coord_flip()

#Centered difference value histogram by file
plots[[2]]<-results %>%
  gather("differences", "time", 16:17) %>%
  filter(time <20000)%>%
  ggplot(aes(scale(time), fill=as.factor(differences))) +
  geom_histogram(alpha=.6, position = "identity",  bins=20) + facet_grid(~file) +
  labs(x="Centered values", y="Frequency", title="Differences of Boundary Alignment between GE-Phrase and IU by file",
       fill="Boundary") +
  coord_flip()

#Centered difference value histogram
plots[[3]]<-results %>%
  gather("differences", "time", 16:17) %>%
  filter(time <20000)%>%
  ggplot(aes(time, fill=as.factor(differences))) +
  geom_histogram(alpha=.6, position = "identity", bins=20) +
  labs(x="Differences (ms)", y="Frequency", title= "Differences of Boundary Alignment between GE-Phrase and IU", fill="Boundary") +
  geom_vline(xintercept = 500)+
  geom_vline(xintercept = -500)
  coord_flip()

#Centered difference value histogram by file
plots[[4]]<-results %>%
  gather("differences", "time", 16:17) %>%
  filter(time <20000)%>%
  ggplot(aes(time, fill=as.factor(differences))) +
  geom_histogram(alpha=.6, position = "identity", bins=25) + facet_grid(~file) +
  labs(x="Differences (ms)", y="Frequency", title= "Differences of Boundary Alignment between GE-Phrase and IU by file", fill="Boundary") +
  coord_flip()

#Overlapping ratio histogram
plots[[5]]<-results %>%
  ggplot(aes(overlapping_ratio)) +
  geom_histogram(alpha=.6, position = "identity", bins=20) + 
  labs(x="Overlapping_ratio", y="Frequency", title = "Overlapping of GE-Phrase with respect to closest IU")

#Barplot of boundaries closer than 200 ms
plots[[6]]<-results %>%
  gather("differences", "time", 16:17) %>%
  mutate(is_close = ifelse(test = abs(time) <= 200, TRUE, FALSE)) %>%
  ggplot(aes(is_close, fill=as.factor(is_close))) + geom_bar() + facet_grid(~differences) +
  labs(x="Are boundaries closer than 200 ms?", y="Count", title = "Boundaries of GE-Phrase and IU that are less than 200 ms one from another",
       fill="< 200 ms")

#Print files to RStudio
for (k in 1:length(plots)) {
  print(plots[[k]])
}

#### Save output ####
#Save graphs
setwd("C:/Users/milab/OneDrive/Documentos/UFMG/-Dados/bgest_data") #choose.dir(caption="Select folder to save graphs:"))
for (j in 1:length(plots)) {
  jpeg(paste("graph", j, ".jpg"), width = 1000, height = 600)
  print(plots[[j]])
  dev.off()
}

#Save tables
write_csv(results, "results.csv")
write_csv(gestures, "gestures.csv")
write_csv(utterances, "utterances.csv")

### taaaaam final countdown
results_total <- results %>%
  mutate(ge_phrase_dur = ge_phrase_end-ge_phrase_begin,
         iu_dur = iu_end-iu_begin) %>%
  select("stroke_begin", "stroke_end",
         "ge_phrase_begin", "ge_phrase_end", "ge_phrase_dur",
         "overlap", "iu_begin", "iu_end","iu_dur", "overlapping_time", "interval_time",
         "utt_id",    "iu_id",  "diff_begin", "diff_end",
         "ge_phrase_dur", "overlapping_ratio") %>%
  droplevels()

disp_sd <- results_total %>%
  summarize_if(is.numeric, sd) 

disp_mean <- results_total %>%
  summarize_if(is.numeric, mean)

disp_median <- results_total %>%
  summarize_if(is.numeric, median)

disp_min <- results_total %>%
  summarize_if(is.numeric, min)

disp_max <- results_total %>%
  summarize_if(is.numeric, max)

(dispersion_all <- rbind("sd" = disp_sd, "mean" = disp_mean, "median" = disp_median, "min" = disp_min, "max" = disp_max))

write.table(dispersion %>%
              select("ge_phrase_dur", "iu_dur", "overlapping_time",
                     "diff_begin", "diff_end", 
                     "overlapping_ratio_gpr", "overlapping_ratio_str"), file = "dipersion_total.txt", sep = ";")

#### Clean workspace ####
rm(list=setdiff(ls(), c("results", "gestures", "utterances", "plots", "dispersion_all",
                        "check_overlapping", "get_strokes", "get_closest",
                        "gestures_raw", "utterances_raw")))

