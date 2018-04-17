# Load Simulation Results Functions ---------------------------------------

# Functions to load simulation results from different simulation engines
# Florian Antretter
# Date: 2018-04-12
# Version: 0.1


# Function to Load WUFI Plus Results Data ---------------------------------

f.GetWPData <- function(v.Path)
{
  require(readr)
  require(lubridate)
  Count.Col            <- as.numeric(strsplit(as.character(read.table(v.Path, nrows=1, skip=4)[[4]]), ":"))
  Col.Names            <- as.character((read.table(v.Path, nrows=Count.Col, skip=5, sep="\t"))[,1])
  
  Dist.Col             <- c(24, rep(14,Count.Col-1))
  
  de_locale            <- locale(date_names = "de", date_format = "%d.%m.%Y %H:%M:%S", time_format = "%H:%M:%S",
                                 decimal_mark = ",", grouping_mark = "", tz = "UTC", encoding = "UTF-8", asciify = FALSE)
  
  myLocale             <- locale(date_names = "en", date_format = "%m/%d/%Y", time_format = "%H", 
                                 decimal_mark = ".", tz = "UTC", encoding = "UTF-8", asciify = FALSE)
  
  v.WufiData           <- read_fwf(file=v.Path, col_positions = fwf_widths(Dist.Col, col_names = Col.Names), 
                                   skip=6+Count.Col, col_types = paste(c("c", "n", rep("d", Count.Col-2)), sep="", collapse=""), locale=de_locale)
  
  v.WufiData[[1]]      <- dmy_hms(v.WufiData[[1]])
  return(v.WufiData)
}