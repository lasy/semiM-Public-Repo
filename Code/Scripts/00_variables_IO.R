
IO = list()

# INPUT AND OUTPUT DATA
source("../../semiM-Restricted-Access-Repo/Scripts/00_variables_restricted_IO.R") # assign value to IO$r_Data

if(par$data_type == "synthetic"){
  folder_name = "synthetic/"
}else if(par$data_type ==  "subset"){
  folder_name = "Kindara_subset/"
}else{
  folder_name = "Kindara/"
}

IO$input_data = paste0(IO$r_Data,folder_name,"input_data/")
IO$tmp_data = paste0(IO$r_Data,folder_name, "tmp_data/")
IO$output_data = paste0(IO$r_Data,folder_name,"output_data/")
if(!dir.exists(IO$tmp_data)){dir.create(IO$tmp_data, recursive = TRUE)}
if(!dir.exists(IO$output_data)){dir.create(IO$output_data, recursive = TRUE)}


# PUBLIC OUTPUT FIGURES AND DATA

IO$panels = paste0("../Figures Tables Media/Figures/panels/",folder_name)
if(!dir.exists(IO$panels)){dir.create(IO$panels, recursive = TRUE)}

IO$tables = paste0("../Figures Tables Media/Tables/",folder_name)
if(!dir.exists(IO$tables)){dir.create(IO$tables, recursive = TRUE)}

IO$public_output_data = paste0("../Data/",folder_name)
IO$out_Rdata = paste0(IO$public_output_data, "Rdata/")
IO$out_csv = paste0(IO$public_output_data, "CSV/")

if(!dir.exists(IO$public_output_data)){dir.create(IO$public_output_data, recursive = TRUE)}
if(!dir.exists(IO$out_Rdata)){dir.create(IO$out_Rdata)}
if(!dir.exists(IO$out_csv)){dir.create(IO$out_csv)}


# Data from the FAM paper
IO$FAM_paper_kindara_data = "../../../FAM/FAM-Restricted-Access-Repo/Data/Kindara/"