if(dir.exists(here("Output")) == FALSE) {   # Create an Output directory?
  dir.create(here("Output"))
}

if(dir.exists(here("Output")) == TRUE) {
  write(Sys.info(),                         # Store system information
        file = here("Output/System Info.txt"))
  write(Sys.getenv(),                       # Store environment variables
        file = here("Output/Environment.txt"))
}

if(dir.exists(here("Images")) == FALSE) {   # Create an Images directory?
  dir.create(here("Images"))
}
