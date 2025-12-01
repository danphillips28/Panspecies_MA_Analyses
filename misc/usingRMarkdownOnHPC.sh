#!/bin/bash

# Rmarkdown needs pandoc to render. Need to install locally because no pandoc module
#module load Anaconda3
#conda create -n pandoc_env -c conda-forge pandoc

# Another thing: Safer to install packages locally. Requires setting up local packages directory where R expects it. 
# In R I typed these to find and make that directory
#dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)
#Sys.getenv("R_LIBS_USER") # Just to check what it is
# [1] "/home/ocdm0351/R/x86_64-pc-linux-gnu-library/4.4"
# Can set this to be the lib path in my Rprofiler so packages are automatically installed there:

RPROFILE="$HOME/.Rprofile"
LINE='.libPaths(Sys.getenv("R_LIBS_USER"))'

# Create the file if it doesn't exist
if [ ! -f "$RPROFILE" ]; then
    echo "Creating $RPROFILE"
    echo "$LINE" > "$RPROFILE"
else
    # Check if the line already exists
    if ! grep -Fxq "$LINE" "$RPROFILE"; then
        echo "Adding line to $RPROFILE"
        echo "$LINE" >> "$RPROFILE"
    else
        echo "Line already exists in $RPROFILE"
    fi
fi

