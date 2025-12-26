#snakemake \
#  -s workflow/Snakefile \
#  -j 4 \
#  --cores 2 \
#  --rerun-incomplete \
#  --printshellcmds

clear

snakemake \
      -s workflow/Snakefile \
      --jobs 4 \
      --cores 2 \
      --keep-going \
      --latency-wait 60 \
      --unlock

snakemake \
      -s workflow/Snakefile \
      --jobs 4 \
      --cores 2 \
      --keep-going \
      --latency-wait 60 \
      --printshellcmds \
      --rerun-incomplete
      
# snakemake -s Snakefile.updated2 -j 32 --rerun-incomplete --printshellcmds
