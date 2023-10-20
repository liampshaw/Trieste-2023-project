# R-M sites in plasmids

This repository was created during a group project as part of the EMBO workshop on 'Plasmids as vehicles of AMR spread', ICTP, Trieste, September 18-22 2023. 

It is currently a work-in-progress. We investigated the density of palindromes in core and accessory genes of five different plasmid taxonomic units (PTUs). 

# Workflow (work-in-progress)

* Run roary on each PTU using [https://github.com/Adalijuanluo/Plasmid_pan](https://github.com/Adalijuanluo/Plasmid_pan)
* Run `scripts/workflow-50.sh`: first, use roary output and .ffn files from Prokka to obtain core/accessory gene fastas (every gene sequence from every plasmid ends up in either of these two files) with `make_core_accessory_fasta_from_roary.py`. Then also calculate palindrome densities and GC contents.
* Analyse these fastas to get avoidance score of motifs with R'MES - see `run-rmes.sh` 


