# VirusIntegrate

Detection of virus integration sites in tumor samples

#Song Cao

Usage: perl VirusIntegrate.pl <run_folder> <step_number>

<run_folder> = full path of the folder holding files for this sequence run

<step_number> run this pipeline step by step. (running the whole pipeline if step number is 0)

[1]  Run bwa realignment using paired-end information againt human+virus genome

[2]  Run annotation for fusion site
