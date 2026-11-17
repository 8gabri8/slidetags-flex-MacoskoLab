version 1.0

# Slide-tags positioning for 10x Flex, on Terra.
#   1. list_lines: which lines of the sample sheet have a spatial library (sb_index)
#   2. run_line:   one VM per line, all in parallel; runs slidetags_flex.py --line N
# Each run_line downloads what its line needs and writes the results to
# <slidetags_output_dir>/<Link>_<sb_index>/ itself. summary.csv files are also returned to Terra.

workflow slidetags_flex {
  input {
    File sample_sheet                   # your Cumulus sample sheet + sb_index, pucks, cellranger_output_dir, slidetags_output_dir
    String docker                       # the image built from the Dockerfile
    String barcode_length = "auto"      # Cell Ranger barcode length (without -1): auto, 16, 24
    Int knn = 36                        # eps = median distance from a bead to its knn-th nearest bead
    Float cmes = 10.0                   # clusters closer than eps*cmes count as one
    Int cpu = 8
    String memory = "32G"
    Int disk_gb = 200                   # must hold the spatial FASTQs of one line
    Int preemptible = 1
  }

  call list_lines { input: sample_sheet = sample_sheet, docker = docker }

  scatter (line in list_lines.lines) {
    call run_line {
      input: sample_sheet = sample_sheet, line = line, docker = docker, barcode_length = barcode_length,
             knn = knn, cmes = cmes, cpu = cpu, memory = memory, disk_gb = disk_gb, preemptible = preemptible
    }
  }

  output {
    Array[File] summaries = run_line.summary
  }
}

task list_lines {
  input {
    File sample_sheet
    String docker
  }
  command <<<
    python3 /opt/slidetags_flex.py "~{sample_sheet}" --count-lines > lines.txt
  >>>
  output {
    Array[String] lines = read_lines("lines.txt")
  }
  runtime {
    docker: docker
    cpu: 1
    memory: "2G"
    disks: "local-disk 10 HDD"
    preemptible: 3
  }
}

task run_line {
  input {
    File sample_sheet
    String line
    String docker
    String barcode_length
    Int knn
    Float cmes
    Int cpu
    String memory
    Int disk_gb
    Int preemptible
  }
  command <<<
    set -euo pipefail
    python3 /opt/slidetags_flex.py "~{sample_sheet}" --line ~{line} --barcode-length ~{barcode_length} \
      --knn ~{knn} --cmes ~{cmes} --cores ~{cpu} --workdir work
    cp work/*/results/summary.csv summary.csv
  >>>
  output {
    File summary = "summary.csv"
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: memory
    disks: "local-disk ~{disk_gb} SSD"
    preemptible: preemptible
  }
}
