version 1.0

# Slide-tags positioning for 10x Flex, on Terra. Cell Ranger (Cumulus) must have run before.
#   1. list_lines: which lines of the sample sheet have a spatial library (sb_index), and the sections
#   2. run_line:   one VM per line, all in parallel: steps 1-7 (slidetags_flex.py --line N)
#                  -> <slidetags_output_dir>/<Link>_<sb_index>/
#   3. summarize:  after ALL run_line jobs, one VM per section, in parallel: steps 8-9
#                  (slidetags_flex.py --summarize SECTION): spatial libraries of the same GEM well
#                  pooled, every nucleus placed once, files.csv, coords.csv, PDF
#                  -> <slidetags_output_dir>/summary/<section>/
# only_summarize = true skips run_line and rebuilds the summaries from the run_line results
# already in the bucket (e.g. to redo the PDFs).

# Input names follow Cumulus cellranger_workflow where both have the same setting
# (input_csv_file, num_cpu, memory, preemptible, zones), so the two JSONs look alike.
workflow slidetags_flex {
  input {
    File input_csv_file                 # the SAME sample sheet as Cumulus + sb_index, pucks, cellranger_output_dir, slidetags_output_dir
    String docker                       # the image built from the Dockerfile
    Boolean only_summarize = false      # true: skip run_line, only redo the per-section summaries
    String barcode_length = "auto"      # Cell Ranger barcode length (without -1): auto, 16, 26
    Int knn = 36                        # eps = median distance from a bead to its knn-th nearest bead
    Float cmes = 10.0                   # clusters closer than eps*cmes count as one
    Int num_cpu = 8                     # run_line VM
    String memory = "32G"
    Int disk_space = 100                # GB; must hold the spatial FASTQs of one line
    Int preemptible = 1
    Int summarize_num_cpu = 8           # summarize VM: reads all RNA matrices of one section
    String summarize_memory = "64G"
    Int summarize_disk_space = 100      # GB; Cell Ranger matrices of one section + run_line results
    String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
  }

  call list_lines { input: input_csv_file = input_csv_file, docker = docker, zones = zones }

  if (!only_summarize) {
    scatter (line in list_lines.lines) {
      call run_line {
        input: input_csv_file = input_csv_file, line = line, docker = docker, barcode_length = barcode_length,
               knn = knn, cmes = cmes, num_cpu = num_cpu, memory = memory, disk_space = disk_space,
               preemptible = preemptible, zones = zones
      }
    }
  }

  scatter (section in list_lines.sections) {
    call summarize {
      input: input_csv_file = input_csv_file, section = section, docker = docker, barcode_length = barcode_length,
             knn = knn, cmes = cmes, num_cpu = summarize_num_cpu, memory = summarize_memory,
             disk_space = summarize_disk_space, preemptible = preemptible, zones = zones,
             wait_for = run_line.summary   # summaries start only after every run_line has finished
    }
  }

  output {
    Array[File]? line_summaries = run_line.summary
    Array[File] section_pdfs = summarize.pdf
    Array[File] section_files = summarize.files
    Array[File] section_coords = summarize.coords
  }
}

task list_lines {
  input {
    File input_csv_file
    String docker
    String zones
  }
  command <<<
    set -euo pipefail
    python3 /opt/slidetags_flex.py "~{input_csv_file}" --count-lines > lines.txt
    python3 /opt/slidetags_flex.py "~{input_csv_file}" --list-sections > sections.txt
  >>>
  output {
    Array[String] lines = read_lines("lines.txt")
    Array[String] sections = read_lines("sections.txt")
  }
  runtime {
    docker: docker
    zones: zones
    cpu: 1
    memory: "2G"
    disks: "local-disk 10 HDD"
    preemptible: 3
  }
}

task run_line {
  input {
    File input_csv_file
    String line
    String docker
    String barcode_length
    Int knn
    Float cmes
    Int num_cpu
    String memory
    Int disk_space
    Int preemptible
    String zones
  }
  command <<<
    set -euo pipefail
    python3 /opt/slidetags_flex.py "~{input_csv_file}" --line ~{line} --barcode-length ~{barcode_length} \
      --knn ~{knn} --cmes ~{cmes} --cores ~{num_cpu} --workdir work
    cp work/*/results/summary.csv summary.csv
  >>>
  output {
    File summary = "summary.csv"
  }
  runtime {
    docker: docker
    zones: zones
    cpu: num_cpu
    memory: memory
    disks: "local-disk ~{disk_space} SSD"
    preemptible: preemptible
  }
}

task summarize {
  input {
    File input_csv_file
    String section
    String docker
    String barcode_length
    Int knn
    Float cmes
    Int num_cpu
    String memory
    Int disk_space
    Int preemptible
    String zones
    Array[File]? wait_for
  }
  command <<<
    set -euo pipefail
    python3 /opt/slidetags_flex.py "~{input_csv_file}" --summarize "~{section}" --barcode-length ~{barcode_length} \
      --knn ~{knn} --cmes ~{cmes} --cores ~{num_cpu} --workdir work
    cp "work/summary/~{section}/~{section}_summary.pdf" "~{section}_summary.pdf"
    cp "work/summary/~{section}/files.csv" files.csv
    cp "work/summary/~{section}/coords.csv" coords.csv
  >>>
  output {
    File pdf = "~{section}_summary.pdf"
    File files = "files.csv"
    File coords = "coords.csv"
  }
  runtime {
    docker: docker
    zones: zones
    cpu: num_cpu
    memory: memory
    disks: "local-disk ~{disk_space} SSD"
    preemptible: preemptible
  }
}
