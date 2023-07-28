process FASTQC {
  publishDir "${label}_fastqc_out"

  input:
    tuple path(read_1), path(read_2)
    val(label)
  
  output:
    path("*.zip"), emit: fastqc_zip
    path("*.html"), emit: fastqc_html
  
  script:
  """
  fastqc ${read_1} ${read_2}
  """
}