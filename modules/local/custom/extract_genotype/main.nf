process EXTRACT_GENOTYPE {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(csv)

    output:
    tuple val(meta), stdout, emit: genotype

    script:
    """
    col=\$(awk -F';' 'NR==1{for(i=1;i<=NF;i++){if(\$i=="clade"){print i; exit}}}' "$csv")
    awk -F';' -v c="\$col" 'NR==2{print \$c; exit}' "$csv"
    """

    stub:
    """
    echo "default"
    """
}
