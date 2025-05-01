process tabix {

    input:
    tuple val(meta), path(vcfz)

    output:
    tuple val(meta), path("*.tbi")

    script:
    """
    tabix $vcfz
    """
}