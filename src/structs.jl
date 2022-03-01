"""
QTL parameters.
- **pos**: QTL location indices of the SNPs
- **effect**: effects of allele `1`, of each QTL
"""
struct QTL
    pos
    effect
end
