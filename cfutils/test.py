import ssw

aligner = ssw.Aligner()
alignment = aligner.align(
    reference="AAAAACGCCTCCCTCGCGCCATCAGTTTT",
    query="GGGGGCGCCTCCCTCGCGCCATCAGAAATTT",
)
print(alignment.alignment_report())
print(dir(alignment))
