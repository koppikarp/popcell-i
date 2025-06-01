# popcell-i
current workflow
- select positions from scores csv
- run ensembl on gene list
- run annotate with selected positions, cds summary, and exons json file

note:
original popcell just tries to insert tag as close as possible to termini
popcell-i will only generate pegrna if tag can be inserted at exact location
therefore popcell-i is much more strict
(then)
- find pam