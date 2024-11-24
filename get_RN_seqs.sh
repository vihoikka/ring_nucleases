# This script extracts all RN protein sequences from the path ~/scratch/private/runs/ring_nucleases/run1/backup_150824/71_ring_nucleases_analysis/*/<rn>.faa
# The RNs are crn1, crn2, crn3, csx15, csx16, csx20

# Extract RN protein sequences
for rn in crn1 crn2 crn3 csx15 csx16 csx20
do
    cat ~/scratch/private/runs/ring_nucleases/run1/backup_150824/71_ring_nucleases_analysis/*/$rn.faa >> $rn.faa
done

#take each created faa file and add the RN name to the header as prefix followed by _
for rn in crn1 crn2 crn3 csx15 csx16 csx20
do
    sed -i "s/^>/>$rn\_/" $rn.faa
done