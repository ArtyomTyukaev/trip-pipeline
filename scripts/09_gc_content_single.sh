#!/bin/bash

filename="$1"

echo "Processing file: $filename"
# Increase Coordinates by 500 bp Upstream and Downstream
echo "Increasing Coordinates by 500 bp Upstream and Downstream"
bedtools slop -b 5000 -g /home/projects/ulitsky/shared/ucsctracks/hg38/hg38.chrom.sizes -i "$filename.bed" > "${filename}_tmp_extended.bed"

# Extract FASTA Sequences with
echo "Extracting FASTA "
bedtools getfasta -fi /home/projects/ulitsky/lidiag/genomic_data/hg38.fa -bed "${filename}_tmp_extended.bed" -name | \
awk '/^>/ {sub(/::.*/, "", $0); print; next} {print}' > "${filename}_sequences.fasta"

# Deleting if exist
rm "${filename}_gc_content.txt"

#Checking that fasta seq in one line
echo "Number of seq in not one line" | grep -v "^>" "${filename}_sequences.fasta" | awk 'length($0) > 0 && length($0) <= 60' | wc -l

# Calculate GC Content
echo "Calculating GC Content"
while read -r line; do
    if [[ $line == ">"* ]]; then
        name=$(echo $line | sed 's/>//')
    else
        gc_count=$(echo $line | grep -o "[GCgc]" | wc -l)
        total_count=$(echo -n $line | wc -c)
        gc_content=$(echo "scale=2; ($gc_count/$total_count)*100" | bc)

        echo -e "$name\t$gc_content" >> "${filename}_gc_content.txt"
    fi
done < "${filename}_sequences.fasta"

paste "$filename.bed" "${filename}_gc_content.txt" > "${filename}_with_gc.bed";

echo "Checking discrepancies between the 4th and 7th columns"
awk '$4 != $7' "${filename}_with_gc.bed"


#Deleting tmp files
#rm "${filename}_tmp_extended.bed" "${filename}_sequences.fasta" "${filename}_gc_content.txt"

echo "Done processing $filename"
echo "--------------------------"
