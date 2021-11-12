#!/bin/bash
cd $1

find .  -name "*\.zip" -exec unzip -q {} \; >/dev/null # unzip the results files

# - .json files have no stochastic content, may be md5sum-checked
# - .txt files are generically named, some stochastic content, may be md5sum-checked after rounding
# - .seg files are processed similarly
# - .pdf files are NOT PROCESSED, there is some stochastic part besides Date info

# Therefore:
# - Check md5sums for all types of files, sort

echo ".json files:"
find . -name "*.json" | xargs md5sum | sort -V

echo "mutations.txt files:"
for f in $(find -name *_mutations.txt);do awk '{print $1,$2,$3,$4,$5,$7,$8,$9,$10}' $f | md5sum;done | sort -V
echo "segments.txt files:"
for f in $(find -name *_segments.txt);do awk '{printf "%s %i %i %.3f %i %.3f %.3f %i %.3f %i %i%i %.3f\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' $f | md5sum;done | sort -V

echo ".seg files:"
for f in $(find . -name *.seg);do awk '{printf "%s %s %i %i %i %.3f\n", $1, $2, $3, $4, $5, $6}' $f | md5sum;done | sort -V

