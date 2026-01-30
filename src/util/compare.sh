#
#fileA=$1
#fileB=$2
#intersections=$(bedtools intersect -u -a $fileA -b $fileB | wc -l)
#
## Total intervals in both sets
#total_a=$(wc -l < $fileA)
#total_b=$(wc -l < $fileB)
#union=$((total_a + total_b - intersections))
#
#echo "scale=4; $intersections / $union" | bc


fileA=$1
fileB=$2

# Number of overlapping intervals between A and B
intersections=$(bedtools intersect -u -a "$fileA" -b "$fileB" | wc -l)

# Create a true union of intervals from both files
#union_intervals=$(bedtools cat -i "$fileA" "$fileB" |  bedtools merge -i - | wc -l)

union_intervals=$(cat "$fileA" "$fileB" \
                  | sort -k1,1 -k2,2n \
                  | bedtools merge -i - \
                  | wc -l)


# Compute the Jaccard-like ratio using total *base pair* coverage, not just interval counts
echo "scale=4; $intersections / $union_intervals" | bc
#echo $union_intervals

