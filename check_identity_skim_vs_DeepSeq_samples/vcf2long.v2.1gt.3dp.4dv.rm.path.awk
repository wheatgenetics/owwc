#!/bin/awk -f

BEGIN{
 OFS="\t"
}

/^#CHROM/ && !header {
 for(i = 1; i <= NF; i++){
  x=$i
  sub(".+/", "", x)
  sub("_sort_rmdup.bam", "", x)
  s[i]=x;
  }
  header=1
}

/^#/ { 
 next
}


{
 for(i = 10; i <= NF; i++){
  split($i, a, ":")
  if (length($4)==1 && length($5)==1) { print $1, $2, $4, $5, $6, s[i], a[1], a[3], a[4] }
 }
}
