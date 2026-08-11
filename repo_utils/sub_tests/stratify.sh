# ------------------------------------------------------------
#                                 stratify
# ------------------------------------------------------------
run stratify $truv stratify \
        $INDIR/beds/include.bed \
        $INDIR/variants/input1.vcf.gz \
        -o $OD/stratify.txt
if [ $stratify ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/stratify/stratify.txt) $(fn_md5 $OD/stratify.txt)
fi

# Regions on purely numeric contigs must not be inferred as integers
run stratify_numeric_chrom $truv stratify \
        $INDIR/beds/numeric_chrom.bed \
        $INDIR/variants/input1.vcf.gz \
        -o $OD/stratify_numeric_chrom.txt
if [ $stratify_numeric_chrom ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/stratify/stratify_numeric_chrom.txt) $(fn_md5 $OD/stratify_numeric_chrom.txt)
fi
