"""
One off unittests
"""
import os
import sys
import unittest

from collections import defaultdict
from intervaltree import IntervalTree

# Use the current truvari, not any installed libraries
sys.path.insert(0, os.getcwd())
import truvari
from truvari import phab
from truvari.region_vcf_iter import (region_filter_stream, region_filter_fetch,
                                     region_filter_stream_overlap)

# Assume we're running in truvari root directory


class TestBoundaries(unittest.TestCase):
    def test_boundary_issues(self):
        """
        Test if variants are correctly classified as within or outside a boundary.
        """
        vcf_fn = "repo_utils/test_files/variants/boundary.vcf.gz"
        region_start = 10
        region_end = 20

        vcf = truvari.VariantFile(vcf_fn)
        for entry in vcf:
            # Separate failure reporting per variant
            with self.subTest(entry=str(entry)):
                state = entry.info['include'] == 'in'
                self.assertEqual(
                    state, entry.within(region_start, region_end),
                    f"Bad Boundary {str(entry)}"
                )

    def test_region_filtering(self):
        """
        Test if regions are correctly filtered based on the provided BED file.
        """
        vcf_fn = "repo_utils/test_files/variants/boundary_cpx.vcf.gz"
        bed_fn = "repo_utils/test_files/beds/boundary_cpx.bed"

        tree = defaultdict(IntervalTree)
        with open(bed_fn, 'r') as fh:
            for line in fh:
                data = line.strip().split()
                tree[data[0]].addi(int(data[1]), int(data[2]) + 1)

        vcf = truvari.VariantFile(vcf_fn)

        # Using fetch_regions()
        for entry in vcf.fetch_regions(tree, True, False):
            with self.subTest(method="fetch_regions in", entry=str(entry)):
                self.assertEqual(
                    entry.info['include'], 'in', f"Bad in {str(entry)}")

        vcf.reset()
        for entry in vcf.fetch_regions(tree, False, False):
            with self.subTest(method="fetch_regions out", entry=str(entry)):
                self.assertEqual(
                    entry.info['include'], 'out', f"Bad out {str(entry)}")

        # Using region_filter_stream()
        vcf = truvari.VariantFile(vcf_fn)
        for entry in region_filter_stream(vcf, tree, True, False):
            with self.subTest(method="region_filter_stream in", entry=str(entry)):
                self.assertEqual(
                    entry.info['include'], 'in', f"Bad in {str(entry)}")

        vcf.reset()
        for entry in region_filter_stream(vcf, tree, False, False):
            with self.subTest(method="region_filter_stream out", entry=str(entry)):
                self.assertEqual(
                    entry.info['include'], 'out', f"Bad out {str(entry)}")

        # Using region_filter_fetch()
        vcf = truvari.VariantFile(vcf_fn)
        for entry in region_filter_fetch(vcf, tree, False):
            with self.subTest(method="region_filter_fetch in", entry=str(entry)):
                self.assertEqual(
                    entry.info['include'], 'in', f"Bad in {str(entry)}")


class TestOverlapRegionFiltering(unittest.TestCase):
    """
    Overlap region filtering keeps entries intersecting a region instead of
    entries contained by one.
    """
    vcf_fn = "repo_utils/test_files/variants/boundary_cpx.vcf.gz"
    bed_fn = "repo_utils/test_files/beds/boundary_cpx.bed"

    def setUp(self):
        self.tree = defaultdict(IntervalTree)
        with open(self.bed_fn, 'r') as fh:
            for line in fh:
                data = line.strip().split()
                self.tree[data[0]].addi(int(data[1]), int(data[2]) + 1)

    def vcf(self):
        return truvari.VariantFile(self.vcf_fn)

    @staticmethod
    def keys(entries):
        return sorted((e.chrom, e.pos, e.ref, e.alts[0]) for e in entries)

    def bruteforce(self, minsize):
        """
        Every entry sharing minsize positions with a region, tested directly
        """
        ret = []
        for entry in self.vcf():
            qstart, qend = entry.boundaries()
            needed = max(1, min(minsize, qend - qstart))
            if any(truvari.overlap_size(qstart, qend, i.begin, i.end - 1) >= needed
                   for i in self.tree[entry.chrom]):
                ret.append(entry)
        return self.keys(ret)

    def test_agrees_with_bruteforce(self):
        """
        The filter emits each qualifying entry exactly once
        """
        for minsize in (1, 4, 7, 1000):
            with self.subTest(minsize=minsize):
                self.assertEqual(self.bruteforce(minsize),
                                 self.keys(region_filter_stream_overlap(
                                     self.vcf(), self.tree, minsize=minsize)))

    def test_superset_of_containment(self):
        """
        Anything contained by a region necessarily intersects it
        """
        contained = set(self.keys(
            region_filter_stream(self.vcf(), self.tree, True, False)))
        overlapping = set(self.keys(
            region_filter_stream_overlap(self.vcf(), self.tree)))
        self.assertTrue(contained.issubset(overlapping),
                        f"contained entries dropped: {contained - overlapping}")
        self.assertTrue(overlapping - contained,
                        "fixture has no boundary spanning entries")

    def test_minsize_narrows_but_keeps_short_entries(self):
        """
        Raising minsize can only drop entries, never entries too short to reach it
        """
        counts = [len(self.keys(region_filter_stream_overlap(
            self.vcf(), self.tree, minsize=m))) for m in (1, 5, 7, 1000)]
        self.assertEqual(counts, sorted(counts, reverse=True))
        self.assertGreater(counts[0], counts[-1], "minsize never narrowed")
        self.assertGreater(counts[-1], 0, "one position entries were dropped")

    def test_fetch_regions_overlap(self):
        """
        VariantFile.fetch_regions routes overlap to the overlap filter
        """
        self.assertEqual(
            self.keys(region_filter_stream_overlap(self.vcf(), self.tree, minsize=3)),
            self.keys(self.vcf().fetch_regions(self.tree, overlap=3)))
        with self.assertRaises(ValueError):
            self.vcf().fetch_regions(self.tree, inside=False, overlap=1)

    def test_with_region(self):
        """
        The returned region is one the entry actually intersects
        """
        for entry, (chrom, intv) in region_filter_stream_overlap(
                self.vcf(), self.tree, with_region=True):
            self.assertEqual(entry.chrom, chrom)
            qstart, qend = entry.boundaries()
            self.assertGreater(
                truvari.overlap_size(qstart, qend, intv.begin, intv.end - 1), 0)

    def test_overlaps_tree(self):
        """
        VariantRecord.overlaps_tree applies the same rule as the filter
        """
        for minsize in (1, 7):
            with self.subTest(minsize=minsize):
                self.assertEqual(self.bruteforce(minsize),
                                 self.keys(e for e in self.vcf()
                                           if e.overlaps_tree(self.tree, minsize)))
        self.assertFalse(next(iter(self.vcf())).overlaps_tree(defaultdict(IntervalTree)))


class TestFilteringLogic(unittest.TestCase):
    """
    Filtering logic
    """

    def test_filtering_logic(self):
        p = truvari.VariantParams(sizemin=0, sizefilt=0, passonly=True)
        vcf = truvari.VariantFile(
            "repo_utils/test_files/variants/filter.vcf", params=p)

        for entry in vcf:
            # Ensures all entries are tested separately
            with self.subTest(entry=str(entry)):
                try:
                    self.assertTrue(entry.filter_call(),
                                    f"Didn't filter {str(entry)}")
                except ValueError as e:
                    self.assertTrue(
                        e.args[0].startswith("Cannot compare multi-allelic"),
                        f"Unknown exception {str(entry)}"
                    )


class TestPhab(unittest.TestCase):
    """
    Manually running phab to ensure coverage. This doesn't actually check anything,
    just lets coverage.py see the code being executed. We rely on the phab/refine functional sub_tests 
    for accuracy
    """

    def test_phab(self):
        vcfs = [
            "repo_utils/test_files/variants/phab_base.vcf.gz",
            "repo_utils/test_files/variants/phab_comp2.vcf.gz"
        ]
        vcfhap = phab.VCFtoHaplotypes(
            "repo_utils/test_files/references/phab_ref.fa",
            vcfs, ["HG002", "HG005"]
        )

        reg_fn = truvari.make_temp_filename()
        with open(reg_fn, 'w') as fout:
            fout.write("chr1\t700\t900\n")

        regions = phab.parse_regions(reg_fn)
        vcfhap.set_regions(regions, buff=0)
        haps = vcfhap.build_all()
        key = 'chr1:700-900'
        result = phab.align_wrap(haps[key], phab.run_wfa)
        self.assertFalse(result.startswith("ERROR"), f"WFA failed: {result}")

        result = phab.align_wrap(haps[key], phab.run_poa, dedup=True)
        self.assertFalse(result.startswith("ERROR"), f"POA failed: {result}")

        result = phab.run_mafft(haps[key])
        self.assertTrue(isinstance(result, dict), f"MAFFT failed: {result}")

        with self.assertRaises(TypeError):
            phab.align_wrap("", phab.run_wfa)

        with self.assertRaises(ValueError):
            phab.get_align_method("bob")

        d = {'O':'A', 'ref_':'A'}
        self.assertTrue(phab.deduplicate_haps(d)[0] == {'ref_':'A'}, "Phab dedup failed")

        import pickle
        import multiprocessing.shared_memory as shm
        vcf_info = FakeVCFI()
        data = pickle.dumps(vcf_info)
        shared_info = shm.SharedMemory(create=True, size=len(data))
        shared_info.buf[:len(data)] = data
        mem_vcf_info = (shared_info.name, shared_info.size)
        job = phab.PhabJob("fake", mem_vcf_info)
        self.assertTrue(phab.status_marker(1, {}, "", job) == "ERROR: 'str' object is not callable", "PhabJob")
        phab.cleanup_shared_memory(shared_info)


class FakeVCFI():
    def __init__(self):
        pass
    def get_haplotypes(self, _):
        return "A"

if __name__ == "__main__":
    unittest.main()
