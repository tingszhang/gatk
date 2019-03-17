package org.broadinstitute.hellbender.utils.genotyper;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.Object2IntArrayMap;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import org.apache.commons.collections.ListUtils;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.hellbender.utils.fragments.FragmentCollection;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Read-likelihoods container implementation based on integer indexed arrays.
 *
 * @param <A> the type of the allele the likelihood makes reference to.
 *
 * Note: this class uses FastUtil collections for speed.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class ReadLikelihoods<A extends Allele> extends MoleculeLikelihoods<GATKRead, A> {
    /**
     * Constructs a new read-likelihood collection.
     *
     * <p>
     *     The initial likelihoods for all allele-read combinations are
     *     0.
     * </p>
     *
     * @param samples all supported samples in the collection.
     * @param alleles all supported alleles in the collection.
     * @param reads reads stratified per sample.
     *
     * @throws IllegalArgumentException if any of {@code allele}, {@code samples}
     * or {@code reads} is {@code null},
     *  or if they contain null values.
     */
    @SuppressWarnings({"rawtypes", "unchecked"})
    public ReadLikelihoods(final SampleList samples,
                               final AlleleList<A> alleles,
                               final Map<String, List<GATKRead>> reads) {
        super(samples, alleles, reads);
    }

    // Internally used constructor.
    @SuppressWarnings({"unchecked", "rawtypes"})
    ReadLikelihoods(final AlleleList alleles,
                        final SampleList samples,
                        final List<List<GATKRead>> readsBySampleIndex,
                        final Object2IntMap<GATKRead>[] readIndex,
                        final double[][][] values) {
        super(alleles, samples, readsBySampleIndex, readIndex, values);
    }

    /**
     * Removes those read that the best possible likelihood given any allele is just too low.
     *
     * <p>
     *     This is determined by a maximum error per read-base against the best likelihood possible.
     * </p>
     *
     * @param maximumErrorPerBase the minimum acceptable error rate per read base, must be
     *                            a positive number.
     *
     * @throws IllegalStateException is not supported for read-likelihood that do not contain alleles.
     *
     * @throws IllegalArgumentException if {@code maximumErrorPerBase} is negative.
     */
    public void filterPoorlyModeledReads(final double maximumErrorPerBase) {
        Utils.validateArg(alleles.numberOfAlleles() > 0, "unsupported for read-likelihood collections with no alleles");
        Utils.validateArg(!Double.isNaN(maximumErrorPerBase) && maximumErrorPerBase > 0.0, "the maximum error per base must be a positive number");

        new IndexRange(0, samples.numberOfSamples()).forEach(s -> {
            final List<GATKRead> sampleReads = readsBySampleIndex.get(s);
            final List<GATKRead> readsToRemove = IntStream.range(0, sampleReads.size())
                    .filter(r -> readIsPoorlyModelled(s, r, sampleReads.get(r), maximumErrorPerBase))
                    .mapToObj(sampleReads::get)
                    .collect(Collectors.toList());

            removeSampleReads(s, readsToRemove, alleles.numberOfAlleles());
        });
    }

    private boolean readIsPoorlyModelled(final int sampleIndex, final int readIndex, final GATKRead read, final double maxErrorRatePerBase) {
        final double maxErrorsForRead = Math.min(2.0, Math.ceil(read.getLength() * maxErrorRatePerBase));
        final double log10QualPerBase = -4.0;
        final double log10MaxLikelihoodForTrueAllele = maxErrorsForRead * log10QualPerBase;

        final int alleleCount = alleles.numberOfAlleles();
        final double[][] sampleValues = valuesBySampleIndex[sampleIndex];
        for (int a = 0; a < alleleCount; a++) {
            if (sampleValues[a][readIndex] >= log10MaxLikelihoodForTrueAllele) {
                return false;
            }
        }
        return true;
    }

    /**
     * Downsamples reads based on contamination fractions making sure that all alleles are affected proportionally.
     *
     * @param perSampleDownsamplingFraction contamination sample map where the sample name are the keys and the
     *                                       fractions are the values.
     *
     * @throws IllegalArgumentException if {@code perSampleDownsamplingFraction} is {@code null}.
     */
    public void contaminationDownsampling(final Map<String, Double> perSampleDownsamplingFraction) {
        Utils.nonNull(perSampleDownsamplingFraction);

        final int alleleCount = alleles.numberOfAlleles();
        for (int s = 0; s < samples.numberOfSamples(); s++) {
            final String sample = samples.getSample(s);
            final Double fractionDouble = perSampleDownsamplingFraction.get(sample);
            if (fractionDouble == null) {
                continue;
            }
            final double fraction = fractionDouble;
            if (Double.isNaN(fraction) || fraction <= 0.0) {
                continue;
            }
            if (fraction >= 1.0) {
                removeSampleReads(s, readsBySampleIndex.get(s), alleleCount);
            } else {
                final Map<A,List<GATKRead>> readsByBestAllelesMap = readsByBestAlleleMap(s);
                removeSampleReads(s, AlleleBiasedDownsamplingUtils.selectAlleleBiasedReads(readsByBestAllelesMap, fraction),alleleCount);
            }
        }
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    public MoleculeLikelihoods<Fragment, A> combineMates() {
        final int sampleCount = samples.numberOfSamples();
        final double[][][] newLikelihoodValues = new double[sampleCount][][];
        final int alleleCount = alleles.numberOfAlleles();

        final Object2IntMap<Fragment>[] fragmentIndexBySampleIndex = new Object2IntMap[sampleCount];
        final List<List<Fragment>> fragmentsBySampleIndex = new ArrayList<>(sampleCount);

        for (int s = 0; s < sampleCount; s++) {
            final FragmentCollection<GATKRead> fragmentCollection = FragmentCollection.create(sampleReads(s).stream().sorted(Comparator.comparingInt(GATKRead::getStart)).collect(Collectors.toList()));

            final List<Fragment> sampleFragments = Stream.concat(fragmentCollection.getOverlappingPairs().stream().map(Fragment::new),
                    fragmentCollection.getSingletonReads().stream().map(Fragment::new)).collect(Collectors.toList());

            final int sampleFragmentCount = sampleFragments.size();

            final double[][] oldSampleValues = valuesBySampleIndex[s];
            final double[][] newSampleValues = newLikelihoodValues[s] = new double[alleleCount][sampleFragmentCount];

            // For each old allele and read we update the new table keeping the maximum likelihood.
            fragmentIndexBySampleIndex[s] = new Object2IntArrayMap<>(sampleFragmentCount);
            for (int f = 0; f < sampleFragmentCount; f++) {
                fragmentIndexBySampleIndex[s].put(sampleFragments.get(f), f);
                for (int a = 0; a < alleleCount; a++) {
                    for (final GATKRead read : sampleFragments.get(f).getReads()) {
                        final int oldReadIndex = readIndex(s, read);
                        newSampleValues[a][f] += oldSampleValues[a][oldReadIndex];
                    }
                }
            }
            fragmentsBySampleIndex.add(sampleFragments);
        }

        // Finally we create the new read-likelihood
        return new MoleculeLikelihoods<Fragment, A>(
                alleles,
                samples,
                fragmentsBySampleIndex,
                fragmentIndexBySampleIndex, newLikelihoodValues);
    }


    /**
     * Create an independent copy of this read-likelihoods collection
     */
    @VisibleForTesting
    public ReadLikelihoods<A> copy() {

        final int sampleCount = samples.numberOfSamples();
        final int alleleCount = alleles.numberOfAlleles();

        final double[][][] newLikelihoodValues = new double[sampleCount][alleleCount][];

        @SuppressWarnings({"unchecked", "rawtypes"})
        final Object2IntMap<GATKRead>[] newReadIndexBySampleIndex = new Object2IntMap[sampleCount];
        final List<List<GATKRead>> newReadsBySampleIndex = new ArrayList<>(sampleCount);

        for (int s = 0; s < sampleCount; s++) {
            newReadsBySampleIndex.add(new ArrayList<>(readsBySampleIndex.get(s)));
            for (int a = 0; a < alleleCount; a++) {
                newLikelihoodValues[s][a] = valuesBySampleIndex[s][a].clone();
            }
        }

        // Finally we create the new read-likelihood
        return new ReadLikelihoods<>(
                alleles,
                samples,
                newReadsBySampleIndex,
                newReadIndexBySampleIndex,
                newLikelihoodValues);
    }

    public ReadLikelihoods(final MoleculeLikelihoods<GATKRead, A> moleculeLikelihoods) {
        this(moleculeLikelihoods.alleles, moleculeLikelihoods.samples, moleculeLikelihoods.readsBySampleIndex,
                moleculeLikelihoods.readIndexBySampleIndex, moleculeLikelihoods.valuesBySampleIndex);
    }

    @Override
    @SuppressWarnings({"unchecked", "rawtypes"})
    public <B extends Allele> ReadLikelihoods<B> marginalize(final Map<B, List<A>> newToOldAlleleMap) {
        return new ReadLikelihoods<B>(super.marginalize(newToOldAlleleMap));
    }

    @Override
    @SuppressWarnings({"unchecked", "rawtypes"})
    public <B extends Allele> ReadLikelihoods<B> marginalize(final Map<B, List<A>> newToOldAlleleMap, final Locatable overlap) {
        return new ReadLikelihoods<B>(super.marginalize(newToOldAlleleMap, overlap));
    }
}