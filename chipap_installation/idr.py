import os, sys

import gzip, io

import math

import numpy

from scipy.stats.stats import rankdata

from collections import namedtuple, defaultdict, OrderedDict
from itertools import chain

def mean(items):
    items = list(items)
    if len(items)==0: return 0.0
    return sum(items)/float(len(items))

import idr

import idr.optimization
from idr.optimization import estimate_model_params, old_estimator
from idr.utility import calc_post_membership_prbs, compute_pseudo_values

Peak = namedtuple(
    'Peak', ['chrm', 'strand', 'start', 'stop', 'signal', 'summit', 'signalValue', 'pValue', 'qValue'])
MergedPeak = namedtuple(
    'Peak', ['chrm', 'strand', 'start', 'stop', 'summit', 
             'merged_signal', 'signals', 'pks'])

def load_gff(fp):
    """
    chr20   GRIT    TSS     36322438        36322468        44      +       .       gene_id 'chr20_plus_36322407_36500530'; gene_name 'chr20_plus_36322407_36500530'; tss_id 'TSS_chr20_plus_36322407_36500530_pk1'; peak_cov '7,0,11,0,0,0,0,0,3,0,1,0,0,0,6,0,0,0,0,0,3,0,4,0,0,0,8,0,0,1';
    """
    grpd_peaks = defaultdict(list)
    for line in fp:
        if line.startswith("#"): continue
        if line.startswith("track"): continue
        data = line.split()
        signal = float(data[5])
        peak = Peak(data[0], data[6], 
                    int(float(data[3])), int(float(data[4])), 
                    signal, None, 
                    None, None, None )
        grpd_peaks[(peak.chrm, peak.strand)].append(peak)
    return grpd_peaks

def load_bed(fp, signal_index, peak_summit_index=None):
    grpd_peaks = defaultdict(list)
    for line in fp:
        if line.startswith("#"): continue
        if line.startswith("track"): continue
        data = line.split()
        signal = float(data[signal_index])
        if idr.ONLY_ALLOW_NON_NEGATIVE_VALUES and signal < 0:
            raise ValueError("Invalid Signal Value: {:e}".format(signal))
        if peak_summit_index == None or int(data[peak_summit_index]) == -1:
            summit = None
        else:
            summit = int(data[peak_summit_index]) + int(float(data[1]))
        assert summit == None or summit >= 0
        peak = Peak(data[0], data[5], 
                    int(float(data[1])), int(float(data[2])), 
                    signal, summit, 
                    float(data[6]), float(data[7]), float(data[8]) 
        )
        grpd_peaks[(peak.chrm, peak.strand)].append(peak)
    return grpd_peaks

def correct_multi_summit_peak_IDR_values(idr_values, merged_peaks):
    assert len(idr_values) == len(merged_peaks)
    new_values = idr_values.copy()
    # find the maximum IDR value for each peak
    pk_idr_values = defaultdict(lambda: float('inf')) 
    for i, pk in enumerate(merged_peaks):
        pk_idr_values[(pk.chrm, pk.strand, pk.start, pk.stop)] = min(
            pk_idr_values[(pk.chrm, pk.strand, pk.start, pk.stop)], 
            idr_values[i]
        )
    # store the indices best peak indices, and update the values
    best_indices = []
    for i, pk in enumerate(merged_peaks):
        region = (pk.chrm, pk.strand, pk.start, pk.stop)
        if new_values[i] == pk_idr_values[region]:
            best_indices.append(i)
        else:
            new_values[i] = pk_idr_values[region]
    return numpy.array(best_indices), new_values

def iter_merge_grpd_intervals(
        intervals, n_samples, pk_agg_fn,
        use_oracle_pks, use_nonoverlapping_peaks):
    # grp peaks by their source, and calculate the merged
    # peak boundaries
    grpd_peaks = OrderedDict([(i+1, []) for i in range(n_samples)])
    pk_start, pk_stop = 1e12, -1
    for interval, sample_id in intervals:
        # if we've provided a unified peak set, ignore any intervals that 
        # don't contain it for the purposes of generating the merged list
        if (not use_oracle_pks) or sample_id == 0:
            pk_start = min(interval.start, pk_start)
            pk_stop = max(interval.stop, pk_stop)
        # if this is an actual sample (ie not a merged peaks)
        if sample_id > 0:
            grpd_peaks[sample_id].append(interval)
    
    # if there are no identified peaks, continue (this can happen if 
    # we have a merged peak list but no merged peaks overlap sample peaks)
    if pk_stop == -1:
        return None

    # skip regions that dont have a peak in all replicates
    if not use_nonoverlapping_peaks:
        if any(0 == len(peaks) for peaks in grpd_peaks.values()):
            return None

    # find the merged peak summit
    # note that we can iterate through the values because 
    # grpd_peaks is an ordered dict
    replicate_summits = []
    for sample_id, pks in grpd_peaks.items():
        # if an oracle peak set is specified, skip the replicates
        if use_oracle_pks and sample_id != 0: 
            continue

        # initialize the summit to the first peak
        try: replicate_summit, summit_signal = pks[0].summit, pks[0].signal
        except IndexError: replicate_summit, summit_signal =  None, -1e9
        # if there are more peaks, take the summit that corresponds to the 
        # replicate peak with the highest signal value
        for pk in pks[1:]:
            if pk.summit != None and pk.signal > summit_signal:
                replicate_summit, summit_signal = pk.summit, pk.signal
        # make sure a peak summit was specified
        if replicate_summit != None:
            replicate_summits.append( replicate_summit )

    summit = ( int(mean(replicate_summits)) 
               if len(replicate_summits) > 0 else None )

    # note that we can iterate through the values because 
    # grpd_peaks is an ordered dict
    signals = [pk_agg_fn(pk.signal for pk in pks) if len(pks) > 0 else 0
              for pks in grpd_peaks.values()]
    merged_pk = (pk_start, pk_stop, summit, 
                 pk_agg_fn(signals), signals, grpd_peaks)

    yield merged_pk
    return

def iter_matched_oracle_pks(
        pks, n_samples, pk_agg_fn, use_nonoverlapping_peaks=False ):
    """Match each oracle peak to it nearest replicate peaks.

    """
    oracle_pks = [pk for pk, sample_id in pks
                  if sample_id == 0]
    # if there are zero oracle peaks in this 
    if len(oracle_pks) == 0: return None
    # for each oracle peak, find score the replicate peaks
    for oracle_pk in oracle_pks:
        peaks_and_scores = OrderedDict([(i+1, []) for i in range(n_samples)])
        for pk, sample_id in pks:
            # skip oracle peaks
            if sample_id == 0: continue
            
            # calculate the distance between summits, setting it to a large
            # value in case the peaks dont have summits
            summit_distance = sys.maxsize
            if oracle_pk.summit != None and pk.summit != None:
                summit_distance = abs(oracle_pk.summit - pk.summit)
            # calculate the fraction overlap witht he oracle peak
            overlap = (1 + min(oracle_pk.stop, pk.stop) 
                       - max(oracle_pk.start, pk.start) ) 
            overlap_frac = overlap/(oracle_pk.stop - oracle_pk.start + 1)
            
            peaks_and_scores[sample_id].append(
                ((summit_distance, -overlap_frac, -pk.signal), pk))
                
        # skip regions that dont have a peak in all replicates. 
        if not use_nonoverlapping_peaks and any(
                0 == len(peaks) for peaks in peaks_and_scores.values()):
            continue
        
        # build the aggregated signal value, which is jsut the signal value
        # of the replicate peak witgh the closest match
        signals = []
        rep_pks = []
        for rep_id, scored_pks in peaks_and_scores.items():
            scored_pks.sort()
            if len(scored_pks) == 0:
                assert use_nonoverlapping_peaks
                signals.append(0)
                rep_pks.append(None)
            else:
                signals.append(scored_pks[0][1].signal)
                rep_pks.append( [scored_pks[0][1],] )
        all_peaks = [oracle_pk,] + rep_pks
        new_pk = (oracle_pk.start, oracle_pk.stop, oracle_pk.summit, 
                  pk_agg_fn(signals), 
                  signals, 
                  OrderedDict(zip(range(len(all_peaks)), all_peaks)))
        yield new_pk
    
    return


def merge_peaks_in_contig(all_s_peaks, pk_agg_fn, oracle_pks=None,
                          use_nonoverlapping_peaks=False):
    """Merge peaks in a single contig/strand.
    
    returns: The merged peaks. 
    """
    # merge and sort all peaks, keeping track of which sample they originated in
    oracle_pks_iter = []
    if oracle_pks != None: 
        oracle_pks_iter = oracle_pks
    
    # merge and sort all of the intervals, leeping track of their source
    all_intervals = []
    for sample_id, peaks in enumerate([oracle_pks_iter,] + all_s_peaks):
        all_intervals.extend((pk,sample_id) for pk in peaks)
    all_intervals.sort()
    
    # grp overlapping intervals. Since they're already sorted, all we need
    # to do is check if the current interval overlaps the previous interval
    grpd_intervals = [[],]
    curr_start, curr_stop = all_intervals[0][:2]
    for pk, sample_id in all_intervals:
        if pk.start < curr_stop:
            curr_stop = max(pk.stop, curr_stop)
            grpd_intervals[-1].append((pk, sample_id))
        else:
            curr_start, curr_stop = pk.start, pk.stop
            grpd_intervals.append([(pk, sample_id),])

    # build the unified peak list, setting the score to 
    # zero if it doesn't exist in both replicates
    merged_pks = []
    if oracle_pks == None:
        for intervals in grpd_intervals:
            for merged_pk in iter_merge_grpd_intervals(
                    intervals, len(all_s_peaks), pk_agg_fn,
                    use_oracle_pks=(oracle_pks != None),
                    use_nonoverlapping_peaks = use_nonoverlapping_peaks):
                merged_pks.append(merged_pk)
    else:        
        for intervals in grpd_intervals:
            for merged_pk in iter_matched_oracle_pks(
                    intervals, len(all_s_peaks), pk_agg_fn,
                    use_nonoverlapping_peaks = use_nonoverlapping_peaks):
                merged_pks.append(merged_pk)
    
    return merged_pks

def merge_peaks(all_s_peaks, pk_agg_fn, oracle_pks=None, 
                use_nonoverlapping_peaks=False):
    """Merge peaks over all contig/strands
    
    """
    # if we have reference peaks, use its contigs: otherwise use
    # the union of the replicates contigs
    if oracle_pks != None:
        contigs = sorted(oracle_pks.keys())
    else:
        contigs = sorted(set(chain(*[list(s_peaks.keys()) for s_peaks in all_s_peaks])))

    merged_peaks = []
    for key in contigs:
        # check to see if we've been provided a peak list and, if so, 
        # pass it down. If not, set the oracle peaks to None so that 
        # the callee knows not to use them
        if oracle_pks != None: contig_oracle_pks = oracle_pks[key]
        else: contig_oracle_pks = None
        
        # since s*_peaks are default dicts, it will never raise a key error, 
        # but instead return an empty list which is what we want
        merged_contig_peaks = merge_peaks_in_contig(
            [s_peaks[key] for s_peaks in all_s_peaks], 
            pk_agg_fn, contig_oracle_pks, 
            use_nonoverlapping_peaks=use_nonoverlapping_peaks)
        merged_peaks.extend(
            MergedPeak(*(key + pk)) for pk in merged_contig_peaks)
    
    merged_peaks.sort(key=lambda x:x.merged_signal, reverse=True)
    return merged_peaks

def build_rank_vectors(merged_peaks):
    # allocate memory for the ranks vector
    s1 = numpy.zeros(len(merged_peaks))
    s2 = numpy.zeros(len(merged_peaks))
    # add the signal
    for i, x in enumerate(merged_peaks):
        s1[i], s2[i] = x.signals

    rank1 = numpy.lexsort((numpy.random.random(len(s1)), s1)).argsort()
    rank2 = numpy.lexsort((numpy.random.random(len(s2)), s2)).argsort()
    
    return ( numpy.array(rank1, dtype=numpy.int), 
             numpy.array(rank2, dtype=numpy.int) )

def build_idr_output_line_with_bed6(
        m_pk, IDR, localIDR, output_file_type, signal_type, 
        use_oracle_peak_values=True):
    # initialize the line with the bed6 entires - these are 
    # present in all of the output types
    rv = [m_pk.chrm, str(m_pk.start), str(m_pk.stop), 
          ".", "%i" % (min(1000, int(-125*math.log2(IDR+1e-12)))), m_pk.strand]
    if output_file_type == 'bed':
        # if we just want a bed, there's nothing else to be done
        pass
    # for narrow/broad peak files, we need to add the 3 score fields
    elif output_file_type in ('narrowPeak', 'broadPeak'):
        # if we want to use the oracle peak values for the scores, and an oracle
        # peak is specified
        if use_oracle_peak_values and 0 in m_pk.pks:
            signal_values = [
                m_pk.pks[0].signalValue, m_pk.pks[0].pValue, m_pk.pks[0].qValue]
            signal_values = ["%.5f" % x for x in signal_values]
        else:
            # set the signal values that we didn't use to -1 per the standard 
            signal_values = ["-1", "-1", "-1"]
            signal_values[
                {"signal.value": 0, "p.value": 1, "q.value": 2}[signal_type]
                ] = ("%.5f" % m_pk.merged_signal)
        rv.extend(signal_values)
        # if this is a narrow peak, we also need to add the summit
        if output_file_type == 'narrowPeak':
            rv.append(str(-1 if m_pk.summit == None 
                          else m_pk.summit - m_pk.start))
    else:
        raise ValueError("Unrecognized output format '{}'".format(outputFormat))

    rv.append("%f" % -math.log10(max(1e-5, localIDR)))    
    rv.append("%f" % -math.log10(max(1e-5, IDR)))

    for key, signal in enumerate(m_pk.signals):
        # we add one to the key because key=0 corresponds to the oracle peaks
        key += 1
        # if there is no matching peak for this replicate
        if m_pk.pks[key] == []:
            rv.append("-1")
            rv.append("-1")
            rv.append( "%.5f" % signal )
            if output_file_type == 'narrowPeak':
                rv.append("-1")
        else:
            rv.append( "%i" % min(x.start for x in m_pk.pks[key]))
            rv.append( "%i" % max(x.stop for x in m_pk.pks[key]))
            rv.append( "%.5f" % signal )
            if output_file_type == 'narrowPeak':
                rv.append( "%i" % int(
                    mean(x.summit-x.start for x in m_pk.pks[key] if x.summit != None)
                ))
                                       
            
    return "\t".join(rv)

def build_backwards_compatible_idr_output_line(
        m_pk, IDR, localIDR, output_file_type, signal_type):
    rv = [m_pk.chrm,]
    for key, signal in enumerate(m_pk.signals):
        rv.append( "%i" % min(x.start for x in m_pk.pks[key+1]))
        rv.append( "%i" % max(x.stop for x in m_pk.pks[key+1]))
        rv.append( "." )
        rv.append( "%.5f" % signal )
    
    rv.append("%.5f" % localIDR)
    rv.append("%.5f" % IDR)
    rv.append(m_pk.strand)
        
    return "\t".join(rv)

def calc_local_IDR(theta, r1, r2):
    """
    idr <- 1 - e.z
    o <- order(idr)
    idr.o <- idr[o]
    idr.rank <- rank(idr.o, ties.method = "max")
    top.mean <- function(index, x) {
        mean(x[1:index])
    }
    IDR.o <- sapply(idr.rank, top.mean, idr.o)
    IDR <- idr
    IDR[o] <- IDR.o
    """
    mu, sigma, rho, p = theta
    z1 = compute_pseudo_values(r1, mu, sigma, p, EPS=1e-12)
    z2 = compute_pseudo_values(r2, mu, sigma, p, EPS=1e-12)
    localIDR = 1-calc_post_membership_prbs(numpy.array(theta), z1, z2)
    if idr.FILTER_PEAKS_BELOW_NOISE_MEAN:
        localIDR[z1 + z2 < 0] = 1 

    # it doesn't make sense for the IDR values to be smaller than the 
    # optimization tolerance
    localIDR = numpy.clip(localIDR, idr.CONVERGENCE_EPS_DEFAULT, 1)
    return localIDR

def calc_global_IDR(localIDR):
    local_idr_order = localIDR.argsort()
    ordered_local_idr = localIDR[local_idr_order]
    ordered_local_idr_ranks = rankdata( ordered_local_idr, method='max' )
    IDR = []
    for i, rank in enumerate(ordered_local_idr_ranks):
        IDR.append(ordered_local_idr[:rank].mean())
    IDR = numpy.array(IDR)[local_idr_order.argsort()]
    return IDR

def fit_model_and_calc_local_idr(r1, r2, 
                                 starting_point=None,
                                 max_iter=idr.MAX_ITER_DEFAULT, 
                                 convergence_eps=idr.CONVERGENCE_EPS_DEFAULT, 
                                 fix_mu=False, fix_sigma=False):
    # in theory we would try to find good starting point here,
    # but for now just set it to somethign reasonable
    if starting_point is None:
        starting_point = (idr.DEFAULT_MU, idr.DEFAULT_SIGMA,
                          idr.DEFAULT_RHO, idr.DEFAULT_MIX_PARAM)
    
    idr.log("Initial parameter values: [%s]" % " ".join(
            "%.2f" % x for x in starting_point))
    
    # fit the model parameters    
    idr.log("Fitting the model parameters", 'VERBOSE');
    if idr.PROFILE:
            import cProfile
            cProfile.runctx("""theta, loss = estimate_model_params(
                                    r1,r2,
                                    starting_point, 
                                    max_iter=max_iter, 
                                    convergence_eps=convergence_eps,
                                    fix_mu=fix_mu, fix_sigma=fix_sigma)
                                   """, 
                            {'estimate_model_params': estimate_model_params}, 
                            {'r1':r1, 'r2':r2, 
                             'starting_point': starting_point,
                             'max_iter': max_iter, 
                             'convergence_eps': convergence_eps,
                             'fix_mu': fix_mu, 'fix_sigma': fix_sigma} )
            assert False
    theta, loss = estimate_model_params(
        r1, r2,
        starting_point, 
        max_iter=max_iter, 
        convergence_eps=convergence_eps,
        fix_mu=fix_mu, fix_sigma=fix_sigma)
    
    idr.log("Finished running IDR on the datasets", 'VERBOSE')
    idr.log("Final parameter values: [%s]"%" ".join("%.2f" % x for x in theta))
    
    # calculate the global IDR
    localIDRs = calc_local_IDR(numpy.array(theta), r1, r2)
    return localIDRs

def write_results_to_file(merged_peaks, output_file, 
                          output_file_type, signal_type,
                          max_allowed_idr=1.0,
                          soft_max_allowed_idr=1.0,
                          localIDRs=None, IDRs=None, 
                          useBackwardsCompatibleOutput=False):
    if useBackwardsCompatibleOutput:
        build_idr_output_line = build_backwards_compatible_idr_output_line
    else:
        build_idr_output_line = build_idr_output_line_with_bed6
    
    # write out the result
    idr.log("Writing results to file", "VERBOSE");
    
    if localIDRs is None or IDRs is None:
        assert IDRs is None
        assert localIDRs is None
        localIDRs = numpy.ones(len(merged_peaks))
        IDRs = numpy.ones(len(merged_peaks))

    
    num_peaks_passing_hard_thresh = 0
    num_peaks_passing_soft_thresh = 0
    for localIDR, IDR, merged_peak in zip(
            localIDRs, IDRs, merged_peaks):
        # skip peaks with global idr values below the threshold
        if max_allowed_idr != None and IDR > max_allowed_idr: 
            continue
        num_peaks_passing_hard_thresh += 1
        if IDR <= soft_max_allowed_idr:
            num_peaks_passing_soft_thresh += 1
        opline = build_idr_output_line(
            merged_peak, IDR, localIDR, output_file_type, signal_type)
        print( opline, file=output_file )

    if len(merged_peaks) == 0: return
    
    idr.log(
        "Number of reported peaks - {}/{} ({:.1f}%)\n".format(
            num_peaks_passing_hard_thresh, len(merged_peaks),
            100*float(num_peaks_passing_hard_thresh)/len(merged_peaks))
    )
    
    idr.log(
        "Number of peaks passing IDR cutoff of {} - {}/{} ({:.1f}%)\n".format(
            soft_max_allowed_idr, 
            num_peaks_passing_soft_thresh, len(merged_peaks),
            100*float(num_peaks_passing_soft_thresh)/len(merged_peaks))
    )
    
    return

def parse_args():
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""
Program: IDR (Irreproducible Discovery Rate)
Version: {PACKAGE_VERSION}
Contact: Nathan Boley <npboley@gmail.com>
""".format(PACKAGE_VERSION=idr.__version__))

    def PossiblyGzippedFile(fname):
        if fname.endswith(".gz"):
            return io.TextIOWrapper(gzip.open(fname, 'rb'))
        else:
            return open(fname, 'r')
    
    parser.add_argument( '--samples', '-s', type=PossiblyGzippedFile, nargs=2, 
                         required=True,
                         help='Files containing peaks and scores.')
    parser.add_argument( '--peak-list', '-p', type=PossiblyGzippedFile,
        help='If provided, all peaks will be taken from this file.')
    parser.add_argument( '--input-file-type', default='narrowPeak',
                         choices=['narrowPeak', 'broadPeak', 'bed', 'gff'], 
        help='File type of --samples and --peak-list.')
    
    parser.add_argument( '--rank',
        help="Which column to use to rank peaks."\
            +"\t\nOptions: signal.value p.value q.value columnIndex"\
            +"\nDefaults:\n\tnarrowPeak/broadPeak: signal.value\n\tbed: score")
    
    default_ofname = "idrValues.txt"
    parser.add_argument( '--output-file', "-o", 
                         default=default_ofname, 
        help='File to write output to.\nDefault: {}'.format(default_ofname))
    parser.add_argument( '--output-file-type', 
                         choices=['narrowPeak', 'broadPeak', 'bed'], 
                         default=None, 
        help='Output file type. Defaults to input file type when available, otherwise bed.')

    parser.add_argument( '--log-output-file', "-l", type=argparse.FileType("w"),
                         default=sys.stderr,
                         help='File to write output to. Default: stderr')
    
    parser.add_argument( '--idr-threshold', "-i", type=float, default=None,
        help="Only return peaks with a global idr threshold below this value."\
            +"\nDefault: report all peaks")
    parser.add_argument( '--soft-idr-threshold', type=float, default=None, 
        help="Report statistics for peaks with a global idr below this "\
        +"value but return all peaks with an idr below --idr.\nDefault: %.2f" \
                         % idr.DEFAULT_SOFT_IDR_THRESH)

    parser.add_argument( '--use-old-output-format', 
                         action='store_true', default=False,
                         help="Use old output format.")

    parser.add_argument( '--plot', action='store_true', default=False,
                         help='Plot the results to [OFNAME].png')
        
    parser.add_argument( '--use-nonoverlapping-peaks', 
                         action="store_true", default=False,
        help='Use peaks without an overlapping match and set the value to 0.')
    
    parser.add_argument( '--peak-merge-method', 
                         choices=["sum", "avg", "min", "max"], default=None,
        help="Which method to use for merging peaks.\n" \
              + "\tDefault: 'sum' for signal/score/column indexes, 'min' for p/q-value.")

    parser.add_argument( '--initial-mu', type=float, default=idr.DEFAULT_MU,
        help="Initial value of mu. Default: %.2f" % idr.DEFAULT_MU)
    parser.add_argument( '--initial-sigma', type=float, 
                         default=idr.DEFAULT_SIGMA,
        help="Initial value of sigma. Default: %.2f" % idr.DEFAULT_SIGMA)
    parser.add_argument( '--initial-rho', type=float, default=idr.DEFAULT_RHO,
        help="Initial value of rho. Default: %.2f" % idr.DEFAULT_RHO)
    parser.add_argument( '--initial-mix-param', 
        type=float, default=idr.DEFAULT_MIX_PARAM,
        help="Initial value of the mixture params. Default: %.2f" \
                         % idr.DEFAULT_MIX_PARAM)

    parser.add_argument( '--fix-mu', action='store_true', 
        help="Fix mu to the starting point and do not let it vary.")    
    parser.add_argument( '--fix-sigma', action='store_true', 
        help="Fix sigma to the starting point and do not let it vary.")    

    parser.add_argument( '--dont-filter-peaks-below-noise-mean', 
                         default=False,
                         action='store_true', 
        help="Allow signal points that are below the noise mean (should only be used if you know what you are doing).")    

    parser.add_argument( '--use-best-multisummit-IDR',
                         default=False, action='store_true',
                         help="Set the IDR value for a group of multi summit peaks (a group of peaks with the same chr/start/stop but different summits) to the best value across all of these peaks. This is a work around for peak callers that don't do a good job splitting scores across multi summit peaks (e.g. MACS). If set in conjunction with --plot two plots will be created - one with alternate summits and one without.  Use this option with care.")

    parser.add_argument( '--allow-negative-scores', 
                         default=False,
                         action='store_true', 
        help="Allow negative values for scores. (should only be used if you know what you are doing)")    

    parser.add_argument( '--random-seed', type=int, default=0, 
        help="The random seed value (sor braking ties). Default: 0") 
    parser.add_argument( '--max-iter', type=int, default=idr.MAX_ITER_DEFAULT, 
        help="The maximum number of optimization iterations. Default: %i" 
                         % idr.MAX_ITER_DEFAULT)
    parser.add_argument( '--convergence-eps', type=float, 
                         default=idr.CONVERGENCE_EPS_DEFAULT, 
        help="The maximum change in parameter value changes " \
             + "for convergence. Default: %.2e" % idr.CONVERGENCE_EPS_DEFAULT)
    
    parser.add_argument( '--only-merge-peaks', action='store_true', 
        help="Only return the merged peak list.")    
    
    parser.add_argument( '--verbose', action="store_true", default=False, 
                         help="Print out additional debug information")
    parser.add_argument( '--quiet', action="store_true", default=False, 
                         help="Don't print any status messages")

    parser.add_argument('--version', action='version', 
                        version='IDR %s' % idr.__version__)

    args = parser.parse_args()

    args.output_file = open(args.output_file, "w")
    idr.log_ofp = args.log_output_file

    if args.output_file_type is None:
        if args.input_file_type in ('narrowPeak', 'broadPeak', 'bed'):
            args.output_file_type = args.input_file_type
        else:
            args.output_file_type = 'bed'
    
    if args.verbose: 
        idr.VERBOSE = True 

    global QUIET
    if args.quiet: 
        idr.QUIET = True 
        idr.VERBOSE = False

    if args.dont_filter_peaks_below_noise_mean is True:
        idr.FILTER_PEAKS_BELOW_NOISE_MEAN = False

    if args.allow_negative_scores is True:
        idr.ONLY_ALLOW_NON_NEGATIVE_VALUES = False
        
    assert idr.DEFAULT_IDR_THRESH == 1.0
    if args.idr_threshold == None and args.soft_idr_threshold == None:
        args.idr_threshold = idr.DEFAULT_IDR_THRESH
        args.soft_idr_threshold = idr.DEFAULT_SOFT_IDR_THRESH
    elif args.soft_idr_threshold == None:
        assert args.idr_threshold != None
        args.soft_idr_threshold = args.idr_threshold
    elif args.idr_threshold == None:
        assert args.soft_idr_threshold != None
        args.idr_threshold = idr.DEFAULT_IDR_THRESH

    numpy.random.seed(args.random_seed)

    if args.plot:
        try: 
            import matplotlib
        except ImportError:
            idr.log("WARNING: matplotlib does not appear to be installed and "\
                    +"is required for plotting - turning plotting off.", 
                    level="WARNING" )
            args.plot = False
    
    return args

def load_samples(args):
    # decide what aggregation function to use for peaks that need to be merged
    idr.log("Loading the peak files", 'VERBOSE')
    if args.input_file_type in ['narrowPeak', 'broadPeak']:
        if args.rank == None: signal_type = 'signal.value'
        else: signal_type = args.rank

        try: 
            signal_index = {"score": 4, "signal.value": 6, 
                            "p.value": 7, "q.value": 8}[signal_type]
        except KeyError:
            raise ValueError(
                "Unrecognized signal type for {} filetype: '{}'".format(
                    args.input_file_type, signal_type))

        if args.peak_merge_method != None:
            peak_merge_fn = {
                "sum": sum, "avg": mean, "min": min, "max": max}[
                    args.peak_merge_method]
        elif signal_index in (4,6):
            peak_merge_fn = sum
        else:
            peak_merge_fn = min
        if args.input_file_type == 'narrowPeak':
            summit_index = 9
        else:
            summit_index = None
        f1, f2 = [load_bed(fp, signal_index, summit_index) 
                  for fp in args.samples]
        oracle_pks =  (
            load_bed(args.peak_list, signal_index, summit_index) 
            if args.peak_list != None else None)
    elif args.input_file_type in ['bed', ]:
        # set the default
        if args.rank == None: 
            signal_type = 'score'

        if args.rank == 'score':
            signal_type = 'score'
            signal_index = 4
        else:
            try: 
                signal_index = int(args.rank) - 1
                signal_type = "COL%i" % (signal_index + 1)
            except ValueError:
                raise ValueError("For bed files --signal-type must either "\
                                 +"be set to score or an index specifying "\
                                 +"the column to use.")
        
        if args.peak_merge_method == None:
            peak_merge_fn = sum
        else:
            peak_merge_fn = {
                "sum": sum, "avg": mean, "min": min, "max": max}[
                    args.peak_merge_method]
        
        f1, f2 = [load_bed(fp, signal_index) for fp in args.samples]
        oracle_pks =  (
            load_bed(args.peak_list, signal_index) 
            if args.peak_list != None else None)
    elif args.input_file_type in ['gff', ]:
        # set the default
        if args.rank == None: 
            signal_type = 'score'
        else:
            assert args.rank == 'score'
        
        if args.peak_merge_method == None:
            peak_merge_fn = sum
        else:
            peak_merge_fn = {
                "sum": sum, "avg": mean, "min": min, "max": max}[
                    args.peak_merge_method]
        
        f1, f2 = [load_gff(fp) for fp in args.samples]
        oracle_pks =  (
            load_gff(args.peak_list) 
            if args.peak_list != None else None)
    else:
        raise ValueError( "Unrecognized file type: '{}'".format(
            args.input_file_type))
    # build a unified peak set
    idr.log("Merging peaks", 'VERBOSE')
    merged_peaks = merge_peaks([f1, f2], peak_merge_fn, 
                               oracle_pks, args.use_nonoverlapping_peaks)
    return merged_peaks, signal_type

def plot(args, scores, ranks, IDRs, ofprefix=None):
    assert len(args.samples) == 2
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot

    if ofprefix is None:
        ofprefix = args.output_file.name
    
    colors = []
    for i in range(len(IDRs)):
        colors.append('#FC9272' if IDRs[i]>args.soft_idr_threshold else 'k')
        
    #matplotlib.rc('font', family='normal', weight='bold', size=10)

    fig = matplotlib.pyplot.figure( num=None, figsize=(12, 12))

    matplotlib.pyplot.subplot(221)
    matplotlib.pyplot.axis([0, 1, 0, 1])
    matplotlib.pyplot.xlabel("Sample 1 Rank")
    matplotlib.pyplot.ylabel("Sample 2 Rank")
    matplotlib.pyplot.title(
        "Ranks - (red >= %.2f IDR)" % args.soft_idr_threshold)
    matplotlib.pyplot.scatter((ranks[0]+1)/float(max(ranks[0])+1), 
                              (ranks[1]+1)/float(max(ranks[1])+1), 
                              edgecolor='none',
                              c=colors,
                              alpha=0.25,s=5)

    matplotlib.pyplot.subplot(222)
    matplotlib.pyplot.xlabel("Sample 1 log10 Score")
    matplotlib.pyplot.ylabel("Sample 2 log10 Score")
    matplotlib.pyplot.title(
        "Log10 Scores - (red >= %.2f IDR)" % args.soft_idr_threshold)
    matplotlib.pyplot.scatter(numpy.log10(scores[0]+1),
                              numpy.log10(scores[1]+1), 
                              edgecolor='none',
                              c=colors, 
                              alpha=0.25,s=5)
    
    def make_boxplots(sample_i):
        groups = defaultdict(list)
        norm_ranks = (ranks[sample_i]+1)/float(max(ranks[sample_i])+1)
        for rank, idr_val in zip(norm_ranks, -numpy.log10(IDRs)):
            groups[int(20*rank)].append(float(idr_val))
        group_labels = sorted((x + 2.5)/20 for x in groups.keys())
        groups = [x[1] for x in sorted(groups.items())]

        matplotlib.pyplot.title(
            "Sample %i Ranks vs IDR Values" % (sample_i+1))
        matplotlib.pyplot.axis([0, 21, 
                                0, 0.5-math.log10(idr.CONVERGENCE_EPS_DEFAULT)])
        matplotlib.pyplot.xlabel("Sample %i Peak Rank" % (sample_i+1))
        matplotlib.pyplot.ylabel("-log10 IDR")
        matplotlib.pyplot.xticks([], [])

        matplotlib.pyplot.boxplot(groups, sym="")    

        matplotlib.pyplot.axis([0, 21, 
                                0, 0.5-math.log10(idr.CONVERGENCE_EPS_DEFAULT)])
        matplotlib.pyplot.scatter(20*norm_ranks, 
                                  -numpy.log10(IDRs), 
                                  alpha=0.1, c='black', edgecolor='none', s=5)

    matplotlib.pyplot.subplot(223)
    make_boxplots(0)

    matplotlib.pyplot.subplot(224)
    make_boxplots(1)

    matplotlib.pyplot.savefig(ofprefix + ".png")
    return

def main():
    args = parse_args()

    # load and merge peaks
    merged_peaks, signal_type = load_samples(args)
    s1 = numpy.array([pk.signals[0] for pk in merged_peaks])
    s2 = numpy.array([pk.signals[1] for pk in merged_peaks])

    # build the ranks vector
    idr.log("Ranking peaks", 'VERBOSE')
    r1, r2 = build_rank_vectors(merged_peaks)
    
    if args.only_merge_peaks:
        localIDRs, IDRs = None, None
    else:
        if len(merged_peaks) < 20:
            error_msg = "Peak files must contain at least 20 peaks post-merge"
            error_msg += "\nHint: Merged peaks were written to the output file"
            write_results_to_file(
                merged_peaks, args.output_file,
                args.output_file_type, signal_type)
            raise ValueError(error_msg)

        localIDRs = fit_model_and_calc_local_idr(
            r1, r2, 
            starting_point=(
                args.initial_mu, args.initial_sigma, 
                args.initial_rho, args.initial_mix_param),
            max_iter=args.max_iter,
            convergence_eps=args.convergence_eps,
            fix_mu=args.fix_mu, fix_sigma=args.fix_sigma )    

        # if the use chose to use the best multi summit IDR, then
        # make the correction and plot just the corrected peaks
        if args.use_best_multisummit_IDR:
            update_indices, localIDRs = correct_multi_summit_peak_IDR_values(
                localIDRs, merged_peaks)
            IDRs = calc_global_IDR(localIDRs)        
            if args.plot:
                assert len(args.samples) == 2
                plot(args,
                     [s1[update_indices], s2[update_indices]],
                     [r1[update_indices], r2[update_indices]],
                     IDRs[update_indices],
                     args.output_file.name + ".noalternatesummitpeaks")
        # we wrap this in an else statement to avoid calculating the global IDRs
        # twice
        else:
            IDRs = calc_global_IDR(localIDRs)        
        
        if args.plot:
            assert len(args.samples) == 2
            plot(args, [s1, s2], [r1, r2], IDRs)

    
    num_peaks_passing_thresh = write_results_to_file(
        merged_peaks, 
        args.output_file, 
        args.output_file_type, 
        signal_type,
        localIDRs=localIDRs, 
        IDRs=IDRs,
        max_allowed_idr=args.idr_threshold,
        soft_max_allowed_idr=args.soft_idr_threshold,
        useBackwardsCompatibleOutput=args.use_old_output_format)
    
    args.output_file.close()

if __name__ == '__main__':
    try:
        main()
    finally:
        if idr.log_ofp != sys.stderr: idr.log_ofp.close()
